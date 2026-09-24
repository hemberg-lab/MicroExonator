import csv
import gzip
import importlib.util
import os
import pathlib
import tempfile
import unittest
from unittest import mock


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
FIXTURE_ROOT = REPOSITORY / "tests" / "integration" / "mm10_simulated"
VALIDATOR_PATH = FIXTURE_ROOT / "validate_smoke.py"
RUNNER_PATH = FIXTURE_ROOT / "run_smoke_test.py"


def load_validator():
    spec = importlib.util.spec_from_file_location("mm10_smoke_validator", VALIDATOR_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_runner():
    spec = importlib.util.spec_from_file_location("mm10_smoke_runner", RUNNER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_fastq(path, records):
    with gzip.open(path, "wt") as handle:
        for read_id, sequence, qualities in records:
            handle.write("@{}\n{}\n+\n{}\n".format(read_id, sequence, qualities))


class FastqValidationTests(unittest.TestCase):
    def setUp(self):
        self.validator = load_validator()

    def test_fastq_rejects_duplicate_ids_and_wrong_read_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            fastq = pathlib.Path(directory) / "reads.fastq.gz"
            write_fastq(
                fastq,
                [
                    ("read-1", "A" * 100, "I" * 100),
                    ("read-1", "C" * 99, "I" * 99),
                ],
            )

            with self.assertRaisesRegex(
                self.validator.ValidationError, "duplicate read ID.*100 nt"
            ):
                self.validator.inspect_fastq(fastq)

    def test_paired_fastq_rejects_unsynchronized_mates(self):
        with tempfile.TemporaryDirectory() as directory:
            r1 = pathlib.Path(directory) / "r1.fastq.gz"
            r2 = pathlib.Path(directory) / "r2.fastq.gz"
            write_fastq(r1, [("pair-1/1", "A" * 100, "I" * 100)])
            write_fastq(r2, [("different/2", "T" * 100, "I" * 100)])

            with self.assertRaisesRegex(
                self.validator.ValidationError, "not synchronized"
            ):
                self.validator.inspect_paired_fastqs(r1, r2)


class TruthValidationTests(unittest.TestCase):
    def setUp(self):
        self.validator = load_validator()

    def test_truth_rejects_event_sample_without_twelve_valid_isoform_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            truth = pathlib.Path(directory) / "reads.tsv.gz"
            fields = [
                "read_id",
                "sample",
                "kind",
                "event_id",
                "ME",
                "isoform",
                "anchor_up",
                "anchor_down",
                "valid_evidence",
                "error_positions_r1",
                "error_positions_r2",
            ]
            with gzip.open(truth, "wt", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for index in range(12):
                    writer.writerow(
                        {
                            "read_id": "included-{}".format(index),
                            "sample": "A1",
                            "kind": "target",
                            "event_id": "event-1",
                            "ME": "chr1_+_100_103",
                            "isoform": "included",
                            "anchor_up": 20,
                            "anchor_down": 77,
                            "valid_evidence": "1",
                            "error_positions_r1": "",
                            "error_positions_r2": "",
                        }
                    )
                for index in range(11):
                    writer.writerow(
                        {
                            "read_id": "skipped-{}".format(index),
                            "sample": "A1",
                            "kind": "target",
                            "event_id": "event-1",
                            "ME": "chr1_+_100_103",
                            "isoform": "skipped",
                            "anchor_up": 50,
                            "anchor_down": 50,
                            "valid_evidence": "1",
                            "error_positions_r1": "",
                            "error_positions_r2": "",
                        }
                    )

            with self.assertRaisesRegex(
                self.validator.ValidationError,
                "event-1.*A1.*skipped.*11.*minimum is 12",
            ):
                self.validator.inspect_truth_support(
                    truth, {"event-1"}, {"A1"}, minimum_valid_support=12
                )


class CommittedFixtureIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.validator = load_validator()

    def test_single_end_fixture_satisfies_committed_truth_contract(self):
        report = self.validator.validate_fixture(FIXTURE_ROOT, "single_end")

        self.assertEqual(report["events"], 100)
        self.assertEqual(report["samples"], 9)
        self.assertEqual(report["fastq_reads"], report["truth_records"])
        self.assertGreater(report["distinct_anchor_pairs"], 100)

    def test_paired_end_fixture_satisfies_committed_truth_contract(self):
        report = self.validator.validate_fixture(FIXTURE_ROOT, "paired_end")

        self.assertEqual(report["events"], 100)
        self.assertEqual(report["samples"], 9)
        self.assertEqual(report["fastq_pairs"], report["truth_records"])
        self.assertGreater(report["distinct_anchor_pairs"], 100)


class PipelineOutputValidationTests(unittest.TestCase):
    def setUp(self):
        self.validator = load_validator()

    def make_run(self, root):
        truth_dir = root / "truth"
        truth_dir.mkdir()
        with open(truth_dir / "events.tsv", "w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=("event_id", "ME", "profile"), delimiter="\t"
            )
            writer.writeheader()
            writer.writerow(
                {"event_id": "event-1", "ME": "chr1_+_100_103", "profile": "stable"}
            )
            writer.writerow(
                {"event_id": "event-2", "ME": "chr2_-_200_214", "profile": "stable"}
            )

        layout_dir = root / "single_end"
        layout_dir.mkdir()
        with open(layout_dir / "bulk_samples.tsv", "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("sample", "condition"))
            writer.writerow(("A1", "A"))
            writer.writerow(("A2", "A"))

        run_dir = root / "run"
        quant_dir = run_dir / "Report" / "quant" / "corrected" / "PSI_sparse" / "bulk" / "se"
        quant_dir.mkdir(parents=True)
        for sample in ("A1", "A2"):
            with gzip.open(quant_dir / "{}.corrected.PSI.gz".format(sample), "wt", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(("sample", "ME", "ME_coverages", "excluding_covs", "PSI", "CI_Lo", "CI_Hi"))
                writer.writerow((sample, "chr1_+_100_103", 30, 30, 0.5, 0.36, 0.64))
                writer.writerow((sample, "chr2_-_200_214", 45, 15, 0.75, 0.62, 0.85))

        robust = run_dir / "Report" / "out.robustly_detected.txt"
        robust.parent.mkdir(exist_ok=True)
        with open(robust, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("ME", "Transcript"))
            writer.writerow(("chr1_+_100_103", "tx1"))
            writer.writerow(("chr2_-_200_214", "tx2"))
            writer.writerow(("chr3_+_300_306", "extra"))
        return truth_dir / "events.tsv", layout_dir / "bulk_samples.tsv", run_dir

    def test_completed_run_accepts_all_truth_events_and_reports_extras(self):
        with tempfile.TemporaryDirectory() as directory:
            events, samples, run_dir = self.make_run(pathlib.Path(directory))

            report = self.validator.validate_pipeline_outputs(
                "single_end", events, samples, run_dir
            )

            self.assertEqual(report["truth_events"], 2)
            self.assertEqual(report["quantified_event_samples"], 4)
            self.assertEqual(report["additional_discoveries"], ["chr3_+_300_306"])

    def test_completed_run_allows_only_declared_collapsed_truth_event(self):
        """A known repeated 3-nt sequence is not required as a separate ME."""
        with tempfile.TemporaryDirectory() as directory:
            events, samples, run_dir = self.make_run(pathlib.Path(directory))
            collapsed = "chr2_-_200_214"

            robust = run_dir / "Report" / "out.robustly_detected.txt"
            with open(robust, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(("ME", "Transcript"))
                writer.writerow(("chr1_+_100_103", "tx1"))

            quant_dir = run_dir / "Report" / "quant" / "corrected" / "PSI_sparse" / "bulk" / "se"
            for quant in quant_dir.glob("*.corrected.PSI.gz"):
                with gzip.open(quant, "rt") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                with gzip.open(quant, "wt", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writer.writerow(rows[0])

            report = self.validator.validate_pipeline_outputs(
                "single_end", events, samples, run_dir, expected_missing={collapsed}
            )

            self.assertEqual(report["truth_events"], 2)
            self.assertEqual(report["expected_missing_events"], [collapsed])
            self.assertEqual(report["quantified_event_samples"], 2)

    def test_completed_run_rejects_nonnumeric_confidence_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            events, samples, run_dir = self.make_run(pathlib.Path(directory))
            quant = (
                run_dir
                / "Report"
                / "quant"
                / "corrected"
                / "PSI_sparse"
                / "bulk"
                / "se"
                / "A2.corrected.PSI.gz"
            )
            with gzip.open(quant, "wt", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(("sample", "ME", "ME_coverages", "excluding_covs", "PSI", "CI_Lo", "CI_Hi"))
                writer.writerow(("A2", "chr1_+_100_103", 30, 30, 0.5, "NA", 0.64))
                writer.writerow(("A2", "chr2_-_200_214", 45, 15, 0.75, 0.62, 0.85))

            with self.assertRaisesRegex(
                self.validator.ValidationError,
                "A2.*chr1_\\+_100_103.*CI_Lo.*numeric",
            ):
                self.validator.validate_pipeline_outputs(
                    "single_end", events, samples, run_dir
                )


    def write_targets(self, events, target_psi):
        """Rewrite events.tsv with a psi_A target column."""
        with open(events, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=("event_id", "ME", "profile", "psi_A"), delimiter="\t"
            )
            writer.writeheader()
            for index, (microexon, psi) in enumerate(target_psi.items(), start=1):
                writer.writerow(
                    {"event_id": "event-{}".format(index), "ME": microexon, "profile": "stable", "psi_A": psi}
                )

    def test_completed_run_accepts_psi_close_to_simulated_targets(self):
        with tempfile.TemporaryDirectory() as directory:
            events, samples, run_dir = self.make_run(pathlib.Path(directory))
            self.write_targets(events, {"chr1_+_100_103": 0.45, "chr2_-_200_214": 0.75})

            report = self.validator.validate_pipeline_outputs(
                "single_end", events, samples, run_dir
            )

            deviations = {item["ME"]: item["deviation"] for item in report["psi_deviations"]}
            self.assertAlmostEqual(deviations["chr1_+_100_103"], 0.05)
            self.assertAlmostEqual(deviations["chr2_-_200_214"], 0.0)

    def test_completed_run_rejects_psi_far_from_simulated_target(self):
        """An event pinned at high PSI despite simulated skipping must fail."""
        with tempfile.TemporaryDirectory() as directory:
            events, samples, run_dir = self.make_run(pathlib.Path(directory))
            self.write_targets(events, {"chr1_+_100_103": 0.5, "chr2_-_200_214": 0.35})

            with self.assertRaisesRegex(
                self.validator.ValidationError,
                "1 truth events deviate.*chr2_-_200_214 A observed 0.75 target 0.35",
            ):
                self.validator.validate_pipeline_outputs(
                    "single_end", events, samples, run_dir
                )

    def test_psi_accuracy_check_can_be_disabled(self):
        with tempfile.TemporaryDirectory() as directory:
            events, samples, run_dir = self.make_run(pathlib.Path(directory))
            self.write_targets(events, {"chr1_+_100_103": 0.5, "chr2_-_200_214": 0.35})

            report = self.validator.validate_pipeline_outputs(
                "single_end", events, samples, run_dir, psi_tolerance=None
            )

            self.assertEqual(len(report["psi_deviations"]), 2)


class SmokeRunnerTests(unittest.TestCase):
    def test_execution_environment_puts_conda_base_before_system_python(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory) / "miniconda3"
            snakemake = root / "envs" / "snakemake" / "bin" / "snakemake"
            snakemake.parent.mkdir(parents=True)
            snakemake.write_text("")
            with mock.patch.dict(
                runner.os.environ,
                {"PATH": "/usr/bin:{}".format(root / "bin")},
                clear=True,
            ):
                environment = runner._execution_environment(snakemake)

            self.assertEqual(
                pathlib.Path(environment["PATH"].split(os.pathsep)[0]),
                (root / "bin").resolve(),
            )

    def test_samtools_fallback_finds_sibling_conda_environment(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            envs = pathlib.Path(directory) / "miniconda3" / "envs"
            snakemake = envs / "snakemake" / "bin" / "snakemake"
            samtools = envs / "biotools" / "bin" / "samtools"
            snakemake.parent.mkdir(parents=True)
            samtools.parent.mkdir(parents=True)
            snakemake.write_text("")
            samtools.write_text("")

            with mock.patch.dict(runner.os.environ, {}, clear=True), mock.patch.object(
                runner.shutil, "which", return_value=None
            ):
                observed = runner._resolve_samtools(snakemake)

            self.assertEqual(pathlib.Path(observed), samtools.resolve())

    def test_prepare_workdir_stages_absolute_inputs_and_robustness_config(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            fixture = root / "fixture"
            layout = fixture / "single_end"
            fastq = layout / "fastq"
            fastq.mkdir(parents=True)
            write_fastq(fastq / "A1.fastq.gz", [("read-1", "A" * 100, "I" * 100)])
            (fixture / "empty_microexons.bed").write_text("")
            (layout / "local_samples.tsv").write_text("sample\tpath\nA1\tfastq/A1.fastq.gz\n")
            (layout / "bulk_samples.tsv").write_text("sample\tcondition\nA1\tA\n")
            annotation = root / "annotation.bed.gz"
            with gzip.open(annotation, "wt") as handle:
                handle.write("chr1\t0\t200\ttx1\t0\t+\t0\t200\t0\t1\t200,\t0,\n")
            genome = root / "mm10.fa"
            genome.write_text(">chr1\n" + "A" * 200 + "\n")
            repository = root / "repository"
            (repository / "PWM" / "Mouse").mkdir(parents=True)
            (repository / "PWM" / "Mouse" / "mm10_GT_AG_U2_5.good.matrix").write_text("PWM5\n")
            (repository / "PWM" / "Mouse" / "mm10_GT_AG_U2_3.good.matrix").write_text("PWM3\n")
            (repository / "src").mkdir()
            workdir = root / "run"

            config_path = runner.prepare_workdir(
                "single_end", genome, annotation, workdir, fixture, repository
            )

            config = runner.read_simple_yaml(config_path)
            self.assertEqual(config["filter_method"], "robustness")
            self.assertEqual(config["min_detected_samples"], "2")
            self.assertEqual(config["min_reads_PSI"], "5")
            self.assertEqual(config["conservation_bigwig"], "NA")
            self.assertEqual(config["Gene_anontation_GTF"], "NA")
            self.assertEqual(pathlib.Path(config["Genome_fasta"]), genome.resolve())
            self.assertTrue(pathlib.Path(config["Gene_anontation_bed12"]).is_file())
            self.assertTrue((workdir / "data" / "Genome").is_symlink())
            self.assertTrue((workdir / "data" / "Genome").samefile(genome))
            self.assertTrue((workdir / "FASTQ" / "A1.fastq.gz").is_symlink())
            self.assertTrue((workdir / "FASTQ" / "A1.fastq.gz").samefile(fastq / "A1.fastq.gz"))
            download_script = workdir / "download" / "A1.download.sh"
            self.assertTrue(download_script.is_file())
            self.assertLessEqual(
                download_script.stat().st_mtime_ns,
                (workdir / "FASTQ" / "A1.fastq.gz").lstat().st_mtime_ns,
            )
            with open(workdir / "local_samples.tsv") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(pathlib.Path(rows[0]["path"]), (fastq / "A1.fastq.gz").resolve())


if __name__ == "__main__":
    unittest.main()
