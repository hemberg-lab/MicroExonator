"""A real Snakemake run of the circRNA modules on the synthetic fixture (cached_only).

Kept caches are made first with the real tools (CIRIquant counting on a HISAT2 alignment;
STAR-format chimeric junctions emulated, since STAR's macOS builds cannot read gzipped
input). Snakemake then executes, not dry-runs, the catalogue build with an hg19 source
lifted to hg38, CIRCexplorer2 calls from the kept STAR caches, the CIRCexplorer2 and DCC
comparisons, the CIRIquant comparison and the module inventory, with the tools taken from
the two environments on PATH (CIRCRNA_DISCOVERY_ENV first, then CIRCRNA_CIRIQUANT_ENV).
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import circrna_fixture as fixture

REPOSITORY = Path(__file__).resolve().parent.parent
CIRIQUANT_ENV = os.environ.get("CIRCRNA_CIRIQUANT_ENV")
DISCOVERY_ENV = os.environ.get("CIRCRNA_DISCOVERY_ENV")
SNAKEMAKE = os.environ.get("MICROEXONATOR_SNAKEMAKE")
# optional: an environment from envs/umbrella-circrna-diff.yaml (numpy, R, edgeR) runs the
# differential rule too, on junction shards built here by the pipeline's own extractor
DIFF_ENV = os.environ.get("CIRCRNA_DIFF_ENV")
RUNS = (("run_a1", "ctrl"), ("run_a2", "ctrl"), ("run_b1", "case"), ("run_b2", "case"), ("run_c1", "other"))


def read_tsv(path):
    with gzip.open(path, "rt") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in stream]


@unittest.skipUnless(CIRIQUANT_ENV and DISCOVERY_ENV and SNAKEMAKE,
                     "set CIRCRNA_CIRIQUANT_ENV, CIRCRNA_DISCOVERY_ENV and MICROEXONATOR_SNAKEMAKE")
class CircrnaWorkflowNativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        work = cls.work = Path(cls.tmp.name)
        for name in ("rules", "src", "envs"):
            shutil.copytree(REPOSITORY / name, work / name, ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copy2(REPOSITORY / "MicroExonator.smk", work)
        fx = cls.fx = fixture.build(work / "fx")
        cls.env = dict(os.environ, PATH=os.pathsep.join([str(Path(DISCOVERY_ENV) / "bin"),
                                                         str(Path(CIRIQUANT_ENV) / "bin")] +
                                                        ([str(Path(DIFF_ENV) / "bin")] if DIFF_ENV else []) +
                                                        [os.environ["PATH"]]))
        # catalogue sources: the fixture's hg38 BED, and a "hg19" copy lifted through an identity chain
        shutil.copy(fx["catalogue"], work / "circatlas.bed")
        (work / "circbase_hg19.bed").write_text("chr2\t2200\t2400\tcircBase_C2\t0\t-\nchr1\t5\t50\tunlifted\t0\t+\n")
        with gzip.open(work / "hg19ToHg38.over.chain.gz", "wt") as stream:
            stream.write("chain 1 chr2 4000 + 1000 4000 chr2 4000 + 1000 4000 1\n3000\n")
        for name in ("u2_5.matrix", "u2_3.matrix", "conservation.bw", "microexons.bed12", "annotation.bed12",
                     "transcripts.fa", "decoys.txt", "fixed.gtf", "salmon.gtf", "splice_sites.txt"):
            (work / name).write_text("")
        shutil.copy(fx["genome"], work / "genome.fa")
        shutil.copy(fx["gtf"], work / "annotation.gtf")
        with open(work / "umbrella.tsv", "w") as stream:
            stream.write("sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\tsource_type\t"
                         "source_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\n")
            for run, group in RUNS:
                stream.write("s_{0}\t{0}\trep_{0}\tproject\tbatch\t{1}\tfastq\tfx/R1.fastq.gz\tfx/R2.fastq.gz\tPE\t"
                             "unstranded\tauto\ttrue\n".format(run, group))
        (work / "comparisons.yaml").write_text(
            "comparisons:\n  - comparison_id: pilot\n    project_id: project\n    group_a: case\n    group_b: ctrl\n"
            "  - comparison_id: unselected\n    project_id: project\n    group_a: other\n    group_b: ctrl\n")
        (work / "selection.tsv").write_text(
            "comparison_id\ttool\trole\tlibrary_evidence\tdecided_by\tapproved_by\treason\n" + "".join(
                "pilot\t{}\ttest\tRibo-Zero both sides (fixture)\tguillermo\tguillermo\tE2E\n".format(tool)
                for tool in ("ciriquant", "circexplorer2", "dcc")))
        (work / "config.yaml").write_text("\n".join([
            "Genome_fasta: genome.fa", "Gene_anontation_bed12: annotation.bed12", "Gene_anontation_GTF: annotation.gtf",
            "GT_AG_U2_5: u2_5.matrix", "GT_AG_U2_3: u2_3.matrix", "conservation_bigwig: conservation.bw",
            "ME_DB: microexons.bed12", "working_directory: {}/".format(work), "min_number_files_detected: 1",
            "Single_Cell: F", "umbrella_manifest: umbrella.tsv", "umbrella_comparisons: comparisons.yaml",
            "whippet_bin_folder: /opt/whippet/bin", "julia: julia", "umbrella_reference:",
            "  genome_fasta: genome.fa", "  annotation_gtf: annotation.gtf", "  whippet_gtf: fixed.gtf",
            "  me_db: microexons.bed12", "  transcriptome_fasta: transcripts.fa", "  decoys: decoys.txt",
            "  salmon_gtf: salmon.gtf", "  splice_sites: splice_sites.txt",
            "umbrella_modules: [ciriquant, circexplorer2, dcc]", "umbrella_modules_mode: cached_only",
            "umbrella_module_selection: selection.tsv",
            "circrna_catalogue_sources:", "  - {path: circatlas.bed, name: circAtlas}",
            "  - {path: circbase_hg19.bed, assembly: hg19, role: annotate, name: circBase}",
            "circrna_liftover_chain: hg19ToHg38.over.chain.gz", "circrna_star_sa_index_nbases: 5"] +
            (["circrna_differential: true", "circrna_hisat_delta: true", "circrna_delta_draws: 2000"]
             if DIFF_ENV else [])) + "\n")
        # the reference ID the workflow computes from umbrella_reference (reference_id: auto)
        probe = subprocess.run([SNAKEMAKE, "-s", "MicroExonator.smk", "-n", "quant_umbrella_modules", "--configfile",
                                "config.yaml", "--nocolor"], cwd=work, env=cls.env, capture_output=True, text=True)
        import re
        cls.ref = re.search(r"Umbrella reference_id: (\w+)", probe.stdout + probe.stderr).group(1)
        cls.write_caches()
        cls.result = subprocess.run(
            [SNAKEMAKE, "-s", "MicroExonator.smk", "quant_umbrella_modules", "--cores", "2", "--configfile",
             "config.yaml", "--nocolor"], cwd=work, env=cls.env, capture_output=True, text=True)

    @classmethod
    def write_caches(cls):
        """Kept caches for the pilot's runs, as an ingest leaves them."""
        work, fx = cls.work, cls.fx
        python = str(Path(CIRIQUANT_ENV) / "bin" / "python")
        # CIRIquant on a HISAT2 alignment of the fixture reads, against this workflow's catalogue
        catalogue = work / "prep_catalogue"
        subprocess.run([python, "src/circrna_catalogue.py", "build", "--input", str(fx["catalogue"]),
                        "--genome", str(fx["genome"]), "--out-dir", str(catalogue)], cwd=work, env=cls.env,
                       check=True, capture_output=True)
        subprocess.run("hisat2-build -q fx/genome.fa hisat_fx && hisat2 -x hisat_fx -1 fx/R1.fastq.gz -2 fx/R2.fastq.gz "
                       "2> /dev/null | samtools sort -o fx.bam - && samtools index fx.bam", shell=True, cwd=work,
                       env=cls.env, check=True)
        subprocess.run([python, "src/run_ciriquant.py", "prep", "--reads", "fx/R1.fastq.gz", "fx/R2.fastq.gz",
                        "--bam", "fx.bam", "--catalogue", str(catalogue / "catalogue.bed"), "--genome", "fx/genome.fa",
                        "--out-dir", "prep_ciriquant", "--run-id", "x", "--threads", "2", "--tmp-root", "tmp"],
                       cwd=work, env=cls.env, check=True, capture_output=True)
        if DIFF_ENV:
            cls.write_shards()
        text, cls.truth = fixture.star_chimeric_lines()
        sys.path.insert(0, str(work))
        cwd = os.getcwd()
        os.chdir(work)
        try:
            import yaml
            from src.umbrella_manifest import load_umbrella_manifest
            from src.umbrella_modules import cache_id, module_reference_id
            config = yaml.safe_load((work / "config.yaml").read_text())
            references = {tool: module_reference_id(cls.ref, tool, config) for tool in ("ciriquant", "circrna_star")}
            for record in load_umbrella_manifest(work / "umbrella.tsv").with_reference_id(cls.ref).included_runs():
                if record.group not in ("ctrl", "case"):
                    continue
                for tool, reference in references.items():
                    folder = work / "umbrella/modules/cache" / cls.ref / "project/batch" / record.run_id / tool / cache_id(
                        record, reference)
                    folder.mkdir(parents=True)
                    if tool == "ciriquant":
                        shutil.copy(work / "prep_ciriquant" / "counts.tsv.gz", folder)
                        record_json = json.loads((work / "prep_ciriquant" / "cache.json").read_text())
                        record_json.update(run_id=record.run_id, module_reference_id=reference)
                        (folder / "cache.json").write_text(json.dumps(record_json))
                    else:
                        with gzip.open(folder / "Chimeric.out.junction.gz", "wt") as stream:
                            stream.write(text)
                        with gzip.open(folder / "SJ.out.tab.gz", "wt") as stream:
                            stream.write("")
                        (folder / "Log.final.out").write_text("")
                        (folder / "cache.json").write_text(json.dumps({"run_id": record.run_id}))
        finally:
            os.chdir(cwd)
            sys.path.remove(str(work))

    @classmethod
    def write_shards(cls):
        """Retained junction shards as the core pipeline keeps them, one per group."""
        work = cls.work
        subprocess.run("samtools view fx.bam | python3 src/summarize_junctions.py extract --table fx.junctions.tsv.gz "
                       "--summary fx.summary.json", shell=True, cwd=work, env=cls.env, check=True)
        from tests import circrna_fixture as fixture
        hints = []
        for chrom, strand, exons in (("chr1", "+", fixture.G1_EXONS), ("chr2", "-", fixture.G2_EXONS)):
            ordered = sorted(exons)
            hints += ["{}\t{}\t{}\t{}\n".format(chrom, a - 1, b - 1, strand)
                      for (_, a), (b, _) in zip(ordered, ordered[1:])]
        (work / "hints.txt").write_text("".join(hints))
        subprocess.run(["python3", "src/summarize_junctions.py", "catalog", "--whippet-gtf", "fx/annotation.gtf",
                        "--annotation-gtf", "fx/annotation.gtf", "--splice-sites", "hints.txt",
                        "--output", "catalog.tsv.gz"], cwd=work, env=cls.env, check=True)
        for group, runs in (("ctrl", ("run_a1", "run_a2")), ("case", ("run_b1", "run_b2"))):
            folder = work / "junctions" / cls.ref / "project" / group
            folder.mkdir(parents=True, exist_ok=True)
            subprocess.run(["python3", "src/summarize_junctions.py", "shard", "--catalog", "catalog.tsv.gz",
                            "--prefix", str(folder / "batch"), "--reference-id", cls.ref, "--manifest-sha256", "x" * 64]
                           + ["--table={}=fx.junctions.tsv.gz".format(run) for run in runs],
                           cwd=work, env=cls.env, check=True)

    @unittest.skipUnless(DIFF_ENV, "set CIRCRNA_DIFF_ENV to run the differential rule as well")
    def test_differential_rule_runs_from_kept_counts_and_shards(self):
        for tool in ("ciriquant", "circexplorer2", "dcc"):
            folders = list(self.analysis(tool).glob("differential/*"))
            self.assertEqual(len(folders), 1, tool)
            status = dict(line.split("\t", 1) for line in (folders[0] / "status.tsv").read_text().splitlines())
            self.assertIn(status["status"], ("tested", "insufficient_data"), tool)
            rows = list(read_tsv(folders[0] / "hisat_delta.tsv.gz"))
            # both circles, donor and acceptor; identical reads on both sides: no confident change
            self.assertEqual(len(rows), 4, tool)
            for row in rows:
                self.assertEqual(row["status"], "experimental", (tool, row))
                self.assertLess(max(float(row["probability_increase"]), float(row["probability_decrease"])), 0.95)
            evidence = list(read_tsv(folders[0] / "hisat_evidence.tsv.gz"))
            self.assertTrue(all(int(row["linear_unique"]) > 0 for row in evidence), tool)
        self.assertFalse(list(self.work.glob("comparisons/{}/project/unselected/modules/*/differential".format(self.ref))))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def analysis(self, tool):
        folders = list((self.work / "comparisons" / self.ref / "project/pilot/modules" / tool).iterdir())
        self.assertEqual(len(folders), 1)
        return folders[0]

    def test_workflow_completes(self):
        self.assertEqual(self.result.returncode, 0, self.result.stderr[-6000:])

    def test_catalogue_lifts_hg19_and_merges_sources(self):
        reports = list(self.work.glob("umbrella/modules/reference/{}/circrna_catalogue/*/report.json".format(self.ref)))
        self.assertEqual(len(reports), 1)
        report = json.loads(reports[0].read_text())
        self.assertEqual(report["kept"], 3)
        lift = json.loads(next(reports[0].parent.glob("sources/1.lift.json")).read_text())
        self.assertEqual((lift["lifted"], lift["failed"]), (1, {"contig_not_in_chain": 1}))
        with gzip.open(reports[0].parent / "annotation.tsv.gz", "rt") as stream:
            names = {line.split("\t")[0]: line.split("\t")[6] for line in stream}
        self.assertEqual(names["chr2:2201|2400"], "circAtlas_C2,circBase_C2")

    def test_ciriquant_comparison_reports_both_circles(self):
        rows = {row["feature_id"]: row for row in read_tsv(self.analysis("ciriquant") / "results.tsv.gz")}
        self.assertEqual(set(rows), set(self.fx["expected_bsj"]))
        self.assertEqual(rows["chr1:1001|1800"]["effect"], "0.000000")     # same reads on both sides
        status = json.loads((self.analysis("ciriquant") / "status.json").read_text())
        self.assertEqual((status["status"], status["selection"]["role"]), ("no_tests", "test"))

    def test_discovery_from_kept_star_caches(self):
        support = {row["circ_id"]: row for row in read_tsv(self.analysis("circexplorer2") / "support.tsv.gz")}
        self.assertEqual(set(support), set(self.truth))
        self.assertEqual(support["chr1:1001|1800"]["catalogue_supported"], "yes")
        dcc = {row["circ_id"]: row for row in read_tsv(self.analysis("dcc") / "support.tsv.gz")}
        self.assertEqual(set(dcc), set(self.truth))
        self.assertEqual(dcc["chr2:2201|2400"]["strand"], "-")

    def test_only_selected_comparison_and_runs(self):
        self.assertFalse((self.work / "comparisons" / self.ref / "project/unselected/modules").exists())
        self.assertFalse(list(self.work.glob("umbrella/modules/cache/{}/project/batch/run_c1/*".format(self.ref))))
        inventory = json.loads(next(self.work.glob("umbrella/modules/inventory/*.json")).read_text())
        self.assertEqual(set(inventory["comparisons"]), {"pilot"})
        self.assertEqual(inventory["scope"], {"ciriquant": "selected", "circexplorer2": "selected", "dcc": "selected"})
        self.assertEqual(inventory["shared_evidence"], {"circexplorer2": "circrna_star", "dcc": "circrna_star"})


if __name__ == "__main__":
    unittest.main()
