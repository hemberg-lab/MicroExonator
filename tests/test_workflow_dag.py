import os
import shutil
import subprocess
import tempfile
import unittest


REPOSITORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find_snakemake():
    configured = os.environ.get("MICROEXONATOR_SNAKEMAKE")
    if configured:
        return configured
    return shutil.which("snakemake")


class WorkflowSelectionTests(unittest.TestCase):
    def run_quant_dry_run(
        self,
        filter_method=None,
        single_cell=False,
        clusters=None,
        include_bulk_manifest=True,
        skip_discovery_and_quant=False,
        explicit_cluster_columns=True,
        extra_config=(),
        print_shell=False,
        target="quant",
        delta_comparisons=None,
        run_metadata=None,
        annotation_gtf=True,
    ):
        snakemake = find_snakemake()
        if not snakemake:
            self.skipTest("set MICROEXONATOR_SNAKEMAKE to run workflow DAG tests")

        with tempfile.TemporaryDirectory() as temp_dir:
            for name in ("MicroExonator.smk",):
                shutil.copy2(os.path.join(REPOSITORY, name), temp_dir)
            for name in ("rules", "src", "envs"):
                shutil.copytree(os.path.join(REPOSITORY, name), os.path.join(temp_dir, name))

            if single_cell and clusters is not None:
                sample_names = list(clusters)
            else:
                sample_names = ["cell_a", "cell_b"] if single_cell else ["sample_a"]

            for name in (
                "genome.fa",
                "annotation.bed12",
                "annotation.gtf",
                "u2_5.matrix",
                "u2_3.matrix",
                "conservation.bw",
                "microexons.bed12",
            ):
                open(os.path.join(temp_dir, name), "w").close()

            for sample in sample_names:
                open(os.path.join(temp_dir, "{}.fastq.gz".format(sample)), "w").close()

            with open(os.path.join(temp_dir, "local_samples.tsv"), "w") as handle:
                handle.write("sample\tpath\n")
                for sample in sample_names:
                    handle.write(
                        "{}\t{}\n".format(
                            sample,
                            os.path.join(temp_dir, "{}.fastq.gz".format(sample)),
                        )
                    )
            if not single_cell and include_bulk_manifest:
                with open(os.path.join(temp_dir, "bulk_samples.tsv"), "w") as handle:
                    handle.write("sample\tcondition\n")
                    handle.write("sample_a\tcontrol\n")
            if clusters is not None:
                with open(os.path.join(temp_dir, "clusters.tsv"), "w") as handle:
                    cell_column = "cell" if explicit_cluster_columns else "sample"
                    handle.write("{}\tcluster\n".format(cell_column))
                    for sample in sample_names:
                        handle.write("{}\t{}\n".format(sample, clusters[sample]))
            if skip_discovery_and_quant:
                existing_files = ["Round2/TOTAL.ME_centric.txt"]
                for sample in sample_names:
                    existing_files.extend(
                        [
                            "Report/quant/{}.out_filtered_ME.PSI.uncorrected.gz".format(
                                sample
                            ),
                            "Report/quant/corrected/counts/{}.ME.adj_counts.gz".format(
                                sample
                            ),
                            "Round2/ME_reads/{}.counts.tsv".format(sample),
                        ]
                    )
                for relative_path in existing_files:
                    absolute_path = os.path.join(temp_dir, relative_path)
                    os.makedirs(os.path.dirname(absolute_path), exist_ok=True)
                    open(absolute_path, "w").close()
            if run_metadata is not None:
                with open(os.path.join(temp_dir, "run_metadata.tsv"), "w") as handle:
                    handle.write(run_metadata)
            if delta_comparisons is not None:
                with open(os.path.join(temp_dir, "whippet.delta.yaml"), "w") as handle:
                    handle.write(delta_comparisons)
            with open(os.path.join(temp_dir, "config.yaml"), "w") as handle:
                handle.write("Genome_fasta: genome.fa\n")
                handle.write("Gene_anontation_bed12: annotation.bed12\n")
                if annotation_gtf:
                    handle.write("Gene_anontation_GTF: annotation.gtf\n")
                handle.write("GT_AG_U2_5: u2_5.matrix\n")
                handle.write("GT_AG_U2_3: u2_3.matrix\n")
                handle.write("conservation_bigwig: conservation.bw\n")
                handle.write("ME_DB: microexons.bed12\n")
                handle.write("working_directory: {}/\n".format(temp_dir))
                if not single_cell and include_bulk_manifest:
                    handle.write("bulk_samples: bulk_samples.tsv\n")
                if clusters is not None:
                    handle.write("cluster_metadata: clusters.tsv\n")
                    if explicit_cluster_columns:
                        handle.write("file_basename: cell\n")
                        handle.write("cluster_name: cluster\n")
                handle.write("min_number_files_detected: 1\n")
                handle.write("Single_Cell: {}\n".format("T" if single_cell else "F"))
                if skip_discovery_and_quant:
                    handle.write("skip_discovery_and_quant: T\n")
                if filter_method is not None:
                    handle.write("filter_method: {}\n".format(filter_method))
                if delta_comparisons is not None:
                    handle.write("whippet_delta: whippet.delta.yaml\n")
                if run_metadata is not None:
                    handle.write("run_metadata: run_metadata.tsv\n")
                for line in extra_config:
                    handle.write(line + "\n")

            result = subprocess.run(
                [snakemake, "-s", "MicroExonator.smk", "-n", "-j", "1"]
                + (["-p"] if print_shell else [])
                + [target],
                cwd=temp_dir,
                env=dict(os.environ, XDG_CACHE_HOME=os.path.join(temp_dir, ".cache")),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

        return result

    def test_quant_uses_robustness_output_by_default(self):
        result = self.run_quant_dry_run()

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/out.robustly_detected.txt", result.stdout)
        self.assertNotIn("Report/out.high_quality.txt", result.stdout)
        self.assertIn("Report/ME_junction_genomic_copies.txt", result.stdout)
        self.assertIn("Report/read_lengths.tsv", result.stdout)

    def test_legacy_quant_runs_mixture_preflight(self):
        result = self.run_quant_dry_run(filter_method="legacy_mixture")

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/out.high_quality.txt", result.stdout)
        self.assertIn("Report/.legacy_mixture_input.valid", result.stdout)

    def test_unclustered_single_cell_quant_uses_all_cells_group(self):
        result = self.run_quant_dry_run(single_cell=True)

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/filter/sc/all_cells.detected.txt", result.stdout)

    def test_clustered_single_cell_quant_uses_clusters_as_filter_groups(self):
        result = self.run_quant_dry_run(
            single_cell=True,
            clusters={"cell_a": "Excitatory neuron", "cell_b": "Glia"},
        )

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/filter/sc/Excitatory_neuron.detected.txt", result.stdout)
        self.assertIn("Report/filter/sc/Glia.detected.txt", result.stdout)

    def test_cluster_metadata_uses_default_column_names(self):
        result = self.run_quant_dry_run(
            single_cell=True,
            clusters={"cell_a": "Neuron", "cell_b": "Glia"},
            explicit_cluster_columns=False,
        )

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/filter/sc/Neuron.detected.txt", result.stdout)

    def test_bulk_quant_requires_bulk_samples_manifest(self):
        result = self.run_quant_dry_run(include_bulk_manifest=False)

        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("bulk_samples.tsv", result.stdout)

    def test_quant_runs_with_a_bed12_annotation_only(self):
        result = self.run_quant_dry_run(annotation_gtf=False)

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/out.robustly_detected.txt", result.stdout)

    def test_legacy_mixture_quant_needs_no_bulk_samples(self):
        result = self.run_quant_dry_run(include_bulk_manifest=False, filter_method="legacy_mixture")

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/out.high_quality.txt", result.stdout)
        self.assertNotIn("Report/filter/se/", result.stdout)

    def test_start_up_prints_no_deprecation_noise(self):
        # a deprecated config key gives one plain note; nothing that looks like an error
        result = self.run_quant_dry_run(extra_config=["filter_mode: unbiased"])
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertEqual(out.count("Note: filter_mode is deprecated; use filter_method: robustness"), 1, out[:2000])
        self.assertNotIn("DeprecationWarning", out, out[:2000])
        self.assertNotIn("invalid escape sequence", out, out[:2000])
        self.assertNotIn("raise WorkflowError", out, out[:2000])

    def test_existing_quantifications_can_use_default_robustness_filter(self):
        result = self.run_quant_dry_run(skip_discovery_and_quant=True)

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/out.robustly_detected.txt", result.stdout)
        self.assertNotIn("Report/ME_junction_genomic_copies.txt", result.stdout)
        self.assertNotIn("Report/read_lengths.tsv", result.stdout)


class DeltaMethodTests(unittest.TestCase):
    """delta_method selects Whippet or the Whippet-free model for differential_inclusion."""

    COMPARISON = "control_vs_case:\n  A: sample_a\n  B: sample_a\n"

    def delta_dry_run(self, *extra_config):
        return WorkflowSelectionTests.run_quant_dry_run(
            self,
            extra_config=extra_config,
            print_shell=True,
            target="differential_inclusion",
            delta_comparisons=self.COMPARISON,
        )

    def test_microexonator_delta_needs_no_whippet(self):
        result = self.delta_dry_run("delta_method: microexonator")

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("src/me_delta.py", result.stdout)
        self.assertIn("Delta/control_vs_case.diff.ME.microexons", result.stdout)
        self.assertIn("Report/quant/corrected/PSI_sparse/bulk/se/sample_a.corrected.PSI.gz", result.stdout)
        self.assertNotIn("whippet-", result.stdout)
        self.assertNotIn("julia", result.stdout)

    def test_microexonator_delta_is_the_default(self):
        result = self.delta_dry_run("whippet_bin_folder: /opt/whippet/bin", "julia: julia")

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("src/me_delta.py", result.stdout)
        self.assertNotIn("whippet-delta.jl", result.stdout)

    def test_whippet_delta_on_request(self):
        result = self.delta_dry_run("delta_method: whippet", "whippet_bin_folder: /opt/whippet/bin", "julia: julia")

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("whippet-delta.jl", result.stdout)
        self.assertNotIn("src/me_delta.py", result.stdout)

    def test_downstream_only_alone_runs_whippet_on_the_annotation(self):
        # the MiMB chapter (section 4.1.2) sets downstream_only and nothing else
        result = self.delta_dry_run("downstream_only: T", "whippet_bin_folder: /opt/whippet/bin", "julia: julia")

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("whippet-index.jl --fasta genome.fa --gtf annotation.gtf", result.stdout)
        self.assertIn("whippet-delta.jl", result.stdout)
        self.assertIn("Whippet/Delta/control_vs_case.diff.gz", result.stdout)
        self.assertNotIn("src/me_delta.py", result.stdout)
        self.assertNotIn("Report/out.robustly_detected", result.stdout)

    def test_downstream_only_needs_no_bulk_samples(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self,
            include_bulk_manifest=False,
            extra_config=("downstream_only: T", "whippet_bin_folder: /opt/whippet/bin", "julia: julia"),
            print_shell=True,
            target="differential_inclusion",
            delta_comparisons=self.COMPARISON,
        )

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("whippet-delta.jl", result.stdout)

    def test_downstream_only_keeps_working_with_the_old_settings(self):
        result = self.delta_dry_run("downstream_only: T", "Only_whippet: T", "delta_method: whippet",
                                    "whippet_bin_folder: /opt/whippet/bin", "julia: julia")

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("whippet-delta.jl", result.stdout)

    def test_downstream_only_refuses_the_microexonator_delta(self):
        result = self.delta_dry_run("downstream_only: T", "delta_method: microexonator",
                                    "whippet_bin_folder: /opt/whippet/bin", "julia: julia")

        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("downstream_only skips the quantification", result.stdout)

    def test_unknown_delta_method_is_rejected(self):
        result = self.delta_dry_run("delta_method: rmats")

        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("delta_method", result.stdout)


class QuantificationOnlyRouteTests(unittest.TestCase):
    """skip_discovery + only_db: the fixed-universe route used by MegaSearch."""

    def annotation_command(self, *extra_config):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self,
            extra_config=("skip_discovery: T", "only_db: T") + extra_config,
            print_shell=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        commands = [
            line.strip()
            for line in result.stdout.splitlines()
            if "src/Get_annotated_microexons.py" in line
        ]
        self.assertEqual(len(commands), 1, result.stdout)
        arguments = commands[0].split()
        return arguments[arguments.index("src/Get_annotated_microexons.py"):]

    def test_synthetic_skip_tags_are_on_by_default(self):
        # sys.argv[9] of Get_annotated_microexons.py is the skip-tag source.
        self.assertEqual(self.annotation_command()[9], "Round1/ME_TAGs.fa")

    def test_synthetic_skip_tags_can_be_disabled(self):
        self.assertEqual(
            self.annotation_command("synthetic_skip_tags: F")[9], "NA"
        )


class SingleCellModuleTests(unittest.TestCase):
    """Snakepool.py and pseudo_pool.smk: the single-cell targets of the MiMB chapter."""

    CELLS = {"n1": "Neuron", "n2": "Neuron", "g1": "Glia", "g2": "Glia"}
    WHIPPET = ("whippet_bin_folder: /opt/whippet/bin", "julia: /opt/julia/bin/julia")
    # C1: two pools per group, one repeat; C2: one pool per group, two repeats
    RUN_METADATA = (
        "Compare_ID\tA.cluster_names\tA.number_of_pools\tB.cluster_names\tB.number_of_pools\tRepeat\n"
        "C1\tNeuron\t2\tGlia\t2\t1\n"
        "C2\tNeuron\t1\tGlia\t1\t2\n"
    )

    def dry_run(self, target, run_metadata=RUN_METADATA, extra_config=()):
        return WorkflowSelectionTests.run_quant_dry_run(
            self,
            single_cell=True,
            clusters=self.CELLS,
            extra_config=self.WHIPPET + tuple(extra_config),
            print_shell=True,
            target=target,
            run_metadata=run_metadata,
        )

    def test_every_chapter_target_builds(self):
        for target in ("snakepool", "quant_unpool_single_cell", "collapse_whippet",
                       "cluster_bams", "collapse_pseudo_pools"):
            with self.subTest(target=target):
                result = self.dry_run(target)
                self.assertEqual(result.returncode, 0, result.stdout[-3000:])

    def test_snakepool_uses_configured_julia_and_filter_output(self):
        out = self.dry_run("snakepool").stdout
        self.assertIn("/opt/julia/bin/julia /opt/whippet/bin/whippet-delta.jl", out)
        self.assertIn("/opt/julia/bin/julia /opt/whippet/bin/whippet-quant.jl", out)
        self.assertNotRegex(out, r"(?m)^\s*julia ")
        self.assertIn("python3 src/get_diff_ME_single_cell.py", out)
        self.assertIn("Report/out.robustly_detected.txt", out)
        self.assertNotIn("Report/out.high_quality.txt", out)

    def test_each_comparison_keeps_its_own_pools_and_repeats(self):
        out = self.dry_run("snakepool").stdout
        delta = [line for line in out.splitlines() if "whippet-delta.jl" in line]
        c1 = [line for line in delta if "-o Whippet/Delta/Single_Cell/C1_rep_" in line]
        c2 = [line for line in delta if "-o Whippet/Delta/Single_Cell/C2_rep_" in line]
        self.assertEqual(len(c1), 1, delta)
        self.assertEqual(len(c2), 2, delta)
        # two pools per group in C1, one in C2
        self.assertEqual(c1[0].split(" -a ")[1].split()[0].count(","), 1, c1[0])
        for line in c2:
            self.assertEqual(line.split(" -a ")[1].split()[0].count(","), 0, line)
        self.assertNotIn("C1_rep_2", out)

    def test_run_metadata_is_only_needed_for_snakepool(self):
        result = self.dry_run("quant_unpool_single_cell", run_metadata=None)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])

    def test_start_up_does_not_print_cluster_sizes(self):
        out = self.dry_run("snakepool").stdout
        self.assertNotRegex(out, r"(?m)^(Neuron|Glia) 2$")


if __name__ == "__main__":
    unittest.main()
