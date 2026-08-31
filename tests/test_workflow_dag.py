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
    ):
        snakemake = find_snakemake()
        if not snakemake:
            self.skipTest("set MICROEXONATOR_SNAKEMAKE to run workflow DAG tests")

        with tempfile.TemporaryDirectory() as temp_dir:
            for name in ("MicroExonator.smk",):
                shutil.copy2(os.path.join(REPOSITORY, name), temp_dir)
            for name in ("rules", "src", "envs"):
                shutil.copytree(os.path.join(REPOSITORY, name), os.path.join(temp_dir, name))

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
            with open(os.path.join(temp_dir, "config.yaml"), "w") as handle:
                handle.write("Genome_fasta: genome.fa\n")
                handle.write("Gene_anontation_bed12: annotation.bed12\n")
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

            result = subprocess.run(
                [snakemake, "-s", "MicroExonator.smk", "-n", "-j", "1", "quant"],
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

    def test_existing_quantifications_can_use_default_robustness_filter(self):
        result = self.run_quant_dry_run(skip_discovery_and_quant=True)

        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Report/out.robustly_detected.txt", result.stdout)


if __name__ == "__main__":
    unittest.main()
