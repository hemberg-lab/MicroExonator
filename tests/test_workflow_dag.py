import os
import re
import shutil
import subprocess
import sys
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
        umbrella=False,
        umbrella_layout="PE",
        umbrella_prebuilt=True,
        umbrella_derived_reference=False,
        umbrella_extra_rows=(),
        umbrella_comparisons=None,
        rulegraph=False,
        umbrella_reference_id="ref",
        pre_existing=(),
        prepare=None,
        umbrella_extra_column=None,
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
            if delta_comparisons is not None:
                with open(os.path.join(temp_dir, "whippet.delta.yaml"), "w") as handle:
                    handle.write(delta_comparisons)
            if umbrella:
                with open(os.path.join(temp_dir, "umbrella.tsv"), "w") as handle:
                    extra = "\t" + umbrella_extra_column if umbrella_extra_column else ""
                    handle.write("sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\tsource_type\tsource_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude" + extra + "\n")
                    source_2 = "r2.fastq.gz" if umbrella_layout == "PE" else ""
                    handle.write("sample_a\trun_{}\trep_a\tproject\tbatch\tcontrol\tfastq\tr1.fastq.gz\t{}\t{}\tunstranded\t{}\ttrue{}\n".format(
                        umbrella_layout.lower(), source_2, umbrella_layout, umbrella_reference_id,
                        "\t" if umbrella_extra_column else ""))
                    for row in umbrella_extra_rows:
                        handle.write(row + "\n")
                if umbrella_comparisons is not None:
                    with open(os.path.join(temp_dir, "comparisons.yaml"), "w") as handle:
                        handle.write(umbrella_comparisons)
                for name in ("r1.fastq.gz", "r2.fastq.gz", "transcripts.fa", "decoys.txt", "fixed.gtf", "salmon.gtf", "splice_sites.txt"):
                    open(os.path.join(temp_dir, name), "w").close()
                if umbrella_prebuilt:
                    os.mkdir(os.path.join(temp_dir, "salmon_index"))
                    open(os.path.join(temp_dir, "salmon_index", "hash.bin"), "w").close()
                    for name in ("whippet.jls", "whippet.jls.exons.tab.gz"):
                        open(os.path.join(temp_dir, name), "w").close()
                    for number in range(1, 9):
                        open(os.path.join(temp_dir, "hisat.{}.ht2".format(number)), "w").close()
            for relative_path, text in pre_existing:
                absolute_path = os.path.join(temp_dir, relative_path)
                os.makedirs(os.path.dirname(absolute_path), exist_ok=True)
                if relative_path.endswith(".gz"):
                    import gzip
                    with gzip.open(absolute_path, "wt") as handle:
                        handle.write(text)
                else:
                    with open(absolute_path, "w") as handle:
                        handle.write(text)
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
                if delta_comparisons is not None:
                    handle.write("whippet_delta: whippet.delta.yaml\n")
                if umbrella:
                    handle.write("umbrella_manifest: umbrella.tsv\n")
                    if umbrella_comparisons is not None:
                        handle.write("umbrella_comparisons: comparisons.yaml\n")
                    handle.write("whippet_bin_folder: /opt/whippet/bin\n")
                    handle.write("julia: julia\n")
                    handle.write("umbrella_reference:\n")
                    reference_files = (
                        ("genome_fasta", "genome.fa"),
                        ("annotation_gtf", "annotation.gtf"),
                        ("whippet_gtf", "fixed.gtf"),
                        ("me_db", "microexons.bed12"),
                    )
                    if not umbrella_derived_reference:
                        reference_files += (
                            ("transcriptome_fasta", "transcripts.fa"),
                            ("decoys", "decoys.txt"),
                            ("salmon_gtf", "salmon.gtf"),
                            ("splice_sites", "splice_sites.txt"),
                        )
                    for key, value in reference_files + ((
                        ("hisat2_index_prefix", "hisat"),
                        ("whippet_index", "whippet.jls"),
                        ("salmon_index", "salmon_index"),
                    ) if umbrella_prebuilt else ()):
                        handle.write("  {}: {}\n".format(key, value))
                for line in extra_config:
                    handle.write(line + "\n")

            if prepare is not None:
                prepare(temp_dir)
            command = [snakemake]
            if sys.platform == "darwin":
                # appdirs ignores XDG_CACHE_HOME on macOS; keep Snakemake's
                # runtime source cache inside this isolated fixture.
                command = [os.path.join(os.path.dirname(snakemake), "python"), "-c",
                           "import appdirs, runpy; appdirs.system = 'linux'; runpy.run_module('snakemake', run_name='__main__')"]
            result = subprocess.run(
                command + ["-s", "MicroExonator.smk"]
                + (["--rulegraph"] if rulegraph else ["-n", "-j", "1"])
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


class UmbrellaNativeQuantificationTests(unittest.TestCase):
    def test_minimal_reference_derives_transcripts_decoys_and_splice_hints(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, umbrella_prebuilt=False,
            umbrella_derived_reference=True, print_shell=True,
            target="quant_umbrella"
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        for rule in ("umbrella_transcriptome_fasta", "umbrella_decoy_names",
                     "umbrella_splice_hints", "umbrella_salmon_index"):
            self.assertIn("rule {}:".format(rule), result.stdout)
        self.assertIn("gffread -w umbrella/reference/ref/salmon/transcripts.fa", result.stdout)
        self.assertIn("-d umbrella/reference/ref/salmon/decoys.txt", result.stdout)

    def test_built_reference_declares_all_index_outputs(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, umbrella_prebuilt=False, print_shell=True,
            target="quant_umbrella"
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("rule umbrella_hisat2_index:", result.stdout)
        self.assertIn("umbrella/reference/ref/hisat2/genome.8.ht2", result.stdout)
        self.assertIn("rule umbrella_salmon_index:", result.stdout)
        self.assertIn("-d decoys.txt", result.stdout)
        self.assertIn("rule umbrella_whippet_index:", result.stdout)
        self.assertIn("--gtf umbrella/reference/ref/whippet.microexons.gtf.gz", result.stdout)
        self.assertNotIn("Report/out.robustly_detected.gtf", result.stdout)

    def test_pe_uses_native_mates_and_never_legacy_se_salmon_route(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target="quant_umbrella"
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("whippet-quant.jl umbrella/work/ref/project/batch/run_pe/R1.fastq.gz umbrella/work/ref/project/batch/run_pe/R2.fastq.gz", result.stdout)
        self.assertIn("salmon quant", result.stdout)
        self.assertIn("-1 umbrella/work/ref/project/batch/run_pe/R1.fastq.gz", result.stdout)
        self.assertIn("-2 umbrella/work/ref/project/batch/run_pe/R2.fastq.gz", result.stdout)
        self.assertIn("--validateMappings", result.stdout)
        self.assertNotIn("salmon/SE/run_pe", result.stdout)
        self.assertNotIn("lengthScaledTPM", result.stdout)
        self.assertNotIn("Report/out.robustly_detected.gtf", result.stdout)
        self.assertNotIn("rule umbrella_hisat2_index:", result.stdout)
        self.assertNotIn("salmon index", result.stdout)
        self.assertNotIn("whippet-index.jl", result.stdout)

    def test_shared_alignment_pe_uses_native_mates_hints_and_fragment_counts(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target="quant_umbrella"
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertIn("hisat2 -p 1 -x hisat -1 umbrella/work/ref/project/batch/run_pe/R1.fastq.gz "
                      "-2 umbrella/work/ref/project/batch/run_pe/R2.fastq.gz "
                      "--known-splicesite-infile umbrella/reference/ref/hisat2_splice_sites.txt", out)
        # ME_DB microexons are added to derived GTF copies, never to the configured files
        self.assertIn("src/insert_microexons_gtf.py --gtf annotation.gtf", out)
        self.assertIn("src/insert_microexons_gtf.py --gtf fixed.gtf", out)
        self.assertNotIn("--dta", out)
        self.assertIn("-p --countReadPairs", out)
        self.assertIn("-t paired", out)
        for shard in ("junctions/ref/project/control/batch.junctions.tsv.gz",
                      "junctions/ref/project/control/batch.capture.tsv.gz",
                      "genes/ref/project/control/batch.featurecounts.tsv.gz",
                      "qc/ref/project/control/batch.qc.tsv.gz",
                      "rmats/ref/project/control/batch.rmats_inventory.tsv",
                      "coverage/ref/project/control/batch.sum_cpm.bw"):
            self.assertIn(shard, out)
        # Every BAM consumer is a declared rule, so the BAM can be temporary.
        for rule in ("umbrella_junctions", "umbrella_featurecounts",
                     "umbrella_rmats_prep", "umbrella_coverage_run"):
            self.assertIn("rule {}:".format(rule), out)

    def test_shared_alignment_se_and_optional_analyses_can_be_disabled(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, umbrella_layout="SE", print_shell=True,
            target="quant_umbrella",
            extra_config=("umbrella_optional:", "  rmats: false", "  coverage: false"),
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertIn("-U umbrella/work/ref/project/batch/run_se/R1.fastq.gz", out)
        self.assertNotIn("--countReadPairs", out)
        self.assertNotIn("rmats.py", out)
        self.assertNotIn("genomecov", out)
        self.assertIn("junctions/ref/project/control/batch.junctions.tsv.gz", out)

    def test_legacy_quant_does_not_reach_umbrella_rules(self):
        result = WorkflowSelectionTests.run_quant_dry_run(self, print_shell=True)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertNotIn("umbrella", result.stdout)
        self.assertNotIn("hisat2 -p", result.stdout)

    def test_umbrella_reports_detection_without_blocking_delta(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
                                        self, umbrella=True, target="quant_umbrella", rulegraph=True,
                                        umbrella_comparisons=(
                                            "comparisons:\n  - comparison_id: case_vs_control\n"
                                            "    project_id: project\n    group_a: case\n    group_b: control\n"),
                                        umbrella_extra_rows=[
                                            "sample_b\trun_b\trep_b\tproject\tbatch\tcase\tfastq\tr1.fastq.gz\t"
                                            "r2.fastq.gz\tPE\tunstranded\tref\ttrue"])
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn('label = "umbrella_fixed_microexons"', result.stdout)
        self.assertIn('label = "umbrella_microexonator_delta"', result.stdout)
        self.assertIn('label = "detection_filter"', result.stdout)
        detection = re.search(r'^\s*(\d+)\[label = "detection_filter"', result.stdout, re.M).group(1)
        target = re.search(r'^\s*(\d+)\[label = "quant_umbrella"', result.stdout, re.M).group(1)
        self.assertEqual(re.findall(r'^\s*' + detection + r' -> (\d+)$', result.stdout, re.M), [target])

    def test_reference_identity_gates_index_builds_and_sra_downloads_are_capped(self):
        rows = ["sample_s\tSRR000001\trep_s\tproject\tbatch\tcontrol\tsra\tSRR000001\t\tPE\tunstranded\tref\ttrue"]
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, umbrella_prebuilt=False, print_shell=True,
            target="quant_umbrella", umbrella_extra_rows=rows)
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertIn("rule umbrella_reference_identity:", out)
        # the identity check runs before any index is built
        self.assertLess(out.index("rule umbrella_reference_identity:"), out.index("rule umbrella_hisat2_index:"))
        self.assertLess(out.index("rule umbrella_reference_identity:"), out.index("rule umbrella_salmon_index:"))
        staging = [block for block in out.split("\n\n") if "rule umbrella_stage_reads:" in block]
        self.assertTrue(any("SRR000001" in block and "get_data=1" in block for block in staging), staging)
        self.assertTrue(any("run_pe" in block and "get_data=0" in block for block in staging), staging)
        # rMATS reads only an uncompressed GTF
        self.assertIn("gzip -dcf umbrella/reference/ref/annotation.microexons.gtf.gz > umbrella/reference/ref/annotation.rmats.gtf", out)
        self.assertIn("--gtf umbrella/reference/ref/annotation.rmats.gtf", out)
        self.assertNotIn("--gtf umbrella/reference/ref/annotation.microexons.gtf.gz --od", out)

    def test_reference_id_auto_is_computed_by_the_workflow(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, umbrella_prebuilt=False, print_shell=True,
            target="quant_umbrella", umbrella_reference_id="auto")
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        match = re.search(r"Umbrella reference_id: ([0-9a-f]{16}) \(computed from umbrella_reference\)", out)
        self.assertIsNotNone(match, out[-3000:])
        # every umbrella path uses the computed ID, none the placeholder
        self.assertIn("umbrella/reference/{}/".format(match.group(1)), out)
        self.assertNotIn("umbrella/reference/auto/", out)
        self.assertNotIn("umbrella/work/auto/", out)

    def test_microexonator_whippet_delta_is_optional_and_feeds_the_synthesis(self):
        rows = ["sample_{0}\trun_{0}\trep_{0}\tproject\tbatch\t{1}\tfastq\tr1.fastq.gz\t"
                "r2.fastq.gz\tPE\tunstranded\tref\ttrue".format(name, group)
                for name, group in (("b", "control"), ("c", "case"), ("d", "case"))]
        comparisons = ("comparisons:\n  - comparison_id: case_vs_control\n    project_id: project\n"
                       "    group_a: case\n    group_b: control\n")
        runs = {}
        for enabled in (False, True):
            extra = ["umbrella_optional:", "  microexonator_whippet_delta: true"] if enabled else []
            runs[enabled] = WorkflowSelectionTests.run_quant_dry_run(
                self, umbrella=True, print_shell=True, target="quant_umbrella",
                umbrella_extra_rows=rows, umbrella_comparisons=comparisons, extra_config=extra)
            self.assertEqual(runs[enabled].returncode, 0, runs[enabled].stdout)
        self.assertNotIn("rule umbrella_microexonator_whippet_delta:", runs[False].stdout)
        out = runs[True].stdout
        self.assertIn("rule umbrella_microexonator_whippet_delta:", out)
        self.assertIn("run_comparison_tools.py microexonator_whippet", out)
        # the Whippet-free delta still runs, with the empty options kept attached to the flag
        self.assertIn("--options= --outputs", out)
        self.assertIn("logs/microexonator_delta.log", out)

    def test_start_up_prints_no_deprecation_noise(self):
        # a deprecated config key gives one plain note; nothing that looks like an error
        result = WorkflowSelectionTests.run_quant_dry_run(self, extra_config=["filter_mode: unbiased"])
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertEqual(out.count("Note: filter_mode is deprecated; use filter_method: robustness"), 1, out[:2000])
        self.assertNotIn("DeprecationWarning", out, out[:2000])
        self.assertNotIn("invalid escape sequence", out, out[:2000])
        self.assertNotIn("raise WorkflowError", out, out[:2000])

    def test_splicing_shards_and_comparison_joins(self):
        rows = ["sample_{0}\trun_{0}\trep_{0}\tproject\tbatch\t{1}\tfastq\tr1.fastq.gz\t"
                "r2.fastq.gz\tPE\tunstranded\tref\ttrue".format(name, group)
                for name, group in (("b", "control"), ("c", "case"), ("d", "case"))]
        comparisons = ("comparisons:\n  - comparison_id: case_vs_control\n    project_id: project\n"
                       "    group_a: case\n    group_b: control\n")
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target="quant_umbrella",
            umbrella_extra_rows=rows, umbrella_comparisons=comparisons,
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertIn("splicing/ref/project/control/batch.microexonator.tsv.gz", out)
        self.assertIn("splicing/ref/project/case/batch.whippet.tsv.gz", out)
        self.assertIn("comparisons/ref/project/case_vs_control/preflight.json", out)
        self.assertIn("--kind junctions", out)
        self.assertIn("--collapse comparisons/ref/project/case_vs_control/preflight.json", out)
        self.assertIn("comparisons/ref/project/case_vs_control/joined/whippet.tsv.gz", out)

    def test_comparison_tools_synthesis_and_target_scopes(self):
        rows = ["sample_{0}\trun_{0}\trep_{0}\tproject\tbatch\t{1}\tfastq\tr1.fastq.gz\t"
                "r2.fastq.gz\tPE\tunstranded\tref\ttrue".format(name, group)
                for name, group in (("b", "control"), ("c", "case"), ("d", "case"))]
        comparisons = ("comparisons:\n  - comparison_id: case_vs_control\n    project_id: project\n"
                       "    group_a: case\n    group_b: control\n")
        runs = {}
        for target in ("quant_microexonator", "quant_whippet", "quant_umbrella"):
            runs[target] = WorkflowSelectionTests.run_quant_dry_run(
                self, umbrella=True, print_shell=True, target=target,
                umbrella_extra_rows=rows, umbrella_comparisons=comparisons)
            self.assertEqual(runs[target].returncode, 0, runs[target].stdout)
        full = runs["quant_umbrella"].stdout
        root = "comparisons/ref/project/case_vs_control/"
        for path in ("synthesis.tsv", "deseq2_featurecounts.results.tsv.gz",
                     "deseq2_tximport.results.tsv.gz", "microexonator_delta.tsv",
                     "whippet_delta.diff.gz", "rmats/SE.MATS.JC.txt",
                     "leafcutter_cluster_significance.txt", "suppa2/SE.dpsi"):
            self.assertIn(root + path, full)
        self.assertIn("src/expression_deseq2.R --route tximport", full)
        self.assertIn("src/run_comparison_tools.py rmats", full)
        self.assertIn("src/run_leafcutter.py --preflight", full)
        self.assertIn("leafcutter-gtf-to-exons", full)
        self.assertIn(root + "leafcutter_effect_sizes.txt", full)
        self.assertIn(root + "leafcutter_introns.counts.gz", full)
        self.assertNotIn("leafcutter_dir", full)
        self.assertIn("suppa.py generateEvents", full)
        # narrower targets stop before alignment and comparisons
        self.assertIn("splicing/ref/project/case/batch.microexonator.tsv.gz", runs["quant_microexonator"].stdout)
        self.assertNotIn("whippet-quant.jl", runs["quant_microexonator"].stdout)
        self.assertIn("splicing/ref/project/case/batch.whippet.tsv.gz", runs["quant_whippet"].stdout)
        for narrow in ("quant_microexonator", "quant_whippet"):
            self.assertNotIn("hisat2 -p", runs[narrow].stdout)
            self.assertNotIn(root, runs[narrow].stdout)

    SELECTION_ROWS = ["sample_{0}\trun_{0}\trep_{0}\tproject\tbatch\t{1}\tfastq\tr1.fastq.gz\t"
                      "r2.fastq.gz\tPE\tunstranded\tref\ttrue".format(name, group)
                      for name, group in (("b", "control"), ("c", "case"), ("d", "case"))]
    SELECTION_COMPARISONS = ("comparisons:\n  - comparison_id: case_vs_control\n    project_id: project\n"
                             "    group_a: case\n    group_b: control\n")
    # the umbrella's own indexes (MicroExonator's discovery builds its own HISAT2 index)
    INDEX_BUILDS = ("rule umbrella_hisat2_index:", "rule umbrella_salmon_index:",
                    "rule umbrella_whippet_index:", "rule umbrella_reference_manifest:")

    def selection_dry_run(self, target, extra=(), prebuilt=True):
        return WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target=target, umbrella_prebuilt=prebuilt,
            umbrella_extra_rows=self.SELECTION_ROWS, umbrella_comparisons=self.SELECTION_COMPARISONS,
            extra_config=extra)

    def test_differential_inclusion_runs_microexonator_and_its_delta_only(self):
        result = self.selection_dry_run("differential_inclusion", prebuilt=False)
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertIn("comparisons/ref/project/case_vs_control/microexonator_delta.tsv", out)
        self.assertIn("run_comparison_tools.py microexonator ", out)
        self.assertIn("splicing/ref/project/case/batch.microexonator.tsv.gz", out)
        self.assertIn("Report/ME_ambiguous_positions.txt", out)
        for absent in ("hisat2 -p", "whippet-quant.jl", "salmon quant", "whippet_delta.diff.gz",
                       "synthesis.tsv") + self.INDEX_BUILDS:
            self.assertNotIn(absent, out)

    def test_differential_inclusion_whippet_delta_on_request(self):
        result = self.selection_dry_run("differential_inclusion", extra=("delta_method: whippet",))
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("comparisons/ref/project/case_vs_control/whippet_delta.diff.gz", result.stdout)
        self.assertNotIn("microexonator_delta.tsv", result.stdout)

    def test_umbrella_refuses_the_legacy_comparisons_file(self):
        result = self.selection_dry_run("differential_inclusion", extra=("whippet_delta: whippet.delta.yaml",))
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("umbrella_comparisons", result.stdout)

    def test_quant_umbrella_follows_the_tool_selection(self):
        result = self.selection_dry_run(
            "quant_umbrella", extra=("umbrella_optional:", "  whippet: false", "  salmon: false"))
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        root = "comparisons/ref/project/case_vs_control/"
        for present in ("hisat2 -p", root + "microexonator_delta.tsv", root + "deseq2_featurecounts.results.tsv.gz",
                        root + "leafcutter_cluster_significance.txt", root + "synthesis.tsv",
                        root + "joined/junctions.tsv.gz", "multiqc_report.html"):
            self.assertIn(present, out)
        for absent in ("whippet-quant.jl", "salmon quant", root + "whippet_delta.diff.gz",
                       root + "deseq2_tximport", root + "suppa2/", root + "joined/whippet.tsv.gz",
                       root + "joined/salmon_counts.tsv.gz", "--whippet-map", "--salmon "):
            self.assertNotIn(absent, out)
        self.assertIn("--hisat2 ", out)

    def test_microexonator_alone_through_the_selection_needs_no_index(self):
        result = self.selection_dry_run(
            "quant_umbrella", prebuilt=False,
            extra=("umbrella_optional:", "  whippet: false", "  salmon: false", "  hisat2: false",
                   "  multiqc: false"))
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertIn("comparisons/ref/project/case_vs_control/microexonator_delta.tsv", out)
        self.assertIn("comparisons/ref/project/case_vs_control/synthesis.tsv", out)
        for absent in ("hisat2 -p", "whippet-quant.jl", "salmon quant", "fastqc") + self.INDEX_BUILDS:
            self.assertNotIn(absent, out)

    def test_default_selection_keeps_the_existing_qc_rule(self):
        result = self.selection_dry_run("quant_umbrella")
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("--featurecounts umbrella/work", result.stdout)
        self.assertIn("--whippet-map umbrella/work", result.stdout)
        self.assertIn("whippet_delta.diff.gz", result.stdout)

    def test_tool_selection_errors_are_actionable(self):
        conflict = self.selection_dry_run(
            "quant_umbrella", extra=("umbrella_optional:", "  hisat2: false", "  leafcutter: true"))
        self.assertNotEqual(conflict.returncode, 0, conflict.stdout)
        self.assertIn("leafcutter needs hisat2", conflict.stdout)
        unknown = self.selection_dry_run("quant_umbrella", extra=("umbrella_optional:", "  star: true"))
        self.assertNotEqual(unknown.returncode, 0, unknown.stdout)
        self.assertIn("unknown umbrella_optional star", unknown.stdout)
        nothing = self.selection_dry_run(
            "quant_umbrella", extra=("umbrella_optional:", "  microexonator: false", "  whippet: false",
                                     "  salmon: false", "  hisat2: false"))
        self.assertNotEqual(nothing.returncode, 0, nothing.stdout)
        self.assertIn("switches off every tool", nothing.stdout)

    UMBRELLA_ROWS = ["sample_{0}\trun_{0}\trep_{0}\tproject\tbatch\t{1}\tfastq\tr1.fastq.gz\t"
                     "r2.fastq.gz\tPE\tunstranded\tref\ttrue".format(name, group)
                     for name, group in (("b", "control"), ("c", "case"), ("d", "case"))]

    def test_new_runs_get_fastqc_kept_summaries_and_a_project_report(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target="quant_umbrella",
            umbrella_extra_rows=self.UMBRELLA_ROWS)
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertIn("qc/ref/project/case/batch/runs/run_c/fastqc", out)
        self.assertIn("src/umbrella_fastqc.py fastqc --run-id run_c", out)
        self.assertIn("--reads-per-file 2000000", out)
        self.assertIn("qc/ref/project/case/batch/runs/run_c/salmon/run_c/aux_info/meta_info.json", out)
        self.assertIn("multiqc/ref/project/multiqc_report.html", out)
        # MultiQC writes its data folder as <report name>_data
        self.assertIn("-n multiqc_report.html", out)
        self.assertIn("multiqc/ref/project/multiqc_report_data", out)

    def test_runs_processed_before_qc_was_kept_never_trigger_fastqc(self):
        # every shard of the project is already written: its reads are gone
        existing = [("qc/ref/project/{}/batch.qc.tsv.gz".format(group), "run_id\n")
                    for group in ("control", "case")]
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target="quant_umbrella",
            umbrella_extra_rows=self.UMBRELLA_ROWS, pre_existing=existing)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertNotIn("umbrella_fastqc", result.stdout)
        self.assertNotIn("umbrella_multiqc", result.stdout)
        # QC can be switched off altogether
        off = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, target="quant_umbrella", umbrella_extra_rows=self.UMBRELLA_ROWS,
            extra_config=["umbrella_optional:", "  multiqc: false"])
        self.assertEqual(off.returncode, 0, off.stdout)
        self.assertNotIn("umbrella_fastqc", off.stdout)

    def test_selection_comparison_joins_only_selected_runs(self):
        comparisons = ("comparisons:\n  - comparison_id: c_vs_rest\n    project_id: project\n"
                       "    a: {label: only_c, samples: [sample_c, sample_d]}\n"
                       "    b: {label: ctrl, groups: [control]}\n")
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target="quant_umbrella",
            umbrella_extra_rows=self.UMBRELLA_ROWS, umbrella_comparisons=comparisons)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("--select comparisons/ref/project/c_vs_rest/preflight.json", result.stdout)

    def test_relabelling_a_run_with_written_shards_is_refused(self):
        # run_c was stored under control; the manifest now says case
        existing = [("junctions/ref/project/control/batch.junctions.tsv.gz",
                     "chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\t"
                     "short_anchor_total\trun_b\trun_c\n")]
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, target="quant_umbrella",
            umbrella_extra_rows=self.UMBRELLA_ROWS, pre_existing=existing)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("run run_c is stored in shard ref/project/control/batch", result.stdout)
        self.assertIn("select on it in the comparisons file", result.stdout)

    MODULE_COMPARISONS = ("comparisons:\n  - comparison_id: case_vs_control\n    project_id: project\n"
                          "    group_a: case\n    group_b: control\n")

    def module_dry_run(self, target, extra=(), prepare=None, majiq_license=True):
        def setup(temp_dir):
            # MAJIQ finds a majiq_license* file in the working directory
            if majiq_license:
                open(os.path.join(temp_dir, "majiq_license_test.lic"), "w").close()
            if prepare:
                prepare(temp_dir)
        return WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target=target,
            umbrella_extra_rows=self.UMBRELLA_ROWS, umbrella_comparisons=self.MODULE_COMPARISONS,
            extra_config=list(extra), prepare=setup)

    def test_majiq_licence_is_checked_before_any_job_and_can_be_named(self):
        missing = self.module_dry_run("quant_majiq", [
            "umbrella_modules_mode: ingest", "majiq_bin_folder: /opt/majiq/bin",
            "majiq_license: licences/none.lic"], majiq_license=False)
        self.assertNotEqual(missing.returncode, 0, missing.stdout)
        self.assertIn("MAJIQ v3 needs its licence file", missing.stdout)
        self.assertIn("licences/none.lic does not exist", missing.stdout)

        def licence(temp_dir):
            os.makedirs(os.path.join(temp_dir, "licences"))
            open(os.path.join(temp_dir, "licences", "academic.lic"), "w").close()

        named = self.module_dry_run("quant_majiq", [
            "umbrella_modules_mode: ingest", "majiq_bin_folder: /opt/majiq/bin",
            "majiq_license: licences/academic.lic"], prepare=licence, majiq_license=False)
        self.assertEqual(named.returncode, 0, named.stdout)
        exported = re.findall(r"export MAJIQ_LICENSE_FILE=(\S+)/licences/academic.lic; ", named.stdout)
        self.assertGreaterEqual(len(exported), 3, named.stdout)   # reference, sj and compare
        # QAPA alone never asks for it
        qapa = self.module_dry_run("quant_qapa", ["umbrella_modules_mode: ingest"], majiq_license=False)
        self.assertEqual(qapa.returncode, 0, qapa.stdout)

    def test_qapa_alone_needs_no_alignment_whippet_or_full_reference(self):
        result = self.module_dry_run("quant_qapa", ["umbrella_modules_mode: ingest"])
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        for expected in ("rule umbrella_qapa_quant:", "rule umbrella_qapa_index:", "qapa build -N",
                         "rule umbrella_qapa_usage:", "src/qapa_dexseq.R"):
            self.assertIn(expected, out)
        for absent in ("rule umbrella_hisat2:", "whippet-quant.jl", "rule umbrella_salmon_quant:",
                       "rule umbrella_reference_manifest:", "rule umbrella_whippet_index:",
                       "rule umbrella_salmon_index:", "Round2", "rule umbrella_dapars2", "rule umbrella_majiq"):
            self.assertNotIn(absent, out)

    def test_apa_references_can_use_their_own_annotation_and_polya_sites(self):
        def apa_inputs(temp_dir):
            os.makedirs(os.path.join(temp_dir, "apa"), exist_ok=True)
            for name in ("gencode.gtf.gz", "polyAs.gtf.gz", "polyasite.bed.gz"):
                open(os.path.join(temp_dir, "apa", name), "w").close()

        result = self.module_dry_run("quant_umbrella_modules", [
            "umbrella_modules: [dapars2, qapa]", "umbrella_modules_mode: ingest",
            "apa_annotation_gtf: apa/gencode.gtf.gz", "qapa_gencode_polya: apa/polyAs.gtf.gz",
            "qapa_polyasite: apa/polyasite.bed.gz"], prepare=apa_inputs)
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        for expected in ("apa_annotation.py dapars-utr --gtf apa/gencode.gtf.gz",
                         "apa_annotation.py qapa-db --gtf apa/gencode.gtf.gz",
                         "apa_annotation.py qapa-gtf --gtf apa/gencode.gtf.gz",
                         "--transcripts basic", "apa_annotation.py polya-bed --gtf apa/polyAs.gtf.gz",
                         "/gencode_polya_sites.bed -p apa/polyasite.bed.gz", "check-bed"):
            self.assertIn(expected, out)
        self.assertNotIn("qapa build -N", out)
        # the splicing reference (and its identity) is untouched
        self.assertIn("rule umbrella_hisat2:", out)

    def test_megasearch_apa_set_feeds_both_apa_references(self):
        def apa_inputs(temp_dir):
            os.makedirs(os.path.join(temp_dir, "apa"), exist_ok=True)
            for name in ("gencode.gtf.gz", "polyAs.gtf.gz", "polyasite.bed.gz"):
                open(os.path.join(temp_dir, "apa", name), "w").close()

        result = self.module_dry_run("quant_umbrella_modules", [
            "umbrella_modules: [dapars2, qapa]", "umbrella_modules_mode: ingest",
            "apa_annotation_gtf: apa/gencode.gtf.gz", "apa_transcripts: megasearch",
            "qapa_gencode_polya: apa/polyAs.gtf.gz", "qapa_polyasite: apa/polyasite.bed.gz"], prepare=apa_inputs)
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        self.assertEqual(out.count("rule umbrella_apa_set:"), 1)
        self.assertIn("apa_annotation.py apa-set --gtf apa/gencode.gtf.gz --gencode-polya apa/polyAs.gtf.gz "
                      "--polyasite apa/polyasite.bed.gz --slop 50", out)
        set_gtf = re.search(r"--out (umbrella/modules/reference/[^ ]+/apa_set/[0-9a-f]{16}/apa_set.gtf.gz)", out)
        self.assertIsNotNone(set_gtf, out)
        for command in ("dapars-utr", "qapa-db", "qapa-gtf"):
            self.assertIn("apa_annotation.py {} --gtf {}".format(command, set_gtf.group(1)), out)
        self.assertNotIn("--transcripts basic", out)
        self.assertIn("--transcripts all", out)

    def test_megasearch_apa_set_needs_poly_a_evidence(self):
        result = self.module_dry_run("quant_dapars2", ["umbrella_modules_mode: ingest",
                                                       "apa_transcripts: megasearch"])
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("apa_transcripts: megasearch selects 3' ends with poly(A) evidence", result.stdout)

    def test_qapa_polya_tracks_go_together(self):
        result = self.module_dry_run("quant_qapa", [
            "umbrella_modules_mode: ingest", "qapa_gencode_polya: apa/polyAs.gtf.gz"])
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("qapa_gencode_polya and qapa_polyasite go together", result.stdout)

    def test_dapars2_alone_aligns_without_salmon_or_whippet(self):
        result = self.module_dry_run("quant_dapars2", ["umbrella_modules_mode: ingest"])
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        for expected in ("rule umbrella_hisat2:", "run_dapars2.py coverage",
                         "rule umbrella_dapars2_software:", "run_dapars2.py compare"):
            self.assertIn(expected, out)
        for absent in ("rule umbrella_salmon_quant:", "whippet-quant.jl", "rule umbrella_salmon_index:",
                       "rule umbrella_whippet_index:", "rule umbrella_qapa", "rule umbrella_featurecounts:",
                       "rule umbrella_junctions:"):
            self.assertNotIn(absent, out)

    def test_majiq_needs_its_licensed_install(self):
        missing = self.module_dry_run("quant_majiq", ["umbrella_modules_mode: ingest"])
        self.assertNotEqual(missing.returncode, 0, missing.stdout)
        self.assertIn("MAJIQ v3 is licensed", missing.stdout)
        installed = self.module_dry_run("quant_majiq", ["umbrella_modules_mode: ingest",
                                                        "majiq_bin_folder: /opt/majiq/bin"])
        self.assertEqual(installed.returncode, 0, installed.stdout)
        self.assertIn("run_majiq.py sj", installed.stdout)
        self.assertIn("/opt/majiq/bin/majiq-build gff3", installed.stdout)
        self.assertNotIn("rule umbrella_salmon_quant:", installed.stdout)
        built = self.module_dry_run("quant_majiq", ["umbrella_modules_mode: ingest",
                                                    "majiq_source: /licensed/majiq-v3.tar.gz"])
        self.assertNotEqual(built.returncode, 0)   # the source path does not exist here
        self.assertIn("majiq-v3.tar.gz", built.stdout)

    def test_all_modules_share_one_alignment_per_run(self):
        result = self.module_dry_run("quant_umbrella_modules", [
            "umbrella_modules: [majiq, dapars2, qapa]", "umbrella_modules_mode: ingest",
            "majiq_bin_folder: /opt/majiq/bin"])
        self.assertEqual(result.returncode, 0, result.stdout)
        stats = dict(re.findall(r"^(umbrella_\w+)\s+(\d+)$", result.stdout, re.M))
        self.assertEqual(stats.get("umbrella_hisat2"), "4")
        for rule in ("umbrella_majiq_sj", "umbrella_dapars2_coverage", "umbrella_qapa_quant"):
            self.assertEqual(stats.get(rule), "4", rule)
        for rule in ("umbrella_majiq_compare", "umbrella_dapars2_compare", "umbrella_qapa_usage"):
            self.assertEqual(stats.get(rule), "1", rule)
        self.assertIn("umbrella_module_inventory", result.stdout)
        self.assertNotIn("umbrella_whippet_quant", result.stdout)

    def test_module_selection_errors_are_actionable(self):
        for extra, message in (
                (["umbrella_modules: [majiq, tapas]"], "unknown umbrella_modules tapas"),
                (["umbrella_modules: [qapa]", "umbrella_allow_restage: true"], "cached_only never stages"),
                (["umbrella_modules: [qapa]"], "cached_only, but 4 run/tool caches are missing")):
            result = self.module_dry_run("quant_qapa", extra)
            self.assertNotEqual(result.returncode, 0, extra)
            self.assertIn(message, result.stdout)
        empty = self.module_dry_run("quant_umbrella_modules")
        self.assertNotEqual(empty.returncode, 0)
        self.assertIn("quant_umbrella_modules needs umbrella_modules", empty.stdout)

    @staticmethod
    def write_module_caches(temp_dir, tools, names=None, metadata=True):
        """Caches as an ingest leaves them: files plus Snakemake records holding
        the producer's ingest-time input set."""
        import base64
        import json as _json
        sys.path.insert(0, temp_dir)
        from src.umbrella_manifest import load_umbrella_manifest
        from src.umbrella_modules import CACHE_FILES, cache_id, module_reference_id
        producers = {"dapars2": "umbrella_dapars2_coverage", "qapa": "umbrella_qapa_quant",
                     "majiq": "umbrella_majiq_sj"}
        settings = {"majiq": {"majiq_bin_folder": "/opt/majiq/bin"}}
        for tool in tools:
            reference = module_reference_id("ref", tool, settings.get(tool, {}))
            for record in load_umbrella_manifest(os.path.join(temp_dir, "umbrella.tsv")).included_runs():
                work = "umbrella/work/ref/project/{}/{}".format(record.batch_id, record.run_id)
                folder = "umbrella/modules/cache/ref/project/{}/{}/{}/{}".format(
                    record.batch_id, record.run_id, tool, cache_id(record, reference))
                os.makedirs(os.path.join(temp_dir, folder))
                for name in (names or {}).get(tool, CACHE_FILES[tool]):
                    path = folder + "/" + name
                    open(os.path.join(temp_dir, path), "w").write(_json.dumps({}))
                    if metadata:
                        meta = os.path.join(temp_dir, ".snakemake", "metadata",
                                            base64.urlsafe_b64encode(path.encode()).decode())
                        os.makedirs(os.path.dirname(meta), exist_ok=True)
                        open(meta, "w").write(_json.dumps({
                            "rule": producers[tool], "incomplete": False, "starttime": 1.0,
                            "endtime": 2.0, "input": [work + "/align/aligned.bam", work + "/R1.fastq.gz"],
                            "params": ["x"], "code": "c", "shellcmd": "s"}))
        sys.path.remove(temp_dir)

    def test_named_targets_ignore_other_tools_caches(self):
        # GPT-6.1 review #5: a missing MAJIQ cache must not block QAPA targets
        result = self.module_dry_run(
            "quant_qapa", ["umbrella_modules: [majiq, qapa]", "majiq_bin_folder: /opt/majiq/bin"],
            prepare=lambda d: self.write_module_caches(d, ["qapa"]))
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("rule umbrella_qapa_pau:", result.stdout)
        self.assertNotIn("rule umbrella_majiq", result.stdout)
        # ...while the full module target still refuses the missing MAJIQ caches
        both = self.module_dry_run(
            "quant_umbrella_modules", ["umbrella_modules: [majiq, qapa]", "majiq_bin_folder: /opt/majiq/bin"],
            prepare=lambda d: self.write_module_caches(d, ["qapa"]))
        self.assertNotEqual(both.returncode, 0)
        self.assertIn("run/tool caches are missing (first: run", both.stdout)
        self.assertIn("majiq", both.stdout)

    def test_majiq_marker_without_its_sj_is_not_a_cache(self):
        result = self.module_dry_run(
            "prepare_majiq", ["umbrella_modules_mode: ingest", "majiq_bin_folder: /opt/majiq/bin"],
            prepare=lambda d: self.write_module_caches(d, ["majiq"], names={"majiq": ["cache.json"]}))
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("have their marker but not all their files", result.stdout)
        complete = self.module_dry_run(
            "prepare_majiq", ["majiq_bin_folder: /opt/majiq/bin"],
            prepare=lambda d: self.write_module_caches(d, ["majiq"]))
        self.assertEqual(complete.returncode, 0, complete.stdout)
        self.assertNotIn("rule umbrella_majiq_sj:", complete.stdout)

    def test_cached_only_comparisons_never_touch_reads(self):
        # caches with ingest-time Snakemake records (GPT-6.1 review, #1)
        def write_caches(temp_dir):
            self.write_module_caches(temp_dir, ["dapars2", "qapa"])
        for target, tools in (("quant_dapars2", "[dapars2]"), ("quant_qapa", "[qapa]"),
                              ("quant_umbrella_modules", "[dapars2, qapa]")):
            result = self.module_dry_run(target, ["umbrella_modules: " + tools], prepare=write_caches)
            self.assertEqual(result.returncode, 0, result.stdout)
            for absent in ("rule umbrella_stage", "rule umbrella_hisat2:", "rule umbrella_dapars2_coverage:",
                           "rule umbrella_qapa_quant:", "rule umbrella_validate_reads:",
                           "Set of input files has changed"):
                self.assertNotIn(absent, result.stdout, target)
        self.assertIn("rule umbrella_qapa_pau:", result.stdout)
        self.assertIn("rule umbrella_dapars2_compare:", result.stdout)

    def test_quant_umbrella_writes_the_module_inventory(self):
        def write_caches(temp_dir):
            self.write_module_caches(temp_dir, ["dapars2", "qapa"])
        result = self.module_dry_run("quant_umbrella", ["umbrella_modules: [dapars2, qapa]"],
                                     prepare=write_caches)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("rule umbrella_module_inventory:", result.stdout)
        self.assertIn("rule umbrella_qapa_pau:", result.stdout)
        without = self.module_dry_run("quant_umbrella")
        self.assertEqual(without.returncode, 0, without.stdout)
        self.assertNotIn("umbrella_module_inventory", without.stdout)

    def test_a_library_of_two_runs_is_processed_once(self):
        rows = ["sample_b\trun_b\trep_b\tproject\tbatch\tcontrol\tfastq\tr1.fastq.gz\tr2.fastq.gz\tPE\tunstranded\tref\ttrue\t",
                "GSM_c1\tSRR_c1\trep_c\tproject\tbatch\tcase\tfastq\tr1.fastq.gz\tr2.fastq.gz\tPE\tunstranded\tref\ttrue\tLIB_c",
                "GSM_c2\tSRR_c2\trep_c\tproject\tbatch\tcase\tfastq\tr1.fastq.gz\tr2.fastq.gz\tPE\tunstranded\tref\ttrue\tLIB_c",
                "sample_d\trun_d\trep_d\tproject\tbatch\tcase\tfastq\tr1.fastq.gz\tr2.fastq.gz\tPE\tunstranded\tref\ttrue\t"]
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target="quant_umbrella", umbrella_extra_rows=rows,
            umbrella_comparisons=self.MODULE_COMPARISONS, umbrella_extra_column="library_id")
        self.assertEqual(result.returncode, 0, result.stdout)
        out = result.stdout
        stats = dict(re.findall(r"^(\w+)\s+(\d+)$", out, re.M))
        # 4 units from 5 rows: the two runs of LIB_c are one sample everywhere
        for rule in ("umbrella_stage_reads", "umbrella_hisat2", "umbrella_whippet_quant",
                     "umbrella_salmon_quant", "umbrella_fastqc"):
            self.assertEqual(stats.get(rule), "4", rule)
        self.assertIn("umbrella/work/ref/project/batch/LIB_c/R1.fastq.gz", out)
        self.assertIn("FASTQ/LIB_c.fastq.gz", out)
        self.assertNotIn("SRR_c1/", out)
        self.assertNotIn("FASTQ/SRR_c1", out)

    def test_se_salmon_uses_single_read_argument(self):
        result = WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, umbrella_layout="SE", print_shell=True,
            target="quant_umbrella"
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("salmon quant", result.stdout)
        self.assertIn("-r umbrella/work/ref/project/batch/run_se/R1.fastq.gz", result.stdout)
        self.assertNotIn("-2 umbrella/work/ref/project/batch/run_se", result.stdout)


if __name__ == "__main__":
    unittest.main()
