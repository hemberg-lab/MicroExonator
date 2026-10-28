"""Snakemake dry runs of the circRNA modules and of per-comparison module scope."""
import os
import gzip
import re
import unittest

from tests import test_workflow_dag as twd

ROWS = twd.UmbrellaNativeQuantificationTests.UMBRELLA_ROWS + [
    "sample_e\trun_e\trep_e\tproject\tbatch\tse_control\tfastq\tr1.fastq.gz\t\tSE\tunstranded\tref\ttrue",
    "sample_f\trun_f\trep_f\tproject\tbatch\tse_case\tfastq\tr1.fastq.gz\t\tSE\tunstranded\tref\ttrue",
    "sample_g\trun_g\trep_g\tproject\tbatch\tse_case\tfastq\tr1.fastq.gz\t\tSE\tunstranded\tref\ttrue",
    "sample_h\trun_h\trep_h\tproject\tbatch\tse_control\tfastq\tr1.fastq.gz\t\tSE\tunstranded\tref\ttrue"]
# main: case (c, d) vs control (pe, b); small: c vs b only; single: SE runs
COMPARISONS = ("comparisons:\n"
               "  - comparison_id: main\n    project_id: project\n    group_a: case\n    group_b: control\n"
               "  - comparison_id: small\n    project_id: project\n"
               "    a: {label: only_c, samples: [sample_c]}\n    b: {label: only_b, samples: [sample_b]}\n"
               "  - comparison_id: single\n    project_id: project\n    group_a: se_case\n    group_b: se_control\n")
HEADER = "comparison_id\ttool\trole\tlibrary_evidence\tdecided_by\tapproved_by\treason\n"
CATALOGUE = ["circrna_catalogue_sources:", "  - {path: circatlas.bed, name: circAtlas}",
             "  - {path: circbase_hg19.bed, format: circbase, assembly: hg19, role: annotate, name: circBase}",
             "circrna_liftover_chain: hg19ToHg38.over.chain.gz"]


def selection_row(comparison, tool, approved="guillermo", role="test"):
    return "\t".join([comparison, tool, role, "A: Ribo-Zero (GSM1); B: Ribo-Zero (GSM2)", "guillermo", approved,
                      "rRNA-depleted pilot"]) + "\n"


class CircrnaDagTests(unittest.TestCase):
    def test_differential_cached_only_uses_retained_junctions(self):
        def setup(folder):
            twd.UmbrellaNativeQuantificationTests.write_module_caches(folder, ["circrna_star"])
            for group in ('case', 'control'):
                path = os.path.join(folder, 'junctions/ref/project', group, 'batch.junctions.tsv.gz')
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with gzip.open(path, 'wt') as stream:
                    runs = 'run_c\trun_d' if group == 'case' else 'run_pe\trun_b'
                    stream.write('chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\tshort_anchor_total\t' + runs + '\n')
        result = self.dry_run('quant_dcc', ['umbrella_modules: [dcc]', 'circrna_differential: true',
                              'circrna_hisat_delta: true'], [selection_row('main', 'dcc')],
                              catalogue=False, prepare=setup)
        self.assertEqual(result.returncode, 0, result.stdout[-6000:])
        self.assertIn('rule umbrella_circrna_differential:', result.stdout)
        for forbidden in ('rule umbrella_stage', 'rule umbrella_hisat2:', 'rule umbrella_junctions:',
                          'rule umbrella_circrna_star:'):
            self.assertNotIn(forbidden, result.stdout)

    def test_differential_refuses_missing_junction_shards(self):
        result = self.dry_run('quant_dcc', ['umbrella_modules: [dcc]', 'circrna_differential: true',
                              'circrna_hisat_delta: true'], [selection_row('main', 'dcc')],
                              catalogue=False,
                              prepare=lambda d: twd.UmbrellaNativeQuantificationTests.write_module_caches(d, ['circrna_star']))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('requires retained HISAT2 junction shards', result.stdout)

    def dry_run(self, target, extra=(), rows=(), catalogue=True, prepare=None):
        def setup(temp_dir):
            for name in ("circatlas.bed", "circbase_hg19.bed", "hg19ToHg38.over.chain.gz"):
                open(os.path.join(temp_dir, name), "w").close()
            with open(os.path.join(temp_dir, "selection.tsv"), "w") as stream:
                stream.write(HEADER + "".join(rows))
            if prepare:
                prepare(temp_dir)
        config = list(extra) + ["umbrella_module_selection: selection.tsv"] + (CATALOGUE if catalogue else [])
        return twd.WorkflowSelectionTests.run_quant_dry_run(
            self, umbrella=True, print_shell=True, target=target, umbrella_extra_rows=ROWS,
            umbrella_comparisons=COMPARISONS, extra_config=config, prepare=setup)

    @staticmethod
    def stats(result):
        return dict(re.findall(r"^(\S+)\s+(\d+)$", result.stdout, flags=re.M))

    def test_circrna_modules_run_only_on_selected_comparisons_beside_all_scope_modules(self):
        rows = [selection_row("small", tool) for tool in ("ciriquant", "circexplorer2", "dcc")]
        result = self.dry_run("quant_umbrella_modules", [
            "umbrella_modules: [vasttools, psisigma, ciriquant, circexplorer2, dcc]",
            "umbrella_modules_mode: ingest"], rows, prepare=self.retained_shards)
        self.assertEqual(result.returncode, 0, result.stdout[-4000:])
        stats = self.stats(result)
        # all-scope modules: every included run (5 PE + 4 SE) and every comparison
        self.assertEqual(stats.get("umbrella_vasttools_align"), "8")
        self.assertEqual(stats.get("umbrella_vasttools_compare"), "3")
        self.assertEqual(stats.get("umbrella_psisigma_compare"), "3")
        # circRNA modules: only the runs of the selected comparison (run_b and run_c)
        for rule in ("umbrella_ciriquant_prep", "umbrella_circrna_star", "umbrella_circexplorer2_calls"):
            self.assertEqual(stats.get(rule), "2", rule)
        for rule in ("umbrella_ciriquant_compare", "umbrella_circexplorer2_support", "umbrella_dcc_compare",
                     "umbrella_circrna_catalogue", "umbrella_circrna_star_index", "umbrella_circrna_refflat",
                     "umbrella_ciriquant_index"):
            self.assertEqual(stats.get(rule), "1", rule)
        self.assertEqual(stats.get("umbrella_circrna_source"), "2")
        out = result.stdout
        self.assertNotIn("/run_d/ciriquant/", out)
        self.assertIn("circrna_catalogue.py lift", out)                   # circBase hg19 -> hg38
        self.assertIn("--chain hg19ToHg38.over.chain.gz", out)
        self.assertIn("run_ciriquant.py prep", out)
        self.assertIn("--layout PE", out)
        self.assertIn("--index-dir umbrella/modules/reference/ref/ciriquant/", out)   # one shared index
        self.assertIn("--format circbase", out)
        self.assertRegex(out, r"--annotate \S+/sources/1\.bed")
        self.assertIn('"role": "test"', out)                                # selection reaches the comparison
        # one HISAT2 alignment per run shared by PSI-Sigma and CIRIquant
        self.assertEqual(stats.get("umbrella_hisat2"), "8")

    def test_selected_scope_needs_a_selection(self):
        result = self.dry_run("quant_ciriquant", ["umbrella_modules_mode: ingest"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ciriquant runs only on selected comparisons and none is selected", result.stdout)

    def test_selection_rows_are_validated(self):
        for row, message in ((selection_row("nope", "ciriquant"), "comparison nope is not in umbrella_comparisons"),
                             (selection_row("small", "ciriquant", approved=""), "approved_by is empty"),
                             (selection_row("small", "vasttools"), "vasttools runs on every comparison"),
                             (selection_row("small", "ciriquant") * 2, "small is selected twice")):
            result = self.dry_run("quant_ciriquant", ["umbrella_modules_mode: ingest"], [row])
            self.assertNotEqual(result.returncode, 0, row)
            self.assertIn(message, result.stdout)

    def test_ciriquant_refuses_single_end_comparisons_but_star_tools_accept_them(self):
        refused = self.dry_run("quant_ciriquant", ["umbrella_modules_mode: ingest"],
                               [selection_row("single", "ciriquant")])
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn("ciriquant needs paired-end runs", refused.stdout)
        accepted = self.dry_run("quant_dcc", ["umbrella_modules_mode: ingest"], [selection_row("single", "dcc")])
        self.assertEqual(accepted.returncode, 0, accepted.stdout[-3000:])
        self.assertEqual(self.stats(accepted).get("umbrella_circrna_star"), "4")

    def test_ciriquant_needs_a_catalogue(self):
        result = self.dry_run("quant_ciriquant", ["umbrella_modules_mode: ingest"],
                              [selection_row("small", "ciriquant")], catalogue=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("set circrna_catalogue_sources", result.stdout)

    def test_selection_is_provenance_not_identity(self):
        # cache IDs are run identity + module_reference_id; each dry run here has its own temporary
        # folder (absolute FASTQ paths), so compare the module identity and the runs selected
        def identities(result):
            return set(re.findall(r'"module_reference_id": "([0-9a-f]{16})"', result.stdout))

        def runs(result):
            return set(re.findall(r"umbrella/modules/cache/ref/project/batch/(\w+)/ciriquant/", result.stdout))
        one = self.dry_run("prepare_ciriquant", ["umbrella_modules_mode: ingest"],
                           [selection_row("small", "ciriquant")])
        two = self.dry_run("prepare_ciriquant", ["umbrella_modules_mode: ingest"],
                           [selection_row("small", "ciriquant", role="polya_reference"),
                            selection_row("main", "ciriquant")])
        self.assertEqual(one.returncode, 0, one.stdout[-3000:])
        self.assertEqual(two.returncode, 0, two.stdout[-3000:])
        self.assertEqual(len(identities(one)), 1)
        self.assertEqual(identities(one), identities(two))                 # role and rows are not identity
        self.assertEqual(runs(one), {"run_b", "run_c"})
        self.assertEqual(runs(two), {"run_b", "run_c", "run_d", "run_pe"})  # main adds its runs

    def test_scope_all_runs_every_comparison(self):
        result = self.dry_run("quant_circexplorer2", ["umbrella_modules_mode: ingest",
                                                      "umbrella_module_scope: {circexplorer2: all}"])
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        stats = self.stats(result)
        self.assertEqual(stats.get("umbrella_circexplorer2_support"), "3")
        self.assertEqual(stats.get("umbrella_circrna_star"), "8")

    def test_cached_only_discovery_reads_kept_star_caches_without_reads(self):
        rows = [selection_row("small", "circexplorer2"), selection_row("small", "dcc")]
        result = self.dry_run("quant_circexplorer2", ["umbrella_modules: [circexplorer2, dcc]"], rows,
                              prepare=lambda d: twd.UmbrellaNativeQuantificationTests.write_module_caches(
                                  d, ["circrna_star"]))
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        out = result.stdout
        self.assertIn("rule umbrella_circexplorer2_calls:", out)       # derived from the kept STAR cache
        self.assertIn("rule umbrella_circexplorer2_support:", out)
        for absent in ("rule umbrella_stage", "rule umbrella_hisat2:", "rule umbrella_circrna_star:",
                       "rule umbrella_circrna_star_index:"):
            self.assertNotIn(absent, out)
        dcc = self.dry_run("quant_dcc", ["umbrella_modules: [circexplorer2, dcc]"], rows,
                           prepare=lambda d: twd.UmbrellaNativeQuantificationTests.write_module_caches(
                               d, ["circrna_star"]))
        self.assertEqual(dcc.returncode, 0, dcc.stdout[-3000:])
        self.assertIn("rule umbrella_dcc_compare:", dcc.stdout)
        self.assertNotIn("rule umbrella_stage", dcc.stdout)
        missing = self.dry_run("quant_dcc", ["umbrella_modules: [dcc]"], [selection_row("small", "dcc")])
        self.assertNotEqual(missing.returncode, 0)
        self.assertIn("run/tool caches are missing", missing.stdout)

    def test_unselected_circrna_modules_add_nothing(self):
        result = self.dry_run("quant_umbrella_modules", ["umbrella_modules: [vasttools]",
                                                         "umbrella_modules_mode: ingest"],
                              [selection_row("small", "ciriquant")])
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        for absent in ("ciriquant", "circrna_star", "circexplorer2", "dcc_compare", "circrna_catalogue"):
            self.assertNotIn("rule umbrella_" + absent, result.stdout)


    def test_whippet_circ_quantifies_selected_runs_beside_the_core_whippet(self):
        result = self.dry_run("quant_whippet_circ", ["umbrella_modules: [whippet_circ]", "umbrella_modules_mode: ingest",
                                                     "whippet_flags: --biascorrect"],
                              [selection_row("small", "whippet_circ")], catalogue=False)
        self.assertEqual(result.returncode, 0, result.stdout[-4000:])
        stats = self.stats(result)
        self.assertEqual(stats.get("umbrella_whippet_circ_quant"), "2")            # run_b and run_c
        self.assertEqual(stats.get("umbrella_whippet_circ_compare"), "1")
        out = result.stdout
        self.assertIn("run_whippet_circ.py quant", out)
        self.assertIn("--flags --biascorrect", out)
        self.assertRegex(out, r"--index /\S+/whippet\.jls ")                    # the existing index, absolute
        for absent in ("rule umbrella_whippet_index:", "rule umbrella_whippet_quant:", "rule umbrella_hisat2:"):
            self.assertNotIn(absent, out)                                           # core Whippet untouched

    def test_whippet_circ_cached_only_compares_kept_caches(self):
        def setup(folder):
            twd.UmbrellaNativeQuantificationTests.write_module_caches(folder, ["whippet_circ"])
        result = self.dry_run("quant_whippet_circ", ["umbrella_modules: [whippet_circ]"],
                              [selection_row("small", "whippet_circ")], catalogue=False, prepare=setup)
        self.assertEqual(result.returncode, 0, result.stdout[-4000:])
        self.assertIn("run_whippet_circ.py compare", result.stdout)
        self.assertNotIn("rule umbrella_stage", result.stdout)
        self.assertNotIn("rule umbrella_whippet_circ_quant:", result.stdout)


    @staticmethod
    def retained_shards(folder):
        for kind, header in (("junctions", "chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\t"
                                             "short_anchor_total\t"), ("genes", "gene_id\tlength\t")):
            for group, runs in (("case", "run_c\trun_d"), ("control", "run_pe\trun_b")):
                suffix = ".junctions.tsv.gz" if kind == "junctions" else ".featurecounts.tsv.gz"
                path = os.path.join(folder, kind, "ref/project", group, "batch" + suffix)
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with gzip.open(path, "wt") as stream:
                    stream.write(header + runs + "\n")

    def test_native_tests_first_and_edger_only_where_there_is_none(self):
        rows = [selection_row("main", tool) for tool in ("ciriquant", "circexplorer2", "dcc")]
        result = self.dry_run("quant_umbrella_modules", ["umbrella_modules: [ciriquant, circexplorer2, dcc]",
                                                         "umbrella_modules_mode: ingest", "circrna_differential: true",
                                                         "circrna_hisat_delta: true"], rows,
                              prepare=self.retained_shards)
        self.assertEqual(result.returncode, 0, result.stdout[-4000:])
        stats = self.stats(result)
        for rule in ("umbrella_ciriquant_de", "umbrella_dcc_linear", "umbrella_dcc_circtest",
                     "umbrella_circrna_native_software"):
            self.assertEqual(stats.get(rule), "1", rule)
        self.assertEqual(stats.get("umbrella_circrna_differential"), "1")        # CIRCexplorer2 only
        out = result.stdout
        self.assertRegex(out, r"circrna_differential.*\n.*\n.*modules/circexplorer2/")
        self.assertRegex(out, r"--gene-shard /\S+/genes/ref/project/case/batch\.featurecounts\.tsv\.gz")
        self.assertIn("--bam run_b=umbrella/work/ref/project/batch/run_b/align/aligned.bam", out)
        self.assertIn("CIRI_DE.R", out)

    def test_edger_beside_native_on_request_and_dcc_falls_back_in_cached_only(self):
        rows = [selection_row("main", "dcc")]
        def setup(folder):
            self.retained_shards(folder)
            twd.UmbrellaNativeQuantificationTests.write_module_caches(folder, ["circrna_star"])
        cached = self.dry_run("quant_dcc", ["umbrella_modules: [dcc]", "circrna_differential: true",
                                            "circrna_hisat_delta: true"], rows, catalogue=False, prepare=setup)
        self.assertEqual(cached.returncode, 0, cached.stdout[-4000:])
        self.assertIn("rule umbrella_circrna_differential:", cached.stdout)     # no BAMs: edgeR backup
        self.assertNotIn("rule umbrella_dcc_circtest:", cached.stdout)
        both = self.dry_run("quant_dcc", ["umbrella_modules: [dcc]", "umbrella_modules_mode: ingest",
                                          "circrna_differential: true", "circrna_hisat_delta: true",
                                          "circrna_differential_also_native: true"], rows, catalogue=False,
                            prepare=self.retained_shards)
        self.assertEqual(both.returncode, 0, both.stdout[-4000:])
        self.assertIn("rule umbrella_dcc_circtest:", both.stdout)
        self.assertIn("rule umbrella_circrna_differential:", both.stdout)


if __name__ == "__main__":
    unittest.main()
