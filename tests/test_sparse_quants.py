"""The corrected PSI files keep microexons with at least min_reads_PSI corrected reads."""

import gzip
import pathlib
import runpy
import sys
import tempfile
import types
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]


def run_script(name, inputs, output, min_reads, config=None):
    fake = types.SimpleNamespace(input=inputs, output={"corrected_sparse": str(output)},
                                 params={"min_reads": min_reads}, config=config or {})
    # the scripts import snakemake.utils, which is not needed here
    stub = types.ModuleType("snakemake")
    stub.utils = types.ModuleType("snakemake.utils")
    stub.utils.min_version = lambda version: None
    saved = {key: sys.modules.get(key) for key in ("snakemake", "snakemake.utils")}
    sys.modules.update({"snakemake": stub, "snakemake.utils": stub.utils})
    try:
        runpy.run_path(str(REPOSITORY / "src" / name), init_globals={"snakemake": fake})
    finally:
        for key, module in saved.items():
            if module is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = module
    with gzip.open(output, "rt") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle]
    key = rows[0].index("ME")
    return rows[0], {row[key]: row for row in rows[1:]}


def write_counts(path, rows):
    with gzip.open(path, "wt") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


class SingleEndSparseQuantTests(unittest.TestCase):
    ROWS = [("s1", "ME_low", 2.0, 1.0), ("s1", "ME_mid", 5.0, 2.0), ("s1", "ME_high", 9.0, 3.0),
            ("s1", "ME_none", 0.0, 0.0)]

    def run_se(self, min_reads):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            write_counts(tmp / "s1.ME.adj_counts.gz", self.ROWS)
            return run_script("get_sparse_quants_se.py",
                              {"corrected_quant": str(tmp / "s1.ME.adj_counts.gz")},
                              tmp / "s1.corrected.PSI.gz", min_reads)

    def test_default_cut_off_is_five_reads(self):
        header, rows = self.run_se("5")
        self.assertEqual(header, ["sample", "ME", "ME_coverages", "excluding_covs", "PSI", "CI_Lo", "CI_Hi"])
        self.assertEqual(sorted(rows), ["ME_high", "ME_mid"])
        self.assertAlmostEqual(float(rows["ME_high"][4]), 0.75)

    def test_cut_off_follows_min_reads_psi(self):
        self.assertEqual(sorted(self.run_se("10")[1]), ["ME_high"])
        self.assertEqual(sorted(self.run_se("1")[1]), ["ME_high", "ME_low", "ME_mid"])

    def test_zero_cut_off_skips_events_without_reads(self):
        self.assertNotIn("ME_none", self.run_se("0")[1])


class PairedEndSparseQuantTests(unittest.TestCase):
    def test_mates_are_summed_before_the_cut_off(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            (tmp / "pairs.tsv").write_text("s_1\ts_2\n")
            write_counts(tmp / "s_1.gz", [("s_1", "ME_a", 3.0, 1.0), ("s_1", "ME_b", 1.0, 0.0)])
            write_counts(tmp / "s_2.gz", [("s_2", "ME_a", 2.0, 1.0), ("s_2", "ME_b", 1.0, 0.0)])
            inputs = {"corrected_quant_rd1": str(tmp / "s_1.gz"), "corrected_quant_rd2": [str(tmp / "s_2.gz")]}
            config = {"paired_samples": str(tmp / "pairs.tsv")}
            _, rows = run_script("get_sparse_quants_pe.py", inputs, tmp / "s_1.corrected.PSI.gz", "5", config)
            self.assertEqual(sorted(rows), ["ME_a"])
            self.assertEqual(float(rows["ME_a"][2]) + float(rows["ME_a"][3]), 7.0)
            _, rows = run_script("get_sparse_quants_pe.py", inputs, tmp / "s_1.corrected.PSI.gz", "2", config)
            self.assertEqual(sorted(rows), ["ME_a", "ME_b"])


class PseudoBulkSparseQuantTests(unittest.TestCase):
    def run_sp(self, min_reads, pool="E8.5_NMP-2"):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            write_counts(tmp / "c1.gz", [("c1", "ME_on", 3.0, 0.0), ("c1", "ME_skipped", 0.0, 4.0),
                                         ("c1", "ME_rare", 1.0, 0.0)])
            write_counts(tmp / "c2.gz", [("c2", "ME_on", 2.0, 1.0), ("c2", "ME_skipped", 0.0, 3.0),
                                         ("c2", "ME_none", 0.0, 0.0)])
            output = tmp / "single_cell" / "{}.corrected.PSI.gz".format(pool)
            output.parent.mkdir()
            return run_script("get_sparse_quants_sp.py", {"cells": [str(tmp / "c1.gz"), str(tmp / "c2.gz")]},
                              output, min_reads)

    def test_cells_are_summed_and_fully_skipped_microexons_are_kept(self):
        header, rows = self.run_sp("5")
        self.assertEqual(header, ["ME", "pseudo_pool", "cell_type", "ME_coverages", "excluding_covs",
                                  "PSI", "CI_Lo", "CI_Hi"])
        self.assertEqual(sorted(rows), ["ME_on", "ME_skipped"])
        self.assertAlmostEqual(float(rows["ME_on"][5]), 5 / 6)
        # seven exclusion reads and no inclusion: PSI 0, counted by the confidence filter
        self.assertEqual(float(rows["ME_skipped"][5]), 0.0)

    def test_cut_off_follows_min_reads_psi(self):
        self.assertEqual(sorted(self.run_sp("7")[1]), ["ME_skipped"])
        self.assertEqual(sorted(self.run_sp("1")[1]), ["ME_on", "ME_rare", "ME_skipped"])
        self.assertNotIn("ME_none", self.run_sp("0")[1])

    def test_rows_name_the_pseudo_bulk_and_its_cell_type(self):
        _, rows = self.run_sp("5")
        self.assertEqual(rows["ME_on"][1:3], ["E8.5_NMP-2", "E8.5 NMP"])


if __name__ == "__main__":
    unittest.main()
