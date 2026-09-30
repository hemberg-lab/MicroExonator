"""Module selection, identities, 3' UTR resources, and the QAPA/MAJIQ wrappers (mock tools)."""

import gzip
import json
import os
import stat
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src import apa_annotation, run_majiq, run_qapa
from src.umbrella_manifest import load_umbrella_manifest
from src.umbrella_modules import (analysis_id, cache_id, missing_caches, module_reference_id,
                                  parse_selection)

HEADER = ("sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\tsource_type\t"
          "source_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\tsex\n")


def manifest(root, group="ctrl", sex="F", batch="b1"):
    path = root / "m.tsv"
    path.write_text(HEADER + "s1\tSRR1\tr1\tproj\t{}\t{}\tsra\tSRR1\t\tPE\tauto\tref\ttrue\t{}\n".format(
        batch, group, sex))
    return load_umbrella_manifest(path).by_run["SRR1"]


class SelectionTests(unittest.TestCase):
    def test_selection_and_mode(self):
        self.assertEqual(parse_selection({}), ([], "cached_only"))
        self.assertEqual(parse_selection({"umbrella_modules": ["qapa", "majiq"],
                                          "umbrella_modules_mode": "ingest",
                                          "umbrella_allow_restage": True}), (["majiq", "qapa"], "ingest"))
        self.assertEqual(parse_selection({"umbrella_modules": "dapars2, qapa"})[0], ["dapars2", "qapa"])
        for bad, message in (({"umbrella_modules": ["tapas"]}, "unknown"),
                             ({"umbrella_modules": ["qapa", "qapa"]}, "twice"),
                             ({"umbrella_modules_mode": "fast"}, "cached_only or ingest"),
                             ({"umbrella_allow_restage": "true"}, "never stages")):
            with self.assertRaisesRegex(ValueError, message):
                parse_selection(bad)


class IdentityTests(unittest.TestCase):
    def test_cache_ids_ignore_labels_and_follow_processing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = module_reference_id("ref", "qapa", {})
            base = cache_id(manifest(root), reference)
            self.assertEqual(cache_id(manifest(root, group="case", sex="M"), reference), base)
            self.assertNotEqual(cache_id(manifest(root, batch="b2"), reference), base)
            self.assertNotEqual(cache_id(manifest(root), module_reference_id("ref", "dapars2", {})), base)

    def test_module_reference_ids_are_per_tool(self):
        qapa = module_reference_id("ref", "qapa", {})
        self.assertNotEqual(module_reference_id("ref", "qapa", {"qapa_polya_sites": "sites.bed"}), qapa)
        self.assertEqual(module_reference_id("ref", "dapars2", {"qapa_polya_sites": "sites.bed"}),
                         module_reference_id("ref", "dapars2", {}))

    def test_analysis_id_follows_the_selected_runs(self):
        preflight = {"runs": {"a": ["r1", "r2"], "b": ["r3", "r4"]},
                     "collapse": {"r1": "r1", "r2": "r2", "r3": "r3", "r4": "r4"}}
        other = {"runs": {"a": ["r1", "r2"], "b": ["r3"]}, "collapse": {"r1": "r1", "r2": "r2", "r3": "r3"}}
        self.assertNotEqual(analysis_id(preflight, "m", {}), analysis_id(other, "m", {}))
        self.assertEqual(analysis_id(preflight, "m", {}), analysis_id(dict(preflight), "m", {}))

    def test_missing_caches(self):
        with tempfile.TemporaryDirectory() as directory:
            run = manifest(Path(directory))
            folder = Path(directory) / "cache"
            folder.mkdir()
            (folder / "quant.sf.gz").write_text("")
            missing = missing_caches([run], lambda record, tool: str(folder), ["qapa"])
            self.assertEqual([(run_id, tool) for run_id, tool, _ in missing], [("SRR1", "qapa")])


GTF = "\n".join([
    'chr1\tx\ttranscript\t101\t1000\t.\t+\t.\tgene_id "ENSG1.4"; transcript_id "ENST1.2"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "AAA";',
    'chr1\tx\texon\t101\t200\t.\t+\t.\tgene_id "ENSG1.4"; transcript_id "ENST1.2"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "AAA";',
    'chr1\tx\texon\t501\t1000\t.\t+\t.\tgene_id "ENSG1.4"; transcript_id "ENST1.2"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "AAA";',
    'chr1\tx\tCDS\t151\t200\t.\t+\t0\tgene_id "ENSG1.4"; transcript_id "ENST1.2"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "AAA";',
    'chr1\tx\tCDS\t501\t600\t.\t+\t0\tgene_id "ENSG1.4"; transcript_id "ENST1.2"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "AAA";',
    # minus strand, coding end in the last (leftmost) exon
    'chr2\tx\texon\t1001\t1500\t.\t-\t.\tgene_id "ENSG2.1"; transcript_id "ENST2.1"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "BBB";',
    'chr2\tx\texon\t2001\t2100\t.\t-\t.\tgene_id "ENSG2.1"; transcript_id "ENST2.1"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "BBB";',
    'chr2\tx\tCDS\t1301\t1500\t.\t-\t0\tgene_id "ENSG2.1"; transcript_id "ENST2.1"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "BBB";',
    # 3' UTR spanning an intron: coding ends before the last exon
    'chr3\tx\texon\t1\t100\t.\t+\t.\tgene_id "ENSG3.1"; transcript_id "ENST3.1"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "CCC";',
    'chr3\tx\texon\t301\t400\t.\t+\t.\tgene_id "ENSG3.1"; transcript_id "ENST3.1"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "CCC";',
    'chr3\tx\tCDS\t11\t90\t.\t+\t0\tgene_id "ENSG3.1"; transcript_id "ENST3.1"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "CCC";',
    # non-coding, and a synthetic model without types
    'chr4\tx\texon\t1\t100\t.\t+\t.\tgene_id "ENSG4.1"; transcript_id "ENST4.1"; gene_type "lncRNA"; transcript_type "lncRNA"; gene_name "DDD";',
    'chr5\tx\texon\t1\t100\t.\t+\t.\tgene_id "ME1"; transcript_id "ME1.t";',
]) + "\n"


class AnnotationTests(unittest.TestCase):
    def test_qapa_db_and_dapars_utrs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with gzip.open(root / "a.gtf.gz", "wt") as stream:
                stream.write(GTF)
            transcripts = apa_annotation.read_transcripts(root / "a.gtf.gz")
            apa_annotation.qapa_db(transcripts, root / "db.txt")
            lines = (root / "db.txt").read_text().splitlines()
            self.assertEqual(lines[1], "ENSG1\tENST1\tprotein_coding\tprotein_coding\tAAA")
            self.assertFalse(any("ME1" in line for line in lines))
            kept, skipped = apa_annotation.dapars_utr(transcripts, root / "utr.bed", root / "win.bed")
            self.assertEqual((root / "utr.bed").read_text().splitlines(), [
                "chr1\t600\t1000\tENST1|AAA|chr1|+\t0\t+", "chr2\t1000\t1300\tENST2|BBB|chr2|-\t0\t-"])
            self.assertEqual(skipped, {"utr_with_intron": 1, "not_protein_coding": 2})
            self.assertEqual((root / "win.bed").read_text().splitlines(), ["chr1\t600\t1000", "chr2\t1000\t1300"])


def write_quant(path, rows):
    with gzip.open(path, "wt") as stream:
        stream.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        for name, effective, reads in rows:
            stream.write("{}\t500\t{}\t0\t{}\n".format(name, effective, reads))


NAME_P = "ENST1_ENSG1_hsa_chr1_100_1000_+_utr_600_800"
NAME_D = "ENST9_ENSG1_hsa_chr1_100_1500_+_utr_600_1300"


class QapaTests(unittest.TestCase):
    def test_technical_runs_combine_and_sites_are_parsed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_quant(root / "a.gz", [(NAME_P, 400.0, 10), (NAME_D, 900.0, 0)])
            write_quant(root / "b.gz", [(NAME_P, 300.0, 30), (NAME_D, 800.0, 20)])
            rows = run_qapa.combine_runs([root / "a.gz", root / "b.gz"])
            proximal, distal = rows
            self.assertEqual(proximal[4], 40)
            self.assertAlmostEqual(proximal[2], 325.0)          # read-weighted
            self.assertAlmostEqual(distal[2], 800.0)
            self.assertAlmostEqual(proximal[3] + distal[3], 1e6)
            self.assertEqual(run_qapa.site_key(NAME_D), ("ENSG1", "600", "1300"))

    def test_pau_and_normalize_with_a_mock_qapa(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            caches = {}
            for run, reads in (("r1", (10, 30)), ("r2", (12, 28)), ("r3", (30, 10)), ("r4", (28, 12))):
                (root / run).mkdir()
                write_quant(root / run / "quant.sf.gz", [(NAME_P, 400.0, reads[0]), (NAME_D, 900.0, reads[1])])
                caches[run] = str(root / run)
            preflight = root / "preflight.json"
            preflight.write_text(json.dumps({
                "project_id": "p", "comparison_id": "c", "inference_supported": True, "reasons": [],
                "replicates": {"a": ["r1", "r2"], "b": ["r3", "r4"]},
                "collapse": {"r1": "r1", "r2": "r2", "r3": "r3", "r4": "r4"}}))
            pau_header = ("APA_ID\tTranscript\tGene\tGene_Name\tChr\tLastExon.Start\tLastExon.End\tStrand\t"
                          "UTR3.Start\tUTR3.End\tLength\tNum_Events\tr1.PAU\tr2.PAU\tr3.PAU\tr4.PAU\n")

            def fake_qapa(command, stdout, stderr, check, cwd):
                self.assertEqual(command[:4], ["qapa", "quant", "--db", "db.txt"])
                self.assertEqual([Path(p).parent.name for p in command[4:]], ["r1", "r2", "r3", "r4"])
                stdout.write(pau_header + "ENSG1_1_P\tENST1\tENSG1\tAAA\tchr1\t100\t1000\t+\t600\t800\t200\t2\t25\t30\t75\t70\n")

            with patch("src.run_qapa.subprocess.run", side_effect=fake_qapa):
                run_qapa.pau(preflight, caches, "db.txt", root / "out")
            counts = gzip.open(root / "out" / "site_counts.tsv.gz", "rt").read().splitlines()
            self.assertEqual(counts[1], "ENSG1_600_1300\tENSG1\t30\t28\t10\t12")
            with gzip.open(root / "out" / "dexseq.tsv.gz", "wt") as stream:
                stream.write("site_id\tgene_id\tlog2fold_a_b\tpvalue\tpadj\nENSG1_600_800\tENSG1\t-1.5\t0.001\t0.01\n")
            run_qapa.normalize(preflight, root / "out" / "pau.tsv", root / "out" / "dexseq.tsv.gz", root / "out")
            rows = gzip.open(root / "out" / "results.tsv.gz", "rt").read().splitlines()
            header, row = rows[0].split("\t"), rows[1].split("\t")
            record = dict(zip(header, row))
            self.assertEqual(record["effect"], "-45")           # 27.5 - 72.5
            self.assertEqual((record["p_value"], record["q_value"], record["status"]), ("0.001", "0.01", "tested"))


FAKE_MAJIQ = """#!/bin/sh
case "$*" in *--help*) echo "usage"; exit 0;; esac
tool=$(basename "$0")
if [ "$tool" = "majiq-build" ] && [ "$1" = "update" ]; then mkdir -p "$3"; cp "$5" "$3/experiments.tsv"; exit 0; fi
if [ "$tool" = "majiq-build" ] && [ "$1" = "sj" ]; then echo sj > "$4"; exit 0; fi
if [ "$1" = "psi-coverage" ]; then shift 3; echo "$@" > "$(echo $0 | sed 's/majiq$//')last_psicov"; : ; fi
if [ "$1" = "psi-coverage" ] || [ "$tool" = "majiq" ] && [ "$1" = "psi-coverage" ]; then exit 0; fi
if [ "$1" = "heterogen" ]; then
  while [ "$1" != "--output-tsv" ]; do shift; done
  printf 'lsv_id\\tgene_id\\tgene_name\\tjunction_coord\\tmedian_dpsi\\ttnom_score\\nL1\\tENSG1\\tAAA\\t10-20\\t0.25\\t0\\n' > "$2"; exit 0
fi
"""


class MajiqTests(unittest.TestCase):
    def test_compare_runs_update_psicoverage_heterogen_and_normalizes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bin_dir = root / "bin"
            bin_dir.mkdir()
            for name in ("majiq", "majiq-build"):
                (bin_dir / name).write_text(FAKE_MAJIQ)
                (bin_dir / name).chmod(0o755)
            caches = {}
            for run in ("r1", "r1b", "r2", "r3", "r4"):
                (root / run).mkdir()
                (root / run / "run.sj").write_text("sj")
                caches[run] = str(root / run)
            preflight = root / "preflight.json"
            preflight.write_text(json.dumps({
                "project_id": "p", "comparison_id": "c", "inference_supported": True, "reasons": [],
                "replicates": {"a": ["A1", "A2"], "b": ["B1", "B2"]},
                "collapse": {"r1": "A1", "r1b": "A1", "r2": "A2", "r3": "B1", "r4": "B2"}}))
            run_majiq.compare(preflight, caches, root / "sg.zarr", root / "out", bin_dir=str(bin_dir))
            log = (root / "out" / "majiq.log").read_text()
            self.assertIn("--psi1", log)
            # technical runs of A1 are pooled into one PsiCoverage
            psicov = [line for line in log.splitlines() if "psi-coverage" in line and "A1.psicov" in line][0]
            self.assertIn("r1/run.sj", psicov)
            self.assertIn("r1b/run.sj", psicov)
            rows = gzip.open(root / "out" / "results.tsv.gz", "rt").read().splitlines()
            record = dict(zip(rows[0].split("\t"), rows[1].split("\t")))
            self.assertEqual((record["feature_id"], record["effect"], record["native_statistic"]), ("L1", "0.25", "0"))
            self.assertIn("not an FDR", record["statistic_type"])
            self.assertEqual(json.loads((root / "out" / "status.json").read_text())["status"], "ok")
            self.assertTrue((root / "out" / "splicegraph.tar.gz").exists())

    def test_unsupported_design_runs_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = root / "preflight.json"
            preflight.write_text(json.dumps({"project_id": "p", "comparison_id": "c",
                                             "inference_supported": False, "reasons": ["one replicate"]}))
            run_majiq.compare(preflight, {}, root / "sg", root / "out", bin_dir="/nonexistent")
            status = json.loads((root / "out" / "status.json").read_text())
            self.assertEqual((status["status"], status["reasons"]), ("unsupported", ["one replicate"]))


if __name__ == "__main__":
    unittest.main()
