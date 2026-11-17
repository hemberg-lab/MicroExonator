import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from src.whippet_core import core_quant, import_core

KEPT = ("psi.gz", "jnc.gz", "gene.tpm.gz", "isoform.tpm.gz", "map.gz")
HEADER = "Gene\tNode\tCoord\tStrand\tType\tPsi\tCI_Width\tCI_Lo,Hi\tTotal_Reads\tComplexity\tEntropy\n"


class CoreTests(unittest.TestCase):
    def produce(self, root, bs=True):
        source = root / "core"
        def native(command, check):
            self.assertEqual(command.count("--circ"), 1)
            prefix = command[command.index("-o") + 1]
            for suffix in KEPT:
                with gzip.open(prefix + "." + suffix, "wt") as stream:
                    stream.write(HEADER if suffix == "psi.gz" else "kept\n")
                    if suffix == "psi.gz" and bs:
                        stream.write("G\t3-1\tchr1:30-1\t+\tBS\t0.2\t0.1\t0,1\t10\tK1\t0\n")
        with patch("src.whippet_core.subprocess.run", native):
            record = core_quant(["one.fq", "two.fq"], "index.jls", source / "quant", "run1", "ref",
                                whippet_bin="/missing/bin", flags="--biascorrect --circ")
        return source, record

    def test_core_keeps_map_and_import_reuses_outputs_without_native_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, record = self.produce(root)
            self.assertTrue((source / "quant.map.gz").exists())
            self.assertTrue(record["circ"])
            before = (source / "quant.psi.gz").read_bytes()
            with patch("src.whippet_core.subprocess.run", side_effect=AssertionError("must not execute Whippet")):
                imported = import_core(source, root / "cache", "run1", {"reference_id": "ref", "cache_id": "c"})
            self.assertEqual(imported["bs_events"], 1)
            self.assertEqual((root / "cache/quant.psi.gz").read_bytes(), before)
            self.assertEqual((source / "quant.psi.gz").read_bytes(), before)

    def test_zero_bs_run_is_valid_and_legacy_or_changed_input_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, record = self.produce(root, bs=False)
            self.assertEqual(import_core(source, root / "zero", "run1", {"reference_id": "ref"})["bs_events"], 0)
            with self.assertRaisesRegex(ValueError, "run"):
                import_core(source, root / "wrong", "other", {"reference_id": "ref"})
            with self.assertRaisesRegex(ValueError, "reference"):
                import_core(source, root / "wrongref", "run1", {"reference_id": "other"})
            with self.assertRaisesRegex(ValueError, "processing"):
                import_core(source, root / "wrongdata", "run1", {"reference_id": "ref", "processing_sha256": "different"})
            (source / "quant.jnc.gz").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "checksum"):
                import_core(source, root / "changed", "run1", {"reference_id": "ref"})
            (source / "quant.provenance.json").unlink()
            with self.assertRaisesRegex(ValueError, "provenance"):
                import_core(source, root / "legacy", "run1", {"reference_id": "ref"})
