import tempfile
import unittest
from pathlib import Path

from src.shard_guard import (hisat2_members, reference_manifest,
                             checksum_path, validate_fixed_microexons,
                             write_immutable_bundle)


class ShardGuardTests(unittest.TestCase):
    def test_hisat2_declares_all_small_and_large_members(self):
        self.assertEqual(hisat2_members("idx/genome"), [
            "idx/genome.1.ht2", "idx/genome.2.ht2", "idx/genome.3.ht2",
            "idx/genome.4.ht2", "idx/genome.5.ht2", "idx/genome.6.ht2",
            "idx/genome.7.ht2", "idx/genome.8.ht2",
        ])
        self.assertEqual(hisat2_members("idx/genome", large=True), [
            "idx/genome.1.ht2l", "idx/genome.2.ht2l", "idx/genome.3.ht2l",
            "idx/genome.4.ht2l", "idx/genome.5.ht2l", "idx/genome.6.ht2l",
            "idx/genome.7.ht2l", "idx/genome.8.ht2l",
        ])

    def test_reference_identity_includes_complete_index_and_rejects_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = {}
            for name in ("genome", "annotation", "whippet", "salmon", "tx2gene"):
                path = root / name
                path.write_text(name)
                paths[name] = str(path)
            paths["hisat2"] = []
            for member in hisat2_members(str(root / "index")):
                Path(member).write_text(member)
                paths["hisat2"].append(member)
            original = reference_manifest(paths)
            self.assertEqual(len(original["checksums"]["hisat2"]), 8)
            Path(paths["hisat2"][7]).write_text("changed")
            changed = reference_manifest(paths)
            self.assertNotEqual(original["reference_id"], changed["reference_id"])
            with self.assertRaisesRegex(ValueError, "reference_id mismatch"):
                reference_manifest(paths, expected_id=original["reference_id"])

    def test_directory_checksum_ignores_snakemake_timestamp_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "hash.bin").write_bytes(b"index")
            (root / ".snakemake_timestamp").write_text("first")
            original = checksum_path(root)
            (root / ".snakemake_timestamp").write_text("second")
            self.assertEqual(checksum_path(root), original)
            (root / "hash.bin").write_bytes(b"changed")
            self.assertNotEqual(checksum_path(root), original)

    def test_fixed_whippet_annotation_must_contain_me_universe(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            me_db = root / "me.txt"
            gtf = root / "fixed.gtf"
            me_db.write_text("chr1_+_10_15\n")
            gtf.write_text('chr1\tsource\texon\t11\t15\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
            self.assertEqual(validate_fixed_microexons(gtf, me_db), 1)
            gtf.write_text('chr1\tsource\texon\t12\t15\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
            with self.assertRaisesRegex(ValueError, "fixed microexon universe"):
                validate_fixed_microexons(gtf, me_db)

    def test_immutable_bundle_accepts_identical_and_refuses_changed_manifest_or_data(self):
        with tempfile.TemporaryDirectory() as directory:
            shard = Path(directory) / "counts.tsv.gz"
            guard = Path(directory) / "checksums.json"
            outputs = {str(shard): b"same bytes"}
            metadata = {"manifest_sha256": "first", "reference_id": "abc"}
            write_immutable_bundle(outputs, guard, metadata)
            first = guard.read_bytes()
            write_immutable_bundle(outputs, guard, metadata)
            self.assertEqual(guard.read_bytes(), first)
            with self.assertRaisesRegex(ValueError, "immutable"):
                write_immutable_bundle(outputs, guard, {**metadata, "manifest_sha256": "second"})
            with self.assertRaisesRegex(ValueError, "immutable"):
                write_immutable_bundle({str(shard): b"different"}, guard, metadata)
            self.assertEqual(shard.read_bytes(), b"same bytes")


if __name__ == "__main__":
    unittest.main()
