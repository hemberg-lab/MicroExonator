import unittest.mock
import json
import io
import tempfile
import unittest
from pathlib import Path

from src.shard_guard import (hisat2_members, main as shard_guard_main, reference_identity, reference_manifest,
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

    def test_reference_id_needs_only_configured_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = {}
            for key in ("genome_fasta", "annotation_gtf", "whippet_gtf", "me_db", "transcriptome_fasta",
                        "decoys", "salmon_gtf", "splice_sites"):
                (root / key).write_text(key)
                reference[key] = str(root / key)
            reference["versions"] = {"hisat2": "2.2.1"}
            identity, settings = reference_identity(reference)
            before = reference_manifest(identity, versions=reference["versions"], settings=settings)
            # files the workflow builds later are recorded but do not change the ID
            (root / "tx2gene.tsv").write_text("built")
            built = dict(identity, tx2gene=str(root / "tx2gene.tsv"))
            after = reference_manifest(built, versions=reference["versions"], identity=identity,
                                       settings=settings)
            self.assertEqual(before["reference_id"], after["reference_id"])
            self.assertIn("tx2gene", after["checksums"])
            # the command line prints the same ID from a config file
            config = root / "config.json"
            config.write_text(json.dumps({"umbrella_reference": reference}))
            with unittest.mock.patch("sys.stdout", new_callable=io.StringIO) as stdout:
                shard_guard_main(["reference-id", "--configfile", str(config)])
            self.assertEqual(stdout.getvalue().strip(), before["reference_id"])
            # a prebuilt index or turning off the microexon insertion gives another reference
            off = reference_manifest(*reference_identity(reference, insert_microexons=False)[:1],
                                     versions=reference["versions"],
                                     settings=reference_identity(reference, False)[1])
            self.assertNotEqual(off["reference_id"], before["reference_id"])

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
