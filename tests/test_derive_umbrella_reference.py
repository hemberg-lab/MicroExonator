"""Tiny reference fixtures; no genome indexing or cluster tools required."""

import gzip
import tempfile
import unittest
from pathlib import Path

from src.derive_umbrella_reference import decoy_names, splice_sites, validate_transcripts
from src.shard_guard import reference_identity


class DeriveUmbrellaReferenceTest(unittest.TestCase):
    def test_genome_names_become_decoys_in_fasta_order(self):
        with tempfile.TemporaryDirectory() as directory:
            genome = Path(directory) / "genome.fa.gz"
            with gzip.open(genome, "wt") as output:
                output.write(">chr1 description\nACGT\n>chr2\nT\n")
            self.assertEqual(decoy_names(genome), "chr1\nchr2\n")

    def test_annotation_introns_become_hisat2_hints(self):
        with tempfile.TemporaryDirectory() as directory:
            gtf = Path(directory) / "genes.gtf"
            gtf.write_text(
                'chr1\ttest\texon\t11\t20\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
                'chr1\ttest\texon\t31\t40\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
            self.assertEqual(splice_sites(gtf), "chr1\t19\t30\t+\n")

    def test_transcript_fasta_must_match_gene_map(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "transcripts.fa"
            mapping = Path(directory) / "tx2gene.tsv"
            fasta.write_text(">t1 gene=g1\nACGT\n")
            mapping.write_text("TXNAME\tGENEID\nt1\tg1\n")
            validate_transcripts(fasta, mapping)
            mapping.write_text("TXNAME\tGENEID\nt1\tg1\nt2\tg2\n")
            with self.assertRaisesRegex(ValueError, "1 missing"):
                validate_transcripts(fasta, mapping)

    def test_reference_identity_needs_only_supplied_inputs(self):
        reference = dict(genome_fasta="genome.fa", annotation_gtf="genes.gtf",
                         whippet_gtf="whippet.gtf", me_db="me_db.txt")
        paths, _ = reference_identity(reference)
        self.assertEqual(set(paths), {"genome", "annotation", "whippet_annotation", "me_db"})


if __name__ == "__main__":
    unittest.main()
