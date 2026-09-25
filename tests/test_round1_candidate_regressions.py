import pathlib
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]


class CandidateGenerationRegressionTests(unittest.TestCase):
    def test_annotation_tag_context_keeps_every_annotated_path(self):
        # The Python 2 original kept the last transcript's flank per boundary
        # and the Python 3 port the longest; every distinct path is now kept
        # (behaviour is tested in test_synthetic_skip_tags).
        source = (REPOSITORY / "src" / "Get_annotated_microexons.py").read_text()
        self.assertIn("SJ_start_seqs.setdefault((chrom, eend), set()).add(f_seq[-flank:])", source)
        self.assertIn("SJ_end_seqs.setdefault((chrom, estart), set()).add(r_seq[:flank])", source)

    def test_round1_rejects_one_sided_insertion_evidence(self):
        source = (REPOSITORY / "src" / "ME_filter1.py").read_text()
        self.assertIn("start + sum(match_lengths) >= up + 8", source)

    def test_annotated_candidates_take_precedence_over_denovo_candidates(self):
        source = (REPOSITORY / "rules" / "Round2.smk").read_text()
        self.assertIn("awk 'NR==FNR {{seen[$1]=1; print; next}} !seen[$1]'", source)


if __name__ == "__main__":
    unittest.main()
