import pathlib
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]


class CandidateGenerationRegressionTests(unittest.TestCase):
    def test_annotation_tag_context_prefers_longest_sequence(self):
        source = (REPOSITORY / "src" / "Get_annotated_microexons.py").read_text()
        self.assertIn(
            "if len(f_seq[-100:]) > len(SJ_start_seqs[(chrom, eend )]):", source
        )
        self.assertIn(
            "if len(r_seq[:100]) > len(SJ_end_seqs[(chrom, estart )]):", source
        )

    def test_round1_rejects_one_sided_insertion_evidence(self):
        source = (REPOSITORY / "src" / "ME_filter1.py").read_text()
        self.assertIn("start + sum(match_lengths) >= up + 8", source)

    def test_annotated_candidates_take_precedence_over_denovo_candidates(self):
        source = (REPOSITORY / "rules" / "Round2.smk").read_text()
        self.assertIn("awk 'NR==FNR {{seen[$1]=1; print; next}} !seen[$1]'", source)


if __name__ == "__main__":
    unittest.main()
