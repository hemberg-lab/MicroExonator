import unittest

from src.whippet_flags import quant_flags, read_paths


class WhippetFlagsTests(unittest.TestCase):
    def test_circular_output_is_enabled_by_default(self):
        self.assertEqual(quant_flags(), "--circ")

    def test_preserves_other_flags_and_does_not_duplicate_circ(self):
        self.assertEqual(quant_flags("--stranded"), "--stranded --circ")
        self.assertEqual(quant_flags("--stranded --circ"), "--stranded --circ")
        self.assertEqual(quant_flags(None), "--circ")

    def test_backfill_reads_have_a_separate_permanent_processing_path(self):
        from types import SimpleNamespace
        run = SimpleNamespace(work_dir="umbrella/work/ref/P/batch/r", metadata={}, layout="PE")
        self.assertEqual(read_paths(run), [run.work_dir + "/R1.fastq.gz", run.work_dir + "/R2.fastq.gz"])
        run.metadata["whippet_read_source"] = "backfill"
        self.assertEqual(read_paths(run), ["umbrella/whippet_backfill/ref/P/batch/r/R1.fastq.gz",
                                           "umbrella/whippet_backfill/ref/P/batch/r/R2.fastq.gz"])
        run.metadata["whippet_read_source"] = "unknown"
        with self.assertRaises(ValueError):
            read_paths(run)
