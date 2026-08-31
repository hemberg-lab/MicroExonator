import os
import tempfile
import unittest


class FilterMethodTests(unittest.TestCase):
    def test_default_filter_method_is_robustness(self):
        try:
            from src.filtering_config import resolve_filter_method
        except ImportError as exc:
            self.fail("filtering configuration is not implemented: {}".format(exc))

        self.assertEqual(resolve_filter_method({}), "robustness")

    def test_legacy_mixture_requires_explicit_selection(self):
        from src.filtering_config import resolve_filter_method

        self.assertEqual(
            resolve_filter_method({"filter_method": "legacy_mixture"}),
            "legacy_mixture",
        )

    def test_unknown_filter_method_is_rejected(self):
        from src.filtering_config import resolve_filter_method

        with self.assertRaisesRegex(ValueError, "filter_method"):
            resolve_filter_method({"filter_method": "automatic"})

    def test_skip_mixture_legacy_setting_maps_to_robustness(self):
        from src.filtering_config import resolve_filter_method

        with self.assertWarns(DeprecationWarning):
            method = resolve_filter_method({"skip_mixture_model_filter": "T"})

        self.assertEqual(method, "robustness")

    def test_original_filter_mode_maps_to_legacy_mixture(self):
        from src.filtering_config import resolve_filter_method

        with self.assertWarns(DeprecationWarning):
            method = resolve_filter_method({"filter_mode": "original"})

        self.assertEqual(method, "legacy_mixture")

    def test_selected_output_matches_filter_method(self):
        try:
            from src.filtering_config import selected_microexon_output
        except ImportError as exc:
            self.fail("filter output selection is not implemented: {}".format(exc))

        self.assertEqual(
            selected_microexon_output("robustness"),
            "Report/out.robustly_detected.txt",
        )
        self.assertEqual(
            selected_microexon_output("legacy_mixture"),
            "Report/out.high_quality.txt",
        )


class FilterGroupTests(unittest.TestCase):
    def test_bulk_run_requires_bulk_samples_manifest(self):
        try:
            from src.filtering_config import load_filter_groups
        except ImportError as exc:
            self.fail("filter group loading is not implemented: {}".format(exc))

        with self.assertRaisesRegex(ValueError, "bulk_samples.tsv"):
            load_filter_groups({"Single_Cell": "F"}, {"sample_a"})

    def test_unclustered_single_cell_run_uses_one_all_cells_group(self):
        from src.filtering_config import load_filter_groups

        groups = load_filter_groups(
            {"Single_Cell": "T"},
            {"cell_2", "cell_1"},
        )

        self.assertEqual(groups["single_cell"], {"all_cells": ["cell_1", "cell_2"]})

    def test_cluster_metadata_defines_single_cell_filter_groups(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            metadata = os.path.join(temp_dir, "clusters.tsv")
            with open(metadata, "w") as handle:
                handle.write("cell\tcluster\n")
                handle.write("cell_2\tNeuron\n")
                handle.write("cell_1\tNeuron\n")
                handle.write("cell_3\tGlia\n")

            groups = load_filter_groups(
                {
                    "Single_Cell": "T",
                    "cluster_metadata": metadata,
                    "file_basename": "cell",
                    "cluster_name": "cluster",
                },
                {"cell_1", "cell_2", "cell_3"},
            )

        self.assertEqual(
            groups["single_cell"],
            {"Glia": ["cell_3"], "Neuron": ["cell_1", "cell_2"]},
        )

    def test_cluster_group_names_are_normalized_for_output_paths(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            metadata = os.path.join(temp_dir, "clusters.tsv")
            with open(metadata, "w") as handle:
                handle.write("cell\tcluster\n")
                handle.write("cell_a\tExcitatory neuron\n")

            groups = load_filter_groups(
                {
                    "Single_Cell": "T",
                    "cluster_metadata": metadata,
                    "file_basename": "cell",
                    "cluster_name": "cluster",
                },
                {"cell_a"},
            )

        self.assertEqual(groups["single_cell"], {"Excitatory_neuron": ["cell_a"]})

    def test_cluster_group_names_must_not_collide_after_normalization(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            metadata = os.path.join(temp_dir, "clusters.tsv")
            with open(metadata, "w") as handle:
                handle.write("cell\tcluster\n")
                handle.write("cell_a\tA B\n")
                handle.write("cell_b\tA_B\n")

            with self.assertRaisesRegex(ValueError, "collision"):
                load_filter_groups(
                    {
                        "Single_Cell": "T",
                        "cluster_metadata": metadata,
                        "file_basename": "cell",
                        "cluster_name": "cluster",
                    },
                    {"cell_a", "cell_b"},
                )

    def test_bulk_manifest_defines_single_and_paired_end_groups(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            manifest = os.path.join(temp_dir, "bulk_samples.tsv")
            with open(manifest, "w") as handle:
                handle.write("sample\tcondition\n")
                handle.write("single_a\tcontrol\n")
                handle.write("pair_1\ttreated\n")

            groups = load_filter_groups(
                {"Single_Cell": "F", "bulk_samples": manifest},
                {"single_a", "pair_1", "pair_2"},
                paired_dict={"pair_1": "pair_2"},
            )

        self.assertEqual(groups["bulk_se"], {"control": ["single_a"]})
        self.assertEqual(groups["bulk_pe"], {"treated": ["pair_1"]})

    def test_bulk_manifest_requires_sample_and_condition_columns(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            manifest = os.path.join(temp_dir, "bulk_samples.tsv")
            with open(manifest, "w") as handle:
                handle.write("sample\tgroup\n")
                handle.write("sample_a\tcontrol\n")

            try:
                load_filter_groups(
                    {"Single_Cell": "F", "bulk_samples": manifest},
                    {"sample_a"},
                )
            except ValueError as exc:
                self.assertRegex(str(exc), "sample.*condition")
            except Exception as exc:
                self.fail("manifest validation raised the wrong error: {}".format(exc))
            else:
                self.fail("manifest without a condition column was accepted")

    def test_bulk_manifest_rejects_unknown_samples(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            manifest = os.path.join(temp_dir, "bulk_samples.tsv")
            with open(manifest, "w") as handle:
                handle.write("sample\tcondition\n")
                handle.write("not_an_input\tcontrol\n")

            with self.assertRaisesRegex(ValueError, "not_an_input"):
                load_filter_groups(
                    {"Single_Cell": "F", "bulk_samples": manifest},
                    {"sample_a"},
                )

    def test_bulk_manifest_rejects_duplicate_samples(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            manifest = os.path.join(temp_dir, "bulk_samples.tsv")
            with open(manifest, "w") as handle:
                handle.write("sample\tcondition\n")
                handle.write("sample_a\tcontrol\n")
                handle.write("sample_a\ttreated\n")

            with self.assertRaisesRegex(ValueError, "duplicate.*sample_a"):
                load_filter_groups(
                    {"Single_Cell": "F", "bulk_samples": manifest},
                    {"sample_a"},
                )

    def test_bulk_manifest_rejects_empty_sample_or_condition(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            manifest = os.path.join(temp_dir, "bulk_samples.tsv")
            with open(manifest, "w") as handle:
                handle.write("sample\tcondition\n")
                handle.write("sample_a\t\n")

            with self.assertRaisesRegex(ValueError, "empty"):
                load_filter_groups(
                    {"Single_Cell": "F", "bulk_samples": manifest},
                    {"sample_a"},
                )

    def test_bulk_manifest_rejects_condition_names_unsafe_for_paths(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            manifest = os.path.join(temp_dir, "bulk_samples.tsv")
            with open(manifest, "w") as handle:
                handle.write("sample\tcondition\n")
                handle.write("sample_a\tcontrol/treated\n")

            with self.assertRaisesRegex(ValueError, "path-safe"):
                load_filter_groups(
                    {"Single_Cell": "F", "bulk_samples": manifest},
                    {"sample_a"},
                )

    def test_bulk_manifest_rejects_ungrouped_input_samples(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            manifest = os.path.join(temp_dir, "bulk_samples.tsv")
            with open(manifest, "w") as handle:
                handle.write("sample\tcondition\n")
                handle.write("sample_a\tcontrol\n")

            with self.assertRaisesRegex(ValueError, "sample_b"):
                load_filter_groups(
                    {"Single_Cell": "F", "bulk_samples": manifest},
                    {"sample_a", "sample_b"},
                )

    def test_mixed_run_without_cluster_metadata_is_rejected(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            manifest = os.path.join(temp_dir, "bulk_samples.tsv")
            with open(manifest, "w") as handle:
                handle.write("sample\tcondition\n")
                handle.write("bulk_a\tcontrol\n")

            with self.assertRaisesRegex(ValueError, "cluster_metadata"):
                load_filter_groups(
                    {"Single_Cell": "T", "bulk_samples": manifest},
                    {"bulk_a", "cell_a"},
                )

    def test_cluster_metadata_requires_configured_columns(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            metadata = os.path.join(temp_dir, "clusters.tsv")
            with open(metadata, "w") as handle:
                handle.write("cell\tcell_type\n")
                handle.write("cell_a\tNeuron\n")

            try:
                load_filter_groups(
                    {
                        "Single_Cell": "T",
                        "cluster_metadata": metadata,
                        "file_basename": "cell",
                        "cluster_name": "cluster",
                    },
                    {"cell_a"},
                )
            except ValueError as exc:
                self.assertRegex(str(exc), "cell.*cluster")
            except Exception as exc:
                self.fail("cluster validation raised the wrong error: {}".format(exc))
            else:
                self.fail("cluster metadata without the configured columns was accepted")

    def test_cluster_metadata_rejects_unknown_cells(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            metadata = os.path.join(temp_dir, "clusters.tsv")
            with open(metadata, "w") as handle:
                handle.write("cell\tcluster\n")
                handle.write("not_an_input\tNeuron\n")

            with self.assertRaisesRegex(ValueError, "not_an_input"):
                load_filter_groups(
                    {
                        "Single_Cell": "T",
                        "cluster_metadata": metadata,
                        "file_basename": "cell",
                        "cluster_name": "cluster",
                    },
                    {"cell_a"},
                )

    def test_cluster_metadata_rejects_duplicate_cells(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            metadata = os.path.join(temp_dir, "clusters.tsv")
            with open(metadata, "w") as handle:
                handle.write("cell\tcluster\n")
                handle.write("cell_a\tNeuron\n")
                handle.write("cell_a\tGlia\n")

            with self.assertRaisesRegex(ValueError, "duplicate.*cell_a"):
                load_filter_groups(
                    {
                        "Single_Cell": "T",
                        "cluster_metadata": metadata,
                        "file_basename": "cell",
                        "cluster_name": "cluster",
                    },
                    {"cell_a"},
                )

    def test_cluster_metadata_rejects_empty_cluster_names(self):
        from src.filtering_config import load_filter_groups

        with tempfile.TemporaryDirectory() as temp_dir:
            metadata = os.path.join(temp_dir, "clusters.tsv")
            with open(metadata, "w") as handle:
                handle.write("cell\tcluster\n")
                handle.write("cell_a\t\n")

            with self.assertRaisesRegex(ValueError, "empty"):
                load_filter_groups(
                    {
                        "Single_Cell": "T",
                        "cluster_metadata": metadata,
                        "file_basename": "cell",
                        "cluster_name": "cluster",
                    },
                    {"cell_a"},
                )

    def test_paired_group_spanning_reads_include_both_mates(self):
        try:
            from src.filtering_config import paired_read_samples
        except ImportError as exc:
            self.fail("paired read expansion is not implemented: {}".format(exc))

        self.assertEqual(
            paired_read_samples(["pair_a_1", "pair_b_1"], {
                "pair_a_1": "pair_a_2",
                "pair_b_1": "pair_b_2",
            }),
            ["pair_a_1", "pair_a_2", "pair_b_1", "pair_b_2"],
        )


if __name__ == "__main__":
    unittest.main()
