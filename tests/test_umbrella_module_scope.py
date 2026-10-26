import tempfile
import unittest
from pathlib import Path

from src import umbrella_modules as um
from src.umbrella_execution import COMPLEMENT, execution_tools, tool_for

HEADER = "\t".join(um.SELECTION_COLUMNS) + "\n"


class ScopeTests(unittest.TestCase):
    def test_defaults_keep_existing_modules_on_every_comparison(self):
        scope = um.parse_scope({})
        self.assertEqual({m for m, s in scope.items() if s == "selected"}, {"ciriquant", "circexplorer2", "dcc"})
        self.assertEqual(scope["vasttools"], "all")
        self.assertEqual(um.parse_scope({"umbrella_module_scope": {"ciriquant": "all"}})["ciriquant"], "all")

    def test_scope_errors(self):
        for value, message in (({"nope": "all"}, "unknown modules nope"), ({"dcc": "some"}, "all or selected"),
                               ("dcc", "must be a mapping")):
            with self.assertRaisesRegex(ValueError, message):
                um.parse_scope({"umbrella_module_scope": value})

    def test_internal_star_cannot_be_selected(self):
        with self.assertRaisesRegex(ValueError, "unknown umbrella_modules circrna_star"):
            um.parse_selection({"umbrella_modules": ["circrna_star"]})


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "selection.tsv"
        self.scope = um.parse_scope({})

    def tearDown(self):
        self.tmp.cleanup()

    def load(self, *rows):
        self.path.write_text(HEADER + "# agents may annotate\n" + "".join("\t".join(row) + "\n" for row in rows))
        return um.load_module_selection(self.path, {"c1", "c2"}, self.scope)

    def test_rows_are_grouped_by_tool(self):
        selection = self.load(("c1", "ciriquant", "test", "Ribo-Zero", "agent:x", "guillermo", "pilot"),
                              ("c2", "dcc", "polya_reference", "TruSeq mRNA", "guillermo", "guillermo", "ref"))
        self.assertEqual(set(selection), {"ciriquant", "dcc"})
        self.assertEqual(selection["ciriquant"]["c1"]["decided_by"], "agent:x")

    def test_missing_columns_and_no_file(self):
        self.path.write_text("comparison_id\ttool\n")
        with self.assertRaisesRegex(ValueError, "lacks columns"):
            um.load_module_selection(self.path, {"c1"}, self.scope)
        self.assertEqual(um.load_module_selection("", {"c1"}, self.scope), {})


class CatalogueIdentityTests(unittest.TestCase):
    def test_sources_hash_and_liftover_need(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "a.bed"
            source.write_text("chr1\t1\t9\tx\t0\t+\n")
            hg38 = {"circrna_catalogue_sources": [{"path": str(source)}]}
            settings = um.circrna_catalogue_settings(hg38)
            self.assertEqual(settings["liftover_chain"], "")
            self.assertEqual(len(settings["sources"][0]["sha256"]), 64)
            first = um.circrna_catalogue_id("ref", hg38)
            source.write_text("chr1\t1\t10\tx\t0\t+\n")
            self.assertNotEqual(um.circrna_catalogue_id("ref", hg38), first)   # content, not path
            hg19 = {"circrna_catalogue_sources": [{"path": str(source), "assembly": "hg19"}]}
            self.assertEqual(um.circrna_catalogue_settings(hg19)["liftover_chain"], "missing")
            # the catalogue enters CIRIquant's identity only
            self.assertNotEqual(um.module_reference_id("ref", "ciriquant", hg38),
                                um.module_reference_id("ref", "ciriquant", {}))
            self.assertEqual(um.module_reference_id("ref", "dcc", hg38), um.module_reference_id("ref", "dcc", {}))


class ExecutionRegistryTests(unittest.TestCase):
    def test_circrna_tools_are_explicit_only(self):
        for tool in ("ciriquant", "circexplorer2", "dcc"):
            self.assertIn(tool, COMPLEMENT)
            self.assertNotIn(tool, execution_tools("complement", {}))
        self.assertEqual(execution_tools("complement", {"tools": ["dcc"]}), ["dcc"])

    def test_rules_map_to_their_tools(self):
        self.assertEqual(tool_for("umbrella_ciriquant_prep", "x"), "ciriquant")
        self.assertEqual(tool_for("umbrella_circrna_star", "x"), "circexplorer2")
        self.assertEqual(tool_for("umbrella_circexplorer2_calls", "x"), "circexplorer2")
        self.assertEqual(tool_for("umbrella_dcc_compare", "x"), "dcc")


if __name__ == "__main__":
    unittest.main()
