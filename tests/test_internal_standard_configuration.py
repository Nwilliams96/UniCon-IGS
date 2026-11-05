"""Check the website's ordered-list ISD configuration and downstream names."""
from itertools import combinations
from pathlib import Path
import re
from types import SimpleNamespace
import unittest

from snakemake.exceptions import WorkflowError

ROOT = Path(__file__).resolve().parents[1]
COMMON = (ROOT / "workflow/rules/common.smk").read_text()


class InternalStandardConfigurationTests(unittest.TestCase):
    def resolve(self, config):
        block = COMMON[COMMON.index('INTERNAL_STANDARD_IDS = config.get("intstds", [])'):COMMON.index("\nAMPLICON_CONCENTRATIONS_FILE =")]
        scope = dict(config={"projectName": "fixture", **config},
                     USE_INTERNAL_STANDARDS=config.get("use_internal_standards", False),
                     WorkflowError=WorkflowError, re=re, combinations=combinations)
        exec(compile(block, "common.smk:internal-standards", "exec"), scope)
        return scope

    def test_current_list_configs_preserve_ids_order_and_named_outputs(self):
        for ids in (["BP"], ["BP", "DR"], ["BP", "DR", "TT"],
                    ["Custom-A", "Custom_B", "SpikeIn3", "SpikeIn4"]):
            with self.subTest(ids=ids):
                result = self.resolve(dict(use_internal_standards=True, intstds=ids))
                self.assertEqual(result["INTERNAL_STANDARD_IDS"], ids)
                self.assertEqual(result["INTERNAL_STANDARD_SLOTS"],
                                 [f"isd_{i}" for i in range(1, len(ids) + 1)])
                self.assertEqual(result["INTERNAL_STANDARD_ID_PAIRS"], list(combinations(ids, 2)))
                self.assertEqual(result["ALL_INTERNAL_STANDARD_IDS_STEM"], "_".join(ids))
                self.assertEqual(result["ISD_CORRECTED_LONG_TABLE"],
                                 "results/05-internal-std-corrected/fixture."
                                 + "_".join(ids) + ".ISD_corrected_long_data.tsv")
                self.assertEqual(result["INTERNAL_STANDARD_METHOD_STEMS"],
                                 [f"{id}_recovery_ratio" for id in ids]
                                 + [f"mean_{'_'.join(ids)}_recovery_ratio",
                                    f"median_{'_'.join(ids)}_recovery_ratio"]
                                 + [f"mean_{a}_and_{b}_recovery_ratio" for a,b in combinations(ids,2)])

    def test_configured_list_is_used_directly_without_a_normalization_copy(self):
        ids = ["Custom-A", "Custom_B", "001"]
        result = self.resolve(dict(use_internal_standards=True, intstds=ids))
        self.assertIs(result["INTERNAL_STANDARD_IDS"], ids)

    def test_missing_or_empty_list_remains_valid_when_standards_are_disabled(self):
        for config in ({}, {"use_internal_standards": False},
                       {"use_internal_standards": False, "intstds": []}):
            with self.subTest(config=config):
                self.assertEqual(self.resolve(config)["INTERNAL_STANDARD_IDS"], [])

    def test_enabled_standards_still_require_nonempty_unique_safe_ids(self):
        for config, message in (
            ({"use_internal_standards": True}, "does not contain any IDs"),
            ({"use_internal_standards": True, "intstds": []}, "does not contain any IDs"),
            ({"use_internal_standards": True, "intstds": [""]}, "may contain only"),
            ({"use_internal_standards": True, "intstds": ["bad name"]}, "may contain only"),
            ({"use_internal_standards": True, "intstds": [" BP "]}, "may contain only"),
            ({"use_internal_standards": True, "intstds": ["BP", "BP"]}, "must be unique"),
        ):
            with self.subTest(config=config), self.assertRaisesRegex(WorkflowError, message):
                self.resolve(config)

    def test_sample_sheet_amount_columns_are_still_validated_before_processing(self):
        block = COMMON.split("# read sample sheet", 1)[1]
        block = block[block.index("if USE_INTERNAL_STANDARDS:"):block.index("def get_final_output():")]
        for enabled, columns, missing in (
            (True, ["sample", "BP_ng", "DR_ng"], None),
            (True, ["sample", "BP_ng"], "DR_ng"),
            (True, ["sample"], "BP_ng, DR_ng"),
            (False, ["sample"], None),
        ):
            with self.subTest(enabled=enabled, columns=columns):
                scope = dict(
                    USE_INTERNAL_STANDARDS=enabled,
                    INTERNAL_STANDARD_IDS=["BP", "DR"],
                    samples=SimpleNamespace(columns=columns),
                    WorkflowError=WorkflowError,
                )
                if missing:
                    with self.assertRaisesRegex(WorkflowError, re.escape(missing)):
                        exec(compile(block, "common.smk:sample-columns", "exec"), scope)
                else:
                    exec(compile(block, "common.smk:sample-columns", "exec"), scope)

    def test_project_name_uses_the_explicit_config_value_and_keeps_name_validation(self):
        block = COMMON[COMMON.index("PROJECT_NAME ="):COMMON.index("RESULTS_LONG_DATA =")]
        self.assertIn('PROJECT_NAME = config["projectName"]', block)
        config = {"projectName": "Configured-Project_001"}
        scope = dict(config=config, WorkflowError=WorkflowError, re=re)
        exec(compile(block, "common.smk:project-name", "exec"), scope)
        self.assertEqual(scope["PROJECT_NAME"], "Configured-Project_001")
        self.assertEqual(scope["RESULTS_EXPORT_DIR"], "Configured-Project_001-Results-Export")
        with self.assertRaises(KeyError):
            exec(compile(block, "common.smk:project-name", "exec"),
                 dict(config={}, WorkflowError=WorkflowError, re=re))
        for invalid_name in ("", "bad name", " bad", "../bad", "bad/project"):
            with self.subTest(invalid_name=invalid_name), self.assertRaisesRegex(
                WorkflowError, "projectName may contain only"
            ):
                exec(compile(block, "common.smk:project-name", "exec"),
                     dict(config={"projectName": invalid_name}, WorkflowError=WorkflowError, re=re))


if __name__ == "__main__":
    unittest.main()
