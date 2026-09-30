#!/usr/bin/env python3

import copy
import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "adapter_contract", ROOT / "tests/check_adapter_contract.py"
)
CONTRACT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONTRACT)


class AdapterContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel = json.loads((ROOT / ".codeocean/app-panel.json").read_text())
        cls.main_text = (ROOT / "code/main.R").read_text()
        cls.source_text = (ROOT / "OMIX_MODULE_SOURCE.md").read_text()
        cls.datasets = json.loads((ROOT / ".codeocean/datasets.json").read_text())

    def audit(self, panel=None, main_text=None, source_text=None, datasets=None):
        return CONTRACT.audit(
            panel if panel is not None else self.panel,
            main_text if main_text is not None else self.main_text,
            source_text if source_text is not None else self.source_text,
            datasets if datasets is not None else self.datasets,
        )

    def test_repository_contract_passes(self):
        self.assertEqual(self.audit(), [])

    def test_missing_panel_control_is_detected(self):
        panel = copy.deepcopy(self.panel)
        panel["parameters"] = panel["parameters"][:-1]
        self.assertTrue(any("names/order" in item for item in self.audit(panel=panel)))

    def test_choice_drift_is_detected(self):
        panel = copy.deepcopy(self.panel)
        panel["parameters"][12]["extra_data"] = ["red/blue by ES", "green"]
        self.assertTrue(any("choices mismatch" in item for item in self.audit(panel=panel)))

    def test_numeric_bound_drift_is_detected(self):
        panel = copy.deepcopy(self.panel)
        by_name = {item["param_name"]: item for item in panel["parameters"]}
        by_name["max_le_genes_heatmap"]["minimum"] = 5
        by_name["pdf_width"]["maximum"] = 20
        findings = self.audit(panel=panel)
        self.assertTrue(any("max_le_genes_heatmap: App Panel minimum mismatch" in item for item in findings))
        self.assertTrue(any("pdf_width: App Panel maximum mismatch" in item for item in findings))

    def test_duplicate_dataset_mount_is_detected(self):
        datasets = copy.deepcopy(self.datasets)
        datasets["attached_datasets"].append(copy.deepcopy(datasets["attached_datasets"][0]))
        self.assertTrue(any("data mounts" in item for item in self.audit(datasets=datasets)))

    def test_undocumented_source_hash_is_detected(self):
        source = self.source_text.replace(CONTRACT.SCIENTIFIC_SHA256, "0" * 64)
        self.assertTrue(any("export hash" in item for item in self.audit(source_text=source)))


if __name__ == "__main__":
    unittest.main()
