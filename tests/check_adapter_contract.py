#!/usr/bin/env python3
"""Check this Code Ocean adapter against the reviewed canonical interface."""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANONICAL_COMMIT = "eee4433cd74d0ce9bd1daba1571ca6312cba7ea2"
SCIENTIFIC_SHA256 = "5fdae88b3968443a492b654d1d0810292fea4cb8ed70242acb4de17f1ffe4e5a"
SCHEMA_SHA256 = "f322dcfbb48616aba755cd51c94772e83a24f61998fc75b69562732fc36342f4"

# Hidden inputs are fixed Code Ocean data ports. Every remaining public or
# advanced control appears below in canonical schema order.
HIDDEN_INPUTS = ["msigdb_database", "gsea_filter_results"]
CONTRACT = [
    {"name": "deg_table", "r_type": "character", "panel_type": "file", "default": None},
    {"name": "sample_metadata", "r_type": "character", "panel_type": "file", "default": None},
    {"name": "contrast_filter", "r_type": "character", "panel_type": "list", "default": "none", "choices": ["none", "keep", "remove"]},
    {"name": "contrasts", "r_type": "character", "panel_type": "text", "default": None},
    {"name": "top_n_pathways", "r_type": "integer", "panel_type": "text", "default": 20, "minimum": 0},
    {"name": "top_n_by_sign", "r_type": "logical", "panel_type": "list", "default": False, "ui_choices": ["FALSE", "TRUE"]},
    {"name": "pathway_bubble_plots", "r_type": "logical", "panel_type": "list", "default": True, "ui_choices": ["TRUE", "FALSE"]},
    {"name": "pathway_bubble_top_n", "r_type": "integer", "panel_type": "text", "default": 20, "minimum": 0},
    {"name": "pathway_bubble_significance_statistic", "r_type": "character", "panel_type": "list", "default": "padj", "choices": ["padj", "pval"]},
    {"name": "collection_color_scale", "r_type": "character", "panel_type": "list", "default": "independent", "choices": ["independent", "shared"]},
    {"name": "max_plots_in_pdf", "r_type": "integer", "panel_type": "text", "default": 0, "minimum": 0},
    {"name": "plots_to_include", "r_type": "character", "panel_type": "list", "default": "ES+RNK+LE", "choices": ["ES", "ES+RNK", "ES+LE", "ES+RNK+LE", "LE"]},
    {"name": "running_score_line_color", "r_type": "character", "panel_type": "list", "default": "ES sign", "choices": ["ES sign", "green"]},
    {"name": "add_max_deviation_line", "r_type": "character", "panel_type": "list", "default": "both", "choices": ["coordinate", "horizontal", "both", "none"]},
    {"name": "rank_area_color", "r_type": "character", "panel_type": "list", "default": "red/blue by Gene score", "choices": ["grey", "red/blue by Gene score"]},
    {"name": "show_es_rank_bar", "r_type": "logical", "panel_type": "list", "default": False, "ui_choices": ["FALSE", "TRUE"]},
    {"name": "show_rnk_peak_line", "r_type": "logical", "panel_type": "list", "default": True, "ui_choices": ["TRUE", "FALSE"]},
    {"name": "show_rnk_le_highlight", "r_type": "logical", "panel_type": "list", "default": True, "ui_choices": ["TRUE", "FALSE"]},
    {"name": "show_es_le_highlight", "r_type": "logical", "panel_type": "list", "default": True, "ui_choices": ["TRUE", "FALSE"]},
    {"name": "heatmap_transform", "r_type": "character", "panel_type": "list", "default": "z-score", "choices": ["z-score", "center by row mean", "center by row median", "none"]},
    {"name": "max_le_genes_heatmap", "r_type": "integer", "panel_type": "text", "default": 50, "minimum": 1},
    {"name": "heatmap_gene_order", "r_type": "character", "panel_type": "list", "default": "rank", "choices": ["rank", "cluster", "input"]},
    {"name": "heatmap_sample_order", "r_type": "character", "panel_type": "list", "default": "group", "choices": ["group", "cluster", "input"]},
    {"name": "heatmap_gene_clustering_distance", "r_type": "character", "panel_type": "list", "default": "euclidean", "choices": ["euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"]},
    {"name": "heatmap_gene_clustering_method", "r_type": "character", "panel_type": "list", "default": "complete", "choices": ["ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"]},
    {"name": "heatmap_sample_clustering_distance", "r_type": "character", "panel_type": "list", "default": "euclidean", "choices": ["euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"]},
    {"name": "heatmap_sample_clustering_method", "r_type": "character", "panel_type": "list", "default": "complete", "choices": ["ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"]},
    {"name": "show_le_heatmap_gene_names", "r_type": "logical", "panel_type": "list", "default": True, "ui_choices": ["TRUE", "FALSE"]},
    {"name": "show_le_heatmap_sample_names", "r_type": "logical", "panel_type": "list", "default": False, "ui_choices": ["FALSE", "TRUE"]},
    {"name": "show_le_heatmap_rank_labels", "r_type": "logical", "panel_type": "list", "default": True, "ui_choices": ["TRUE", "FALSE"]},
    {"name": "heatmap_gene_names_column", "r_type": "character", "panel_type": "text", "default": "GeneName"},
    {"name": "heatmap_sample_names_column", "r_type": "character", "panel_type": "text", "default": "Sample"},
    {"name": "heatmap_group_column", "r_type": "character", "panel_type": "text", "default": "Group"},
    {"name": "pdf_width", "r_type": "numeric", "panel_type": "text", "default": 8.5},
    {"name": "pdf_height", "r_type": "numeric", "panel_type": "text", "default": 6.5},
]

CLI_ORDER = [
    "msigdb_database", "gsea_filter_results", "deg_table", "sample_metadata",
    "output_dir", *[item["name"] for item in CONTRACT[2:]],
]


def comparable(value):
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if text.lower() in {"true", "false"}:
        return text.lower() == "true"
    if re.fullmatch(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", text):
        return float(text)
    return text


def balanced_calls(text: str, name: str) -> list[str]:
    calls, needle, start = [], name + "(", 0
    while True:
        index = text.find(needle, start)
        if index < 0:
            return calls
        depth, quote, escaped = 0, None, False
        for end in range(index + len(name), len(text)):
            char = text[end]
            if quote:
                if escaped:
                    escaped = False
                elif char == "\\":
                    escaped = True
                elif char == quote:
                    quote = None
                continue
            if char in {"'", '"'}:
                quote = char
            elif char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
                if depth == 0:
                    calls.append(text[index : end + 1])
                    start = end + 1
                    break


def parse_literal(raw: str | None):
    if raw is None or raw.strip() in {"NULL", "NA", "NA_character_"}:
        return None
    raw = raw.strip()
    if raw in {"TRUE", "FALSE"}:
        return raw == "TRUE"
    if re.fullmatch(r"[-+]?\d+[Ll]", raw):
        return float(raw[:-1])
    if re.fullmatch(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", raw):
        return float(raw)
    if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in {"'", '"'}:
        return raw[1:-1]
    return raw


def named_arg(call: str, name: str) -> str | None:
    match = re.search(
        rf"(?:^|,)\s*{re.escape(name)}\s*=\s*((?:'(?:\\.|[^'])*'|\"(?:\\.|[^\"])*\"|[^,)]+))",
        call,
        re.DOTALL,
    )
    return match.group(1).strip() if match else None


def parse_options(text: str) -> list[dict]:
    options = []
    for call in balanced_calls(text, "make_option"):
        match = re.search(r"['\"]--([A-Za-z0-9_]+)['\"]", call)
        if match:
            options.append({
                "name": match.group(1),
                "type": parse_literal(named_arg(call, "type")),
                "default": parse_literal(named_arg(call, "default")),
            })
    return options


def audit(panel: dict, main_text: str, source_text: str, datasets: dict) -> list[str]:
    findings = []
    expected_names = [item["name"] for item in CONTRACT]
    parameters = panel.get("parameters", [])
    panel_names = [item.get("param_name") for item in parameters]

    if panel.get("named_parameters") is not True:
        findings.append("named_parameters must be true")
    if panel_names != expected_names:
        findings.append("App Panel names/order differ from the canonical visible contract")
    if len(panel_names) != len(set(panel_names)):
        findings.append("duplicate App Panel param_name values")

    panel_by_name = {item.get("param_name"): item for item in parameters}
    for expected in CONTRACT:
        observed = panel_by_name.get(expected["name"])
        if observed is None:
            continue
        if observed.get("type") != expected["panel_type"]:
            findings.append(f"{expected['name']}: App Panel type mismatch")
        if observed.get("value_type") != ("number" if expected["panel_type"] == "text" and expected["r_type"] in {"integer", "numeric"} else "string"):
            findings.append(f"{expected['name']}: App Panel value_type mismatch")
        if comparable(observed.get("default_value")) != comparable(expected["default"]):
            findings.append(f"{expected['name']}: App Panel default mismatch")
        expected_choices = expected.get("choices", expected.get("ui_choices", []))
        if observed.get("extra_data", []) != expected_choices:
            findings.append(f"{expected['name']}: App Panel choices mismatch")
        for bound in ("minimum", "maximum"):
            if comparable(observed.get(bound)) != comparable(expected.get(bound)):
                findings.append(f"{expected['name']}: App Panel {bound} mismatch")

    options = parse_options(main_text)
    option_names = [item["name"] for item in options]
    if option_names != CLI_ORDER:
        findings.append("adapter CLI names/order differ from the canonical contract")
    option_by_name = {item["name"]: item for item in options}
    for expected in CONTRACT:
        observed = option_by_name.get(expected["name"])
        if observed is None:
            continue
        if observed["type"] != expected["r_type"]:
            findings.append(f"{expected['name']}: adapter CLI type mismatch")
        if comparable(observed["default"]) != comparable(expected["default"]):
            findings.append(f"{expected['name']}: adapter CLI default mismatch")

    mounts = [item.get("mount") for item in datasets.get("attached_datasets", [])]
    if mounts != ["deg-training", "gsea_filter_results", "msigdb"]:
        findings.append("attached data mounts are missing, duplicated, or out of reviewed order")
    if "GSEA_Visualization_Local_v1.R" in {path.name for path in (ROOT / "code").glob("*.R")}:
        findings.append("retired second scientific implementation is still active in code/")
    if 'source(core_file)' not in main_text or 'functions", "gsea_enrichment_plot.R"' not in main_text:
        findings.append("adapter does not source the managed canonical implementation")

    digest = hashlib.sha256((ROOT / "code/functions/gsea_enrichment_plot.R").read_bytes()).hexdigest()
    if digest != SCIENTIFIC_SHA256 or digest not in source_text:
        findings.append("managed scientific export hash is incorrect or undocumented")
    if CANONICAL_COMMIT not in source_text or SCHEMA_SHA256 not in source_text:
        findings.append("canonical commit or schema hash is missing from source record")
    for hidden in HIDDEN_INPUTS:
        if f'"canonical":"{hidden}"' not in source_text:
            findings.append(f"hidden input binding is undocumented: {hidden}")
    return findings


def audit_repository() -> list[str]:
    return audit(
        json.loads((ROOT / ".codeocean/app-panel.json").read_text()),
        (ROOT / "code/main.R").read_text(),
        (ROOT / "OMIX_MODULE_SOURCE.md").read_text(),
        json.loads((ROOT / ".codeocean/datasets.json").read_text()),
    )


if __name__ == "__main__":
    failures = audit_repository()
    if failures:
        for failure in failures:
            print(f"ERROR: {failure}", file=sys.stderr)
        raise SystemExit(1)
    print("GSEA Visualization adapter contract checks passed")
