# SPDX-License-Identifier: MIT
"""The capability manifest (aipe.yaml) and the README agree with the package metadata, and the evidence files exist."""

from pathlib import Path

import pytest
import yaml

import aipe_thermal

ROOT = Path(__file__).resolve().parents[1]
REQUIRED_KEYS = ["schema_version", "id", "name", "type", "domains", "description", "status", "maturity", "repository", "license",
                 "capabilities", "inputs", "outputs", "dependencies", "compatibility", "tools", "connectors", "validation", "documentation"]


@pytest.fixture(scope="module")
def manifest() -> dict:
    return yaml.safe_load((ROOT / "aipe.yaml").read_text(encoding="utf-8"))


def test_identity_and_required_keys(manifest):
    assert all(key in manifest for key in REQUIRED_KEYS)
    assert manifest["schema_version"] == "0.1.0" and manifest["id"] == "aipe.thermal-design" and manifest["type"] == "tool"
    assert manifest["repository"] == "https://github.com/FulongLi/AIPE-Thermal-Design"
    assert manifest["license"] == "CC-BY-4.0 AND MIT"
    ext = manifest["extensions"]["aipe.thermal-design"]
    assert ext["package"] == f"aipe-thermal {aipe_thermal.__version__}" and ext["python_package"] == "aipe_thermal"
    assert ext["heatsink_model_version"] == aipe_thermal.heatsink.__version__
    assert manifest["compatibility"] == {"core": "0.1.x", "engineering_state": "0.1.x", "integration": "mapped"}


def test_capabilities_cover_the_required_areas(manifest):
    caps = set(manifest["capabilities"])
    assert {"thermal-network-steady-state", "thermal-network-transient-foster-cauer", "electrothermal-iteration",
            "cooling-concept-comparison", "heatsink-correlation-model", "heatsink-pareto-optimisation",
            "heatsink-stage-a-screening", "thermal-design-constraints"} <= caps
    for tool in manifest["tools"]:
        assert set(tool["capabilities"]) <= caps
        assert tool["licence_class"] == "open_source" and not tool["optional"]


def test_evidence_and_documentation_paths_exist(manifest):
    assert manifest["validation"]["status"] == "tested"
    for ref in manifest["validation"]["evidence"] + manifest["documentation"]:
        assert (ROOT / ref).is_file(), ref
    assert any("No CFD, FEM or measured" in limit for limit in manifest["validation"]["limitations"])


def test_readme_states_the_package_version_and_both_licences():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert f"Version: {aipe_thermal.__version__}" in readme
    assert "CC BY 4.0" in readme and "MIT" in readme and "LICENSES.md" in readme
    assert "Developer: Fulong Li" in readme and "Status: Prototype" in readme
