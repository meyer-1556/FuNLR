"""Read-only installation diagnostics; no scientific inputs are required."""
from __future__ import annotations

import json
import platform
from pathlib import Path
import shutil
import subprocess
import sys

from . import __version__
from .runner import TOOLS as ALL_TOOLS, sha256, tool_inventory

TOOLS = {name: flags for name, flags in ALL_TOOLS.items() if name != "Rscript"}


_DEPENDENCIES = {"pandas": "pandas", "numpy": "numpy", "PyYAML": "yaml", "matplotlib": "matplotlib"}
_EXPECTED_VERSIONS = {"pandas": "2.2.3", "numpy": "1.26.4", "PyYAML": "6.0.3", "matplotlib": "3.10.8"}


def _python_dependencies() -> dict:
    """Check imports in this interpreter without writing bytecode caches."""
    script = '''import importlib, importlib.metadata, json
result = {}
for package, module in %r.items():
    try:
        if package == "matplotlib":
            # Importing matplotlib initializes caches. Rendering is checked by the demo.
            result[package] = {"ok": True, "version": importlib.metadata.version(package), "probe": "installed distribution metadata"}
        else:
            imported = importlib.import_module(module)
            result[package] = {"ok": True, "version": str(getattr(imported, "__version__", "unknown"))}
    except Exception as exc:
        result[package] = {"ok": False, "error": str(exc)}
print(json.dumps(result))
''' % _DEPENDENCIES
    try:
        result = subprocess.run([sys.executable, "-B", "-c", script], capture_output=True,
                                text=True, errors="replace", timeout=30, check=False)
        if result.returncode:
            raise ValueError(result.stderr.strip() or f"Python dependency check exited {result.returncode}")
        dependencies = json.loads(result.stdout)
    except (OSError, subprocess.TimeoutExpired, ValueError) as exc:
        dependencies = {name: {"ok": False, "error": str(exc)} for name in _DEPENDENCIES}
    for name, expected in _EXPECTED_VERSIONS.items():
        entry = dependencies[name]
        entry["expected_version"] = expected
        if entry["ok"] and entry.get("version") != expected:
            entry.update(ok=False, error=f"Expected {name}=={expected}; found {entry.get('version', 'unknown')}")
    return dependencies


def _probe_tool(name: str, flags: list[str], configured: dict) -> dict:
    """Give individual diagnostics if the shared all-tools check cannot finish."""
    executable = shutil.which(str(configured.get(name, name)))
    if not executable:
        return {"ok": False, "error": f"Executable not found: {configured.get(name, name)}"}
    entry = {"path": str(Path(executable).absolute()), "ok": False}
    try:
        result = subprocess.run([executable, *flags], stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, errors="replace", timeout=15, check=False)
        entry["version_output"] = result.stdout[:12000].strip()
        # Match the existing runner's accepted Exonerate 2.4.0 version exit.
        accepted_exit = (name == "exonerate" and result.returncode == 1
                         and "exonerate" in result.stdout.lower() and "2.4.0" in result.stdout)
        if result.returncode and not accepted_exit:
            entry["error"] = f"Version check exited {result.returncode}: {result.stdout[:400].strip()}"
        else:
            entry.update(ok=True, sha256=sha256(executable))
    except (OSError, subprocess.TimeoutExpired) as exc:
        entry["error"] = str(exc)
    return entry


def check_installation(*, include_demo: bool = False, configured_tools: dict | None = None, plot_backend: str = "matplotlib") -> dict:
    """Return JSON-serializable diagnostics without creating or changing files.

    ``include_demo`` additionally requires HMMER's ``hmmbuild``, used to create
    the small synthetic HMM databases for the bundled installation test.
    Configured executable paths use the same mapping as the pipeline runner.
    """
    if plot_backend not in ("matplotlib", "r"):
        raise ValueError("plot_backend must be matplotlib or r")
    required_tools = ALL_TOOLS if plot_backend == "r" else TOOLS
    configured = dict(configured_tools or {})
    deps = _python_dependencies()
    try:
        inventory = tool_inventory(configured, {"REPORT_PLOT_BACKEND": plot_backend})
        tools = {name: {"ok": True, **entry} for name, entry in inventory.items()}
    except (OSError, ValueError) as exc:
        tools = {name: _probe_tool(name, flags, configured) for name, flags in required_tools.items()}
        # Preserve a shared-runner failure even if a retry happened to succeed.
        inventory_error = str(exc)
    else:
        inventory_error = None
    if include_demo:
        tools["hmmbuild"] = _probe_tool("hmmbuild", ["-h"], configured)
    python_ok = (3, 11) <= sys.version_info < (3, 13)
    errors = []
    if not python_ok:
        errors.append("Python >=3.11,<3.13 is required")
    errors.extend(f"{name}: {item['error']}" for name, item in deps.items() if not item["ok"])
    errors.extend(f"{name}: {item.get('error', 'check failed')}" for name, item in tools.items() if not item["ok"])
    if inventory_error and not any(not item["ok"] for item in tools.values()):
        errors.append(inventory_error)
    return {
        "ok": not errors,
        "purpose": "Installation diagnostics; this does not validate biological results.",
        "funlr_version": __version__,
        "platform": {"system": platform.system(), "release": platform.release(), "machine": platform.machine()},
        "python": {"executable": sys.executable, "version": platform.python_version(),
                   "requires": ">=3.11,<3.13", "ok": python_ok},
        "python_dependencies": deps,
        "tools": tools,
        "includes_demo_dependencies": include_demo,
        "plot_backend": plot_backend,
        "errors": errors,
    }


def check_self_test(*, configured_tools: dict | None = None,
                    plot_backend: str = "matplotlib", models_dir=None) -> dict:
    """Check the installed runtime and packaged fixture without running analyses.

    Production HMMs are checked only at an explicitly supplied directory. Missing
    or damaged package resources fail on source and wheel installations alike.
    """
    from .public_demo import verified_assets
    from .databases import verify_models

    report = check_installation(include_demo=True, configured_tools=configured_tools,
                                plot_backend=plot_backend)
    report["self_test"] = True
    report["purpose"] = ("Runtime and packaged-fixture checks; no searches or rendering are run. "
                         "This is not a biological benchmark or complete production-input validation.")
    try:
        assets, manifest = verified_assets()
        report["public_fixture"] = {"status": "VERIFIED", "files_checked": len(assets),
                                    "sha256": manifest["sha256"]}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        report["public_fixture"] = {"status": "FAILED", "error": str(exc)}
        report["errors"].append(f"Packaged public fixture: {exc}")
    if models_dir is None:
        report["production_models"] = {
            "status": "NOT_REQUESTED",
            "reason": "Pass --models-dir to verify the five fixed reference HMM libraries. They are separate assets.",
        }
    else:
        try:
            report["production_models"] = verify_models(models_dir)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            report["production_models"] = {"status": "FAILED", "error": str(exc)}
            report["errors"].append(f"Production model snapshot: {exc}")
    report["ok"] = not report["errors"]
    return report
