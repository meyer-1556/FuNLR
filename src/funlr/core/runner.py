"""Argument-list command interface with checked, recorded execution."""
import json
import logging
from pathlib import Path
import shutil
import subprocess
import time

from .errors import ToolNotFound


class CommandRunner:
    def __init__(self, tools, commands_log, logger=None, dry_run=False, inventory=None):
        self.tools = dict(tools)
        self.commands_log = Path(commands_log)
        self.logger = logger or logging.getLogger("funlr.commands")
        self.dry_run = dry_run
        self.inventory = inventory or {}

    def resolve(self, key):
        candidate = self.tools.get(key, key)
        executable = shutil.which(str(candidate))
        if not executable:
            raise ToolNotFound(f"Required tool {key!r} not executable: {candidate}")
        return str(Path(executable).absolute())

    def has_tool(self, key):
        try:
            self.resolve(key)
            return True
        except ToolNotFound:
            return False

    def _record(self, record):
        self.commands_log.parent.mkdir(parents=True, exist_ok=True)
        with self.commands_log.open("a") as handle:
            handle.write(json.dumps(record) + "\n")

    def run(self, argv, *, cwd=None, stdout_path=None, stderr_path=None, check=True, env=None):
        argv = list(map(str, argv))
        started = time.time()
        record = {"kind": "external_tool", "argv": argv, "cwd": str(cwd) if cwd else None,
                  "stdout_path": str(stdout_path) if stdout_path else None,
                  "stderr_path": str(stderr_path) if stderr_path else None,
                  "dry_run": self.dry_run, "started": started, "returncode": None}
        if self.dry_run:
            return subprocess.CompletedProcess(argv, 0, "", "")
        self._record(record)
        self.logger.info("RUN: %s", argv)
        out = open(stdout_path, "w") if stdout_path else subprocess.PIPE
        err = open(stderr_path, "w") if stderr_path else subprocess.PIPE
        try:
            # Inherit the stage worker's process group for whole-stage cleanup.
            result = subprocess.run(argv, cwd=cwd, stdout=out, stderr=err, env=env,
                                    text=True, errors="replace", check=False)
            record["returncode"] = result.returncode
        except BaseException as exc:
            record["error"] = str(exc) or type(exc).__name__
            raise
        finally:
            if stdout_path:
                out.close()
            if stderr_path:
                err.close()
            record["duration_s"] = time.time() - started
            self._record(record)
        if check and result.returncode:
            if result.stderr:
                self.logger.error("Tool stderr: %s", result.stderr[-12000:].rstrip())
            if result.stdout:
                self.logger.error("Tool stdout: %s", result.stdout[-12000:].rstrip())
            raise subprocess.CalledProcessError(result.returncode, argv, output=result.stdout, stderr=result.stderr)
        return result

    def tool_version(self, key):
        if key in self.inventory:
            return self.inventory[key]["version_output"]
        from funlr.runner import TOOLS
        flags = TOOLS.get(key, ["--version"])
        result = subprocess.run([self.resolve(key), *flags], capture_output=True, text=True, errors="replace", timeout=15)
        return (result.stdout or result.stderr).strip() or "unknown"
