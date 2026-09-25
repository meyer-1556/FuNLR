"""Run context: paths and shared services passed to every stage."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from funlr.config import FunlrConfig
    from funlr.core.runner import CommandRunner


@dataclass
class RunPaths:
    """Canonical output layout for one FuNLR run (legacy-compatible names)."""

    output_dir: Path

    @property
    def logs_dir(self) -> Path:
        p = self.output_dir / "logs"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def results_dir(self) -> Path:
        p = self.output_dir / "results"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def state_dir(self) -> Path:
        p = self.output_dir / "state"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def final_dir(self) -> Path:
        p = self.output_dir / "final_results"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def pressed_dir(self) -> Path:
        """Private working area for HMM copies and indexes."""
        p = self.output_dir / "work"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def stage_dir(self, n: int) -> Path:
        p = self.results_dir / f"stage{n}"
        p.mkdir(parents=True, exist_ok=True)
        return p


@dataclass
class RunContext:
    config: "FunlrConfig"
    paths: RunPaths
    runner: "CommandRunner"
    state: object | None
    logger: logging.Logger
    dry_run: bool = False

    @property
    def cpu_threads(self) -> int:
        return int(self.config.get("execution", "cpu_threads", 1))

    @property
    def sample_id(self) -> str:
        return str(self.config.get("sample", "sample_id") or "sample")

    def stage_logger(self, n: int, name: str) -> logging.Logger:
        """The parent runner captures this worker's stream in logs/stageN.log."""
        return self.logger
