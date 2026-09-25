"""Internal isolated stage entry point; users run the validated public CLI."""
import importlib
import json
import logging
from pathlib import Path
import sys

from .config import FunlrConfig
from .core.context import RunContext, RunPaths
from .core.runner import CommandRunner
from .runner import STAGES


def main():
    number = int(sys.argv[1])
    if not 0 <= number < len(STAGES):
        raise ValueError("Stage must be between 0 and 7")
    settings = json.loads(Path(sys.argv[2]).read_text())
    config = FunlrConfig(settings["config"])
    root = Path(config.get("execution", "output_dir"))
    logger = logging.getLogger(f"funlr.stage{number}")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    commands = CommandRunner(config["tools"], root / "logs/commands.jsonl", logger,
                             inventory=settings.get("inventory"))
    ctx = RunContext(config, RunPaths(root), commands, None, logger, False)
    module = importlib.import_module("funlr.stages." + STAGES[number])
    try:
        module.run(ctx)
    except Exception:
        logger.exception("Stage %s failed", number)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
