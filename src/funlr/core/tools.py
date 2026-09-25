"""Database helpers never write beside user-supplied resources."""
import shutil
from pathlib import Path

from .errors import StageError

PRESSED_EXTS = ("h3m", "h3i", "h3f", "h3p")


def is_pressed(path):
    return all(Path(f"{path}.{ext}").is_file() and Path(f"{path}.{ext}").stat().st_size > 0 for ext in PRESSED_EXTS)


def ensure_pressed(hmm_path, workdir, runner, logger):
    from funlr.runner import sha256
    source, work = Path(hmm_path).resolve(), Path(workdir).resolve()
    if source.is_relative_to(work) and is_pressed(source):
        return source
    work.mkdir(parents=True, exist_ok=True)
    digest = sha256(source)
    target = work / (digest + ".hmm")
    if not target.is_file() or sha256(target) != digest:
        shutil.copyfile(source, target)
        for ext in PRESSED_EXTS:
            Path(f"{target}.{ext}").unlink(missing_ok=True)
    if not is_pressed(target):
        runner.run([runner.resolve("hmmpress"), "-f", str(target)], check=True)
    if not is_pressed(target):
        raise StageError(f"hmmpress did not create nonempty indexes for {target}")
    return target


def require_tools(runner, tool_keys, logger):
    missing = [key for key in tool_keys if not runner.has_tool(key)]
    if missing:
        raise StageError("Missing required external tools: " + ", ".join(missing))
    return {key: runner.tool_version(key) for key in tool_keys}
