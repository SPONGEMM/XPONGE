"""Helpers for launching the exact installed/source Xponge worker package."""

from __future__ import annotations

import os
from pathlib import Path
import signal
import subprocess
import sys
from typing import Sequence


def worker_command(module: str, *, python_executable: str | None = None) -> list[str]:
    """Return a subprocess command that prefers this Xponge source tree.

    Editable installs can register a meta-path redirect ahead of normal path
    lookup.  A parent process may still have imported this source checkout
    explicitly, so workers must not silently jump to a different Xponge tree.
    Moving ``PathFinder`` first keeps normal installed dependencies available
    while binding the worker to the package that launched it.
    """

    source_root = str(Path(__file__).resolve().parents[2])
    bootstrap = (
        "import runpy,sys;"
        "from importlib.machinery import PathFinder;"
        "sys.meta_path=[item for item in sys.meta_path if item is not PathFinder];"
        "sys.meta_path.insert(0,PathFinder);"
        f"sys.path.insert(0,{source_root!r});"
        f"runpy.run_module({module!r},run_name='__main__')"
    )
    executable = sys.executable
    if python_executable is not None:
        # Preserve virtual-environment launcher symlinks. Resolving the symlink
        # would execute the base interpreter and silently lose that venv's QM
        # packages.
        candidate = Path(python_executable).expanduser().absolute()
        if not candidate.is_file():
            raise FileNotFoundError(f"worker Python executable does not exist: {candidate}")
        executable = str(candidate)
    return [executable, "-c", bootstrap]


def run_worker_subprocess(
    command: Sequence[str],
    *,
    input_text: str,
    timeout_seconds: float,
) -> subprocess.CompletedProcess[str]:
    """Run one isolated worker and reap its whole POSIX process group on timeout."""

    kwargs = {
        "stdin": subprocess.PIPE,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    if os.name == "posix":
        kwargs["start_new_session"] = True
    process = subprocess.Popen(list(command), **kwargs)
    try:
        stdout, stderr = process.communicate(input=input_text, timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            except OSError:
                process.kill()
        else:
            process.kill()
        process.communicate()
        raise
    return subprocess.CompletedProcess(list(command), process.returncode, stdout, stderr)


__all__ = ["run_worker_subprocess", "worker_command"]
