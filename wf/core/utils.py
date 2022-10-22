import shutil
import subprocess
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

from latch.types import LatchDir


def run(command: List[str], check: bool = True, capture_output: bool = False):
    return subprocess.run(command, check=check, capture_output=capture_output)


def remote_output_dir(custom_output_dir: Optional[LatchDir]) -> str:
    if custom_output_dir is None:
        return "/RNA-Seq Outputs/"
    remote_path = custom_output_dir.remote_path
    assert remote_path is not None
    if remote_path[-1] != "/":
        remote_path += "/"
    if remote_path[:8] == "latch://":
        remote_path = remote_path[8:]
    return remote_path


def unzip_if_needed(path: Path) -> Path:

    is_gzipped = str(path).endswith(".gz")
    if not is_gzipped:
        return path

    run(["gunzip", str(path)])
    return Path(str(path).removesuffix(".gz")).resolve()


def run_and_capture_output(command: List[str]) -> Tuple[int, str]:
    captured_stdout = []

    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True,
    ) as process:
        assert process.stdout is not None
        for line in process.stdout:
            print(line)
            captured_stdout.append(line)
        process.wait()
        returncode = process.returncode

    return returncode, "\n".join(captured_stdout)


def sanitize(value: str) -> str:
    return value.lower().replace(" ", "_")


def concatenate_files(filepaths: Iterable[str], output_path: str) -> Path:
    path = Path(output_path).resolve()
    with path.open("w") as output_file:
        for p in filepaths:
            p = p.removesuffix(".gz")
            with open(p, "r") as f:
                shutil.copyfileobj(f, output_file)
    return path
