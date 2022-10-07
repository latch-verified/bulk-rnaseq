import subprocess
from typing import List, Optional

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
