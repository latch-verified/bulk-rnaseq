import subprocess
from typing import List


def run(command: List[str], check: bool = True, capture_output: bool = False):
    return subprocess.run(command, check=check, capture_output=capture_output)
