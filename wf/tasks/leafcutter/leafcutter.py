from pathlib import Path
from typing import Annotated, List, Optional

from flytekit.core.annotation import FlyteAnnotation
from latch import medium_task, message
from latch.types import LatchDir, LatchFile
from wf.core.models import TrimgaloreSalmonOutput
from wf.core.utils import remote_output_dir, run


@medium_task
def leafcutter(
    run_name: str,
    ts_outputs: List[Optional[TrimgaloreSalmonOutput]],
    output_directory: Optional[LatchDir],
    run_splicing: bool = False,
    manual_conditions: Annotated[
        List[List[str]],
        FlyteAnnotation({"_tmp_hack_deseq2": "manual_design_matrix"}),
    ] = [],
) -> (LatchFile, LatchFile, LatchFile):

    if run_splicing is False:
        # random file noop, hack until boolean conditionals work
        return (
            LatchFile("/root/wf/__init__.py"),
            LatchFile("/root/wf/__init__.py"),
            LatchFile("/root/wf/__init__.py"),
        )

    REMOTE_PATH = f"latch:///{remote_output_dir(output_directory)}{run_name}/Alternative Splicing (LeafCutter)/"
    """Remote path prefix for LatchFiles + LatchDirs"""

    message(
        "info",
        {
            "title": "Generating alternative splicing junction file from all samples",
            "body": "\n".join(f"- {x.sample_name}" for x in ts_outputs),
        },
    )

    all_juncfiles = Path("/root/juncs.txt")
    sample_names = []
    with open(all_juncfiles, "w") as f:
        for tso in ts_outputs:
            junc_file = Path(tso.junction_file.local_path)
            sample_name = junc_file.name.replace(" ", "")
            sanitized_junc_file = Path(tso.junction_file.local_path).rename(sample_name)
            sample_names.append(sample_name)
            print(f"local: {sanitized_junc_file.resolve()}")
            f.write(f"{sanitized_junc_file.resolve()}\n")

    cluster_counts = Path(f"/root/{run_name}_perind_numers.counts.gz")
    run(
        [
            "python",
            "/root/leafcutter/clustering/leafcutter_cluster_regtools.py",
            "-j",
            str(all_juncfiles),
            "-m",
            "50",
            "-o",
            run_name,
            "-l",
            "5000000",
            "-k=True",  # Don't error on weird chromosome names
        ]
    )

    groups = Path("/root/groups.txt")
    existing_conds = {}
    with open(groups, "w") as f:
        for cond in manual_conditions:
            if cond[1] == "__ignore":
                continue
            if cond[1] not in existing_conds:
                if len(existing_conds) > 2:
                    print(
                        "Only running differential splicing for first two"
                        f" conditions: {list(existing_conds.keys())}"
                    )
                    break
                else:
                    existing_conds[cond[1]] = True

            f.write(f"{cond[0].replace(' ', '')}.bam {cond[1]}\n")

    run(
        [
            "/root/leafcutter/scripts/leafcutter_ds.R",
            "--num_threads",
            "30",
            str(cluster_counts),
            str(groups),
            "--min_samples_per_intron=1",
            "--min_samples_per_group=1",
        ]
    )

    cluster_sig = Path("/root/leafcutter_ds_cluster_significance.txt")
    cluster_es = Path("/root/leafcutter_ds_effect_sizes.txt")

    return (
        LatchFile(cluster_counts, REMOTE_PATH + cluster_counts.name),
        LatchFile(cluster_sig, REMOTE_PATH + cluster_sig.name),
        LatchFile(cluster_es, REMOTE_PATH + cluster_es.name),
    )
