import csv
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Tuple

from latch import medium_task, message
from latch.types import LatchDir, LatchFile
from wf.core.models import TrimgaloreSalmonOutput
from wf.core.utils import remote_output_dir

_COUNT_TABLE_GENE_ID_COLUMN = "gene_id"


@medium_task
def count_matrix_and_multiqc(
    run_name: str,
    ts_outputs: List[Optional[TrimgaloreSalmonOutput]],
    output_directory: Optional[LatchDir],
) -> Tuple[Optional[LatchFile], Optional[LatchFile]]:

    count_matrix_file = None
    multiqc_report_file = None

    REMOTE_PATH = f"latch:///{remote_output_dir(output_directory)}{run_name}/"
    """Remote path prefix for LatchFiles + LatchDirs"""

    ts_outputs = [x for x in ts_outputs if x and x.passed_tximport]
    message(
        "info",
        {
            "title": "Generating count matrix from successful samples",
            "body": "\n".join(f"- {x.sample_name}" for x in ts_outputs),
        },
    )

    combined_counts = defaultdict(dict)
    for tso in ts_outputs:
        gene_abundance_file = Path(
            tso.gene_abundance_file.local_path).resolve()
        with gene_abundance_file.open("r") as f:
            for row in csv.DictReader(f, dialect=csv.excel_tab):
                gene_name = row["Name"]
                combined_counts[gene_name][tso.sample_name] = row["NumReads"]

    raw_count_table_path = Path("./counts.tsv").resolve()
    with raw_count_table_path.open("w") as file:
        sample_names = (x.sample_name for x in ts_outputs)
        writer = csv.DictWriter(
            file,
            [_COUNT_TABLE_GENE_ID_COLUMN, *sample_names],
            delimiter="\t",
        )
        writer.writeheader()
        for gene_id, data in combined_counts.items():
            data[_COUNT_TABLE_GENE_ID_COLUMN] = gene_id
            writer.writerow(data)

    count_matrix_file = LatchFile(
        str(raw_count_table_path),
        REMOTE_PATH + "Quantification (salmon)/counts.tsv",
    )

    try:
        aux_paths = [
            str(Path(x.salmon_aux_output.local_path).resolve()) for x in ts_outputs
        ]
        subprocess.run(["multiqc", *aux_paths], check=True)
        multiqc_report_file = LatchFile(
            "/root/multiqc_report.html",
            REMOTE_PATH + "multiqc_report.html",
        )
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while generating MultiQC report -> {e}")
        message(
            "error",
            {
                "title": "Unable to generate MultiQC report",
                "body": "See logs for more information",
            },
        )

    return count_matrix_file, multiqc_report_file
