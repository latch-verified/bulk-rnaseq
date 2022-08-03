from pathlib import Path
from typing import List

import matplotlib.pyplot as plt

SUB_SAMPLE_SIZE = 100


def bs(target, lst):
    first, last = 0, len(lst)
    while first < last:
        mid = (first + last) // 2
        if lst[mid] < target:
            first = mid + 1
        else:
            last = mid
    return first


def gtf_to_gbc(gtf_path: Path, quant_paths: List[Path]):
    ranges_by_transcript_name: dict[str, List] = {}

    total_range = [float("inf"), float("-inf")]

    with open(gtf_path, "r") as g:
        for line in g:
            chrom, _, typ, start, end, _, stran, _, info = line.strip(" \n").split("\t")
            start = int(start)
            end = int(end)
            if start < total_range[0]:
                total_range[0] = start
            if end > total_range[1]:
                total_range[1] = end
            transcript_name = (
                info.split(";")[3].removeprefix(' transcript_name "').removesuffix('"')
            )
            if transcript_name not in ranges_by_transcript_name:
                ranges_by_transcript_name[transcript_name] = []
            ranges_by_transcript_name[transcript_name].append((start, end))

    total_length = total_range[1] - total_range[0]
    if total_length < 0:
        raise ValueError("you done goofed")

    num_bins = min(total_length, SUB_SAMPLE_SIZE)
    bin_size = total_length // num_bins
    excess = total_length % num_bins

    bin_sizes = [bin_size] * num_bins
    for i in range(excess):
        bin_sizes[i] += 1

    partial_sums = [0]
    for bin_size in bin_sizes:
        partial_sums.append(partial_sums[-1] + bin_size)

    for transcript_name, lst in ranges_by_transcript_name.items():
        for i, (start, end) in enumerate(lst):
            lst[i] = bs(start, partial_sums), bs(end, partial_sums)

    gbcs = []
    names = []

    for quant_path in quant_paths:
        gbc = [0 for _ in range(num_bins + 3)]
        with open(quant_path, "r") as q:
            line_no = -1
            for line in q:
                line_no += 1
                if line_no == 0:
                    continue
                transcript_name, _, _, _, num_reads = line.strip(" \n").split("\t")
                num_reads = float(num_reads)
                if transcript_name not in ranges_by_transcript_name:
                    raise ValueError("Gene in Quant File that isn't in GTF File")
                for start, end in ranges_by_transcript_name[transcript_name]:
                    gbc[start] += num_reads
                    gbc[end + 1] -= num_reads

        for i, g in enumerate(gbc):
            curr_val = gbc[i - 1] if i > 0 else 1
            gbc[i] = curr_val + g

        gbcs.append(gbc)
        names.append(quant_path.parent.name)

    plot(gbcs, names, fig_name="gene_body_coverage.png")


def plot(lsts: List[List], names: List[str], fig_name: str):

    import numpy as np
    from scipy.interpolate import make_interp_spline

    for lst, name in zip(lsts, names):
        spline = make_interp_spline(np.arange(0, len(lst)), np.array(lst))

        x = np.linspace(0, len(lst), len(lst) * 10)
        y = spline(x)
        fig = plt.Figure()
        plt.plot(x, y, label=name)

    plt.yscale("log")
    plt.xlabel("5' -> 3' Percentiles")
    plt.ylabel("Coverage")
    plt.legend()
    plt.savefig(fig_name)
    plt.show()


if __name__ == "__main__":
    gtf_to_gbc(
        Path("Genome/phaffii_WT.c2.noTrinity.gtf"),
        [
            Path(
                "Quantification (salmon)/220506LovA_D22_174010/220506LovA_D22_174010_quant.sf"
            ),
            Path(
                "Quantification (salmon)/220506LovA_D22_174001/220506LovA_D22_174001_quant.sf"
            ),
        ],
    )
