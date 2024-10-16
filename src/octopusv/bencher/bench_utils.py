import json
from pathlib import Path
from typing import Dict, List, Union


def calculate_metrics(results: Dict[str, List]):
    tp = len(results["tp_call"])
    fp = len(results["fp"])
    fn = len(results["fn"])

    precision = tp / (tp + fp) if tp + fp > 0 else 0
    recall = tp / (tp + fn) if tp + fn > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "Precision": precision,
        "Recall": recall,
        "F1": f1
    }


def write_vcf(file_path: Path, events: List[Union[tuple, object]]):
    with file_path.open('w') as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for event in events:
            if isinstance(event, tuple):  # TRA event
                chrom, start_chrom, end_chrom, start_pos, end_pos, source_file, bnd_pattern = event
                f.write(f"{start_chrom}\t{start_pos}\t.\tN\t{bnd_pattern}\t.\tPASS\t"
                        f"SVTYPE=TRA;END={end_pos};CHR2={end_chrom};SOURCES={source_file}\n")
            else:  # Other SV events
                f.write(f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
                        f"{event.quality}\t{event.filter}\t{';'.join(f'{k}={v}' for k, v in event.info.items())}\n")


def write_summary(file_path: Path, metrics: Dict):
    with file_path.open('w') as f:
        json.dump(metrics, f, indent=2)