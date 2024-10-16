from pathlib import Path
from typing import Dict, List
from octopusv.utils.svcf_parser import SVCFFileEventCreator
from octopusv.merger.sv_merger import SVMerger
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from .bench_utils import calculate_metrics, write_vcf, write_summary


class SVBencher:
    def __init__(self, truth_file: Path, call_file: Path, output_dir: Path, **kwargs):
        self.truth_file = truth_file
        self.call_file = call_file
        self.output_dir = output_dir
        self.kwargs = kwargs
        self.truth_events = None
        self.call_events = None
        self.results = None

    def run_benchmark(self):
        self._parse_files()
        self._compare_events()
        self._write_results()

    def _parse_files(self):
        truth_parser = SVCFFileEventCreator([str(self.truth_file)])
        call_parser = SVCFFileEventCreator([str(self.call_file)])
        truth_parser.parse()
        call_parser.parse()

        truth_classifier = SVClassifierByType(truth_parser.events)
        call_classifier = SVClassifierByType(call_parser.events)
        truth_classifier.classify()
        call_classifier.classify()

        truth_chromosome_classifier = SVClassifiedByChromosome(truth_classifier.get_classified_events())
        call_chromosome_classifier = SVClassifiedByChromosome(call_classifier.get_classified_events())
        truth_chromosome_classifier.classify()
        call_chromosome_classifier.classify()

        self.truth_events = truth_chromosome_classifier.get_classified_events()
        self.call_events = call_chromosome_classifier.get_classified_events()

    def _compare_events(self):
        sv_merger = SVMerger(self.truth_events, **self.kwargs)
        sv_merger.merge()

        tp_base, tp_call, fp, fn = [], [], [], []

        for sv_type, chromosomes in self.call_events.items():
            if sv_type == 'TRA':
                for (chr1, chr2), events in chromosomes.items():
                    for event in events:
                        matched_events = sv_merger.get_events(sv_type, chr1, event.start_pos, event.end_pos)
                        if matched_events:
                            tp_call.append(('TRA', event.start_chrom, event.end_chrom, event.start_pos, event.end_pos,
                                            event.source_file, event.bnd_pattern))
                            tp_base.extend(matched_events)
                        else:
                            fp.append(('TRA', event.start_chrom, event.end_chrom, event.start_pos, event.end_pos,
                                    event.source_file, event.bnd_pattern))
            else:
                for chromosome, events in chromosomes.items():
                    for event in events:
                        matched_events = sv_merger.get_events(sv_type, chromosome, event.start_pos, event.end_pos)
                        if matched_events:
                            tp_call.append(event)
                            tp_base.extend(matched_events)
                        else:
                            fp.append(event)

        for sv_type, chromosomes in self.truth_events.items():
            if sv_type == 'TRA':
                for (chr1, chr2), events in chromosomes.items():
                    for event in events:
                        if event not in tp_base:
                            fn.append(('TRA', event.start_chrom, event.end_chrom, event.start_pos, event.end_pos,
                                    event.source_file, event.bnd_pattern))
            else:
                for chromosome, events in chromosomes.items():
                    for event in events:
                        if event not in tp_base:
                            fn.append(event)

        self.results = {
            "tp_base": tp_base,
            "tp_call": tp_call,
            "fp": fp,
            "fn": fn
        }

    def _write_results(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)

        write_vcf(self.output_dir / "tp-base.vcf", self.results["tp_base"])
        write_vcf(self.output_dir / "tp-call.vcf", self.results["tp_call"])
        write_vcf(self.output_dir / "fp.vcf", self.results["fp"])
        write_vcf(self.output_dir / "fn.vcf", self.results["fn"])

        metrics = calculate_metrics(self.results)
        write_summary(self.output_dir / "summary.json", metrics)