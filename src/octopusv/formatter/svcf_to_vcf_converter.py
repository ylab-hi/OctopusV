class SVCFtoVCFConverter:
    def __init__(self, events, input_svcf_file):
        self.events = events
        self.input_svcf_file = input_svcf_file

    def convert(self):
        vcf_content = self._generate_vcf_header()
        for event in self.events:
            vcf_content += self._convert_event_to_vcf(event)
        return vcf_content

    def _generate_vcf_header(self):
        # read contig info from SVCF
        contig_lines = ""
        with open(self.input_svcf_file) as f:
            for line in f:
                if line.startswith("##contig"):
                    contig_lines += line
                elif line.startswith("#CHROM"):
                    break

        return f"""##fileformat=VCFv4.2
{contig_lines}##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end coordinate">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting this variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method used to detect SV">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"""

    def _convert_event_to_vcf(self, event):
        chrom = event.chrom
        pos = event.pos
        id = event.sv_id
        ref = event.ref if event.ref != "N" else "."
        alt = self._get_alt(event)
        qual = event.quality if hasattr(event, "quality") else "."
        filter = event.filter if hasattr(event, "filter") else "PASS"

        # 调整插入的 END
        end_pos = event.pos if event.sv_type == "INS" else event.end_pos

        info_fields = [f"SVTYPE={event.sv_type}", f"END={end_pos}"]

        # 调整删除和插入的 SVLEN
        if "SVLEN" in event.info:
            svlen_str = event.info["SVLEN"]
            try:
                svlen = int(svlen_str.strip())
                if event.sv_type == "DEL":
                    svlen = -abs(svlen)
                elif event.sv_type == "INS":
                    svlen = abs(svlen)
                info_fields.append(f"SVLEN={svlen}")
            except ValueError:
                pass  # 如果 SVLEN 不是有效的整数，跳过

        if "CHR2" in event.info:
            info_fields.append(f"CHR2={event.info['CHR2']}")
        if "SUPPORT" in event.info:
            info_fields.append(f"SUPPORT={event.info['SUPPORT']}")
        if "SVMETHOD" in event.info:
            info_fields.append(f"SVMETHOD={event.info['SVMETHOD']}")
        if "STRAND" in event.info:
            info_fields.append(f"STRAND={event.info['STRAND']}")

        info = ";".join(info_fields)

        format = "GT:AD:DP:LN"
        gt = event.sample.get("GT", "./.")
        ad = event.sample.get("AD", ".,.")
        dp = self._calculate_dp(ad)
        ln = event.sample.get("LN", ".")

        sample = f"{gt}:{ad}:{dp}:{ln}"

        return f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{sample}\n"

    def _get_alt(self, event):
        if event.alt and event.alt not in ("N", "."):
            return event.alt
        return f"<{event.sv_type}>"

    def _calculate_dp(self, ad):
        try:
            return sum(int(x) for x in ad.split(",") if x != "." and x.strip().isdigit())
        except ValueError:
            return "."
