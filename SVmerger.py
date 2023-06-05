#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
def main():
    parser = argparse.ArgumentParser(description='SVmerger tool')
    parser.add_argument('mode', type=str, help='Mode of operation')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output converted VCF file')

    args = parser.parse_args()

    if args.mode == 'convert':
        convert_BND_to_INV(args.input, args.output)
    else:
        print(f'Unknown mode: {args.mode}')
        parser.print_help()

def parse_vcf(vcf_file_path):
    """
    Parse VCF file into a list of SVEvent objects and return headers.
    """
    events = []
    headers = []

    with open(vcf_file_path) as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line.strip())
                continue  # Skip header lines

            fields = line.strip().split('\t')
            event = SVEvent(*fields) # Unpack fields and send to SVEvent class
            events.append(event)

    return events # A list of SV objects
        
class SVEvent: # A class to represent each SV event
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format, sample):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = self._parse_info(info) # self.info will be a dict
        self.format = format
        self.sample = sample

    def _parse_info(self, info): # Parse the info field according to your requirement
        info_dict = {}
        for item in info.split(';'):
            parts = item.split('=')
            if len(parts) == 2: # Just in case that the INFO doesn't contain = symbol
                key, value = parts
                info_dict[key] = value
            else:
                info_dict[parts[0]] = None
        return info_dict

    def is_duplication(self):
        # Whether is duplication
        return self.info['SVTYPE'] == 'DUP'

    def is_inversion(self):
        # Whether is inversion
        return self.info['SVTYPE'] == 'INV'
    
    def is_insertion(self):
        #  Whether is insertion
        return self.info['SVTYPE'] == 'INS' 

    def is_TRA(self):
        #  Whether is TRA
        return self.info['SVTYPE'] == 'TRA' 

    def is_BND(self):
        #  Whether is BND
        return self.info['SVTYPE'] == 'BND' 

    def convert_to_INV(self, end, svlen):
        # Change BND event to INV event
        self.alt = '<INV>'
        self.info['SVTYPE'] = 'INV'
        self.info['END'] = end
        self.info['SVLEN'] = svlen

    def __str__(self):
        # enable you to output object as vcf string format
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            self.chrom,
            self.pos,
            self.id,
            self.ref,
            self.alt,
            self.qual,
            self.filter,
            ';'.join(f'{k}={v}' for k, v in self.info.items()),  # SVTYPE=INV;END=69650961;SVLEN=1835141
            self.format,
            self.sample
        )

def get_BND_pattern(alt):
    # Get the pattern of BND: t[p[, t]p], ]p]t, [p[t
    if alt[0] in 'ATCGN' and alt[1] == '[':
        return 't[p['
    elif alt[0] in 'ATCGN' and alt[1] == ']':
        return 't]p]'
    elif alt[-1] in 'ATCGN' and alt[-2] == '[':
        return '[p[t'
    elif alt[-1] in 'ATCGN' and alt[-2] == ']':
        return ']p]t'
    else:
        return None
 
def convert_BND_to_INV(input_file, output_file):
    """
    convert on-standard BND to standard inversion
    """
    try:
        headers, events = parse_vcf(input_file)
        
        with open(output_file, 'w') as f:
            for header in headers:  # Write headers to output file first
                f.write(header + '\n')

            for event in events:
                if event.is_BND():
                    pattern = event.get_BND_pattern()
                    if pattern == 't]p]' or pattern == '[p[t':
                        chrom_alt, pos_alt = re.split(r'[\[\]:]', event.alt)[1:3]  # Use regex to split and get the chrom and pos
                        if event.chrom == chrom_alt: # Do this only when same chromosome 
                            pos_alt = int(pos_alt)
                            if pattern == 't]p]' and event.pos < pos_alt:
                                end = pos_alt
                                svlen = abs(event.pos - pos_alt)
                                event.convert_to_INV(end, svlen)
                            elif pattern == '[p[t' and event.pos > pos_alt:
                                end = event.pos
                                event.pos = pos_alt
                                svlen = abs(event.pos - end)
                                event.convert_to_INV(end, svlen)

                f.write(str(event) + '\n')
    
    except Exception as e:
        print("Failed to convert BND to Inversion: ", e)


if __name__ == '__main__':
    main()


下一步，写BND to TRA, BND to DUP的函数，解决顺序问题