#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re

def main():
    parser = argparse.ArgumentParser(description='SVmerger tool')
    parser.add_argument('mode', type=str, help='Mode of operation')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output converted VCF file')

    args = parser.parse_args()

    if args.mode == 'convert':
        # Parse the input VCF file
        headers, events = parse_vcf(args.input)
        # Initialize the EventTransformer with a list of transform strategies
        transformer = EventTransformer([BND_to_INV_Converter(),BND_to_DUP_Converter(), BND_to_TRA_Forward_Converter(), BND_to_TRA_Reverse_Converter()])
        # Apply all transformation strategies to the events
        transformed_events = transformer.apply_transforms(events)
        # Write the transformed events to the output file
        transformer.write_vcf(headers, transformed_events, args.output)
    else:
        print(f'Unknown mode: {args.mode}')
        parser.print_help()

def check_vcf_format(vcf_file_path):
    """
    Check the format of a VCF file.
    Raise an error and exit if the format is incorrect.
    """
    with open(vcf_file_path) as f:
        lines = f.readlines()

    # Check if there is at least one header line
    if not any(line.startswith('#') for line in lines):
        print("ERROR: Invalid VCF format. The file should contain at least one header line starting with '#'.")
        exit(1)

    for line in lines:
        if line.startswith('#'):
            continue  # Skip header lines

        # Check for space in lines
        if ' ' in line:
            print("ERROR: Invalid VCF format. Non-header lines should not contain spaces.")
            exit(1)
        
        fields = line.split('\t')

        # Check the number of columns
        if len(fields) < 10:
            print(f"ERROR: Invalid VCF format. Expected at least 10 fields, but got {len(fields)}")
            exit(1)

        # Check that the position is a number
        try:
            pos = int(fields[1])
        except ValueError:
            print(f"ERROR: Invalid VCF format. Position (field 2) should be a number, but got {fields[1]}")
            exit(1)

        # Check that the quality score is a number or '.'
        if fields[5] != '.':
            try:
                qual = float(fields[5])
            except ValueError:
                print(f"ERROR: Invalid VCF format. Quality score (field 6) should be a number or '.', but got {fields[5]}")
                exit(1)

def parse_vcf(vcf_file_path):
    """
    Parse VCF file into a list of SVEvent objects and return headers.
    """
    check_vcf_format(vcf_file_path) # Check the format first
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

    return headers, events # A list of SV objects
        
class SVEvent: # A class to represent each SV event
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format, sample):
        self.chrom = chrom
        self.pos = int(pos)
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

    def __str__(self):
        # enable you to output object as vcf string format
        # Define your desired order
        ordered_keys = ['SVTYPE', 'END', 'SVLEN', 'SUPPORT', 'COVERAGE', 'STRAND']
        info_items = []
        
        # First add the ordered keys
        for key in ordered_keys:
            if key in self.info:
                info_items.append(f'{key}={self.info[key]}')

        # Then add the remaining keys
        for key in self.info:
            if key not in ordered_keys:
                info_items.append(f'{key}={self.info[key]}')

        info_str = ';'.join(info_items)  # SVTYPE=INV;END=69650961;SVLEN=1835141

        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            self.chrom,
            self.pos,
            self.id,
            self.ref,
            self.alt,
            self.qual,
            self.filter,
            info_str, 
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

class Converter:
    """
    This is an abstract base class for all converter classes. It provides a common interface for all converters.
    The `convert` method is a placeholder that needs to be overridden in each concrete converter class.
    """
    def convert(self, event):
        """
        This method should be overridden in a subclass. 
        """
        raise NotImplementedError

class BND_to_INV_Converter(Converter):
    """
    This class inherits from the `Converter` base class and implements the conversion logic for BND to INV conversion.
    """
    def convert(self, event):  # Each converter has a convert function
        try:
            if event.is_BND():
                pattern = get_BND_pattern(event.alt)
                if pattern == 't]p]' or pattern == '[p[t':
                    split_result = re.split(r'[\[\]:]', event.alt) # Use regex to split and get the chrom and pos
                    if len(split_result) != 4:
                        print(f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}")
                    else:
                        chrom_alt, pos_alt = split_result[1:3]
                        if event.chrom == chrom_alt: # Do this only when same chromosome 
                            pos_alt = int(pos_alt)
                            if pattern == 't]p]' and event.pos < pos_alt:
                                end = pos_alt
                                svlen = abs(event.pos - pos_alt)
                                self.convert_to_INV(event, end, svlen) # self will be passed to convert_to_INV method.
                            elif pattern == '[p[t' and event.pos > pos_alt:
                                end = event.pos
                                event.pos = pos_alt
                                svlen = abs(event.pos - end)
                                self.convert_to_INV(event, end, svlen) 

        except Exception as e:
           print("Failed to convert BND to Inversion: ", e)
           
    def convert_to_INV(self, event, end, svlen):
        # Detail logic changing BND event to INV event
        event.alt = '<INV>'
        event.info['SVTYPE'] = 'INV'
        event.info['END'] = end
        event.info['SVLEN'] = svlen

class BND_to_DUP_Converter(Converter):
    """
    This class inherits from the `Converter` base class and implements the conversion logic for BND to DUP conversion.
    """
    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_BND_pattern(event.alt)
                if pattern in [']p]t', 't[p[']:
                    split_result = re.split(r'[\[\]:]', event.alt)
                    if len(split_result) != 4:
                        print(f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}")
                    else:
                        chrom_alt, pos_alt = split_result[1:3]
                        if event.chrom == chrom_alt:
                            pos_alt = int(pos_alt)
                            if pattern == ']p]t' and event.pos < pos_alt:
                                end = pos_alt
                                svlen = abs(end - event.pos)
                                self.convert_to_DUP(event, end, svlen)
                            elif pattern == 't[p[' and event.pos > pos_alt:
                                end = event.pos
                                svlen = abs(event.pos - pos_alt)
                                event.pos = pos_alt
                                self.convert_to_DUP(event, end, svlen)
        except Exception as e:
           print("Failed to convert BND to Duplication: ", e)

    def convert_to_DUP(self, event, end, svlen):
        event.alt = '<DUP>'
        event.info['SVTYPE'] = 'DUP'
        event.info['END'] = end
        event.info['SVLEN'] = svlen

class BND_to_TRA_Forward_Converter(Converter):
    """
    This class inherits from the `Converter` base class and implements 
    the conversion logic for BND to TRA Forward translocation cut-paste.
    """
    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_BND_pattern(event.alt)
                if pattern in ['t[p[', ']p]t']:
                    split_result = re.split(r'[\[\]:]', event.alt)
                    if len(split_result) != 4:
                        print(f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}")
                    else:
                        chrom_alt, pos_alt = split_result[1:3]
                        if event.chrom == chrom_alt:
                            pos_alt = int(pos_alt)
                            if pattern == 't[p[' and event.pos < pos_alt:
                                end = pos_alt
                                self.convert_to_TRA_forward(event, end)
                            elif pattern == ']p]t' and event.pos > pos_alt:
                                end = pos_alt
                                self.convert_to_TRA_forward(event, end)
        except Exception as e:
           print("Failed to convert BND to Translocation: ", e)

    def convert_to_TRA_forward(self, event, end):
        event.info['SVTYPE'] = 'TRA'
        event.info['END'] = end
        event.info['SVLEN'] = 0

class BND_to_TRA_Reverse_Converter(Converter):
    """
    This class inherits from the `Converter` base class and implements the conversion logic for BND to TRA reverse conversion.
    """
    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_BND_pattern(event.alt)
                if pattern in ['[p[t', 't]p]']:
                    split_result = re.split(r'[\[\]:]', event.alt)
                    if len(split_result) != 4:
                        print(f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}")
                    else:
                        chrom_alt, pos_alt = split_result[1:3]
                        if event.chrom == chrom_alt:
                            pos_alt = int(pos_alt)
                            if pattern == 't]p]' and event.pos > pos_alt:
                                end = pos_alt
                                self.convert_to_TRA_reverse(event, end)
                            elif pattern == '[p[t' and event.pos < pos_alt:
                                end = pos_alt
                                self.convert_to_TRA_reverse(event, end)
        except Exception as e:
            print("Failed to convert BND to Translocation: ", e)

    def convert_to_TRA_reverse(self, event, end):
        event.info['SVTYPE'] = 'TRA'
        event.info['END'] = end
        event.info['SVLEN'] = 0
 
# You can add more converter classes here...

class EventTransformer:  # The input is lists.
    # The EventTransformer class manages transforming events and output together.
    # It is initialized with a list of transform strategies,
    def __init__(self, transform_strategies):
        # Initialize the EventTransformer with a list of transform strategies.
        self.transform_strategies = transform_strategies

    def apply_transforms(self, events):
        # Apply all transformation strategies to a list of events.
        for event in events:  # Try every converters for each event.
            for strategy in self.transform_strategies:
                strategy.convert(event) # Strategy is a converter instance, like BND_to_INV_Converter 
        return events    # Returns: The transformed list of events.

    def write_vcf(self, headers, events, output_file):
        # Write the transformed events to a VCF file.
        with open(output_file, 'w') as f:
            for header in headers:
                f.write(header + '\n')
            for event in events:
                f.write(str(event) + '\n')


if __name__ == '__main__':
    main()

# 每个软件内部要校正SVLEN的计算方式
# 每个软件转换后要内部去冗余
# 我的标准格式：INV用<INV>,DUP 用<DUP> INS, DEL最好保留真实序列