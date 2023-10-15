# The Structural Variant Call Format (SVCF) Version 1.0 Specification

## 1. Introduction

The Structural Variant Call Format (SVCF) is a text file format for representing structural variant (SV) calls. It is based on the [Variant Call Format (VCF) Version 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf).  

SVCF shares much of the same structure as VCF, with meta-information lines, a header line, and then data lines containing information about structural variants. The key differences lie in the representation of the structural variants and additional optional fields to provide more details about each variant. One of the main goals of SVCF is to provide a more intuitive, unified format that can be seamlessly integrated into the SVoctopus software for various operations including merging, statistical analysis, and benchmarking.

## 2. An Example

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1234567	.	N	<DUP>	.	PASS	SVTYPE=DUP;SVLEN=500	GT	0/1
```

## 3. Meta-Information Lines

Same as VCF 4.2. Provides metadata about the file and declarations for the format and info fields.

## 4. Header Line 

The header line names the 10 mandatory columns:

```
#CHROM  POS	ID	REF	ALT QUAL    FILTER	INFO    FORMAT Sample
``` 

Additional optional columns may follow.

## 5. Data Lines

**5.1 CHROM** - Chromosome name 

**5.2 POS** - Position of the variant (BP). For SVs, this may be the approximate midpoint.

**5.3 ID** - Semi-colon separated list of unique identifiers.

**5.4 REF** - The reference allele. For SVs, this may be `<DUP>` for duplication, `<DEL>` for deletion, etc.

**5.5 ALT** - The alternative allele(s). For SVs, this describes the structural variant.

**5.6 QUAL** - Phred-scaled quality score (-10log10 prob of call being wrong). "." if unknown.

**5.7 FILTER** - PASS if this position has passed filters. Otherwise lists failed filters.

**5.8 INFO** - Additional variant info. Standard fields plus custom SV annotations.

**5.9 FORMAT** - Sample names, genotypes, and supporting data.

**5.10 Sample** - 



## 6. Implementation Notes

- Missing values should be denoted with a dot '.'
- Multiple values in the ALT field should be comma separated.

This covers the key points of the SVCF format. Further implementation details are to be determined. The goal is to build upon the VCF standard to represent structural variants in a scalable and extensible manner.