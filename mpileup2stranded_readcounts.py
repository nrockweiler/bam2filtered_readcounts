#!/usr/bin/python3
"""Convert mpileup format to stranded readcounts format.  Stranded readcounts are printed to STDOUT.

OUTPUT COLUMNS:
   1  chr = chromosome
   2  pos = position (1-based)
   3  ref_allele = reference allele
   4  alt_allele = alternative allele (non-reference allele with highest coverage; if ties, the first allele when sorted alphabetically is used)
   5  ref_count = total ref coverage
   6  alt_count = total alt coverage
   7  alt_pos_count = alt coverage on + strand
   8  alt_neg_count = alt coverage on - strand
   9  alt_on_both_strands = 1 if alt was observed on + and - strand; 0 otherwise
  10  max_alt_pos = most internal position of the alt allele on the reads, e.g., 2 means the alt allele was observed 2bp from the end.  No other read containing the alt read had an alt base closer to the center of the read.
  11  A_count = A coverage on + strand
  12  a_count = A coverage on - strand
  13  C_count = C coverage on + strand
  14  c_count = C coverage on - strand
  15  G_count = G coverage on + strand
  16  g_count = G coverage on - strand
  17  T_count = T coverage on + strand
  18  t_count = T coverage on - strand
  19  N_count = N coverage on + strand
  20  n_count = N coverage on - strand
  21  sample = Sample ID (useful when concatenating multiple readcount files)

NOTES:
 - mpileup2readcounts and mpileup2stranded_readcounts output different columns.  Check the usage

SEE ALSO:
 - mpileup2readcounts
"""
import os
import sys
import argparse
import re
import numpy as np

# Notes: if unspecified, ref_allele and alt_allele is uppercase, and lc_*_allele is lowercase

def main():

    args = _parse_cmd_line()

    poi = _load_poi(args.variant_bed)

    if args.header:
        _print_header()

    if args.mpileup is None:
        mpileup_fh = sys.stdin
    else:
        mpileup_fh = open(args.mpileup, 'rU')
    for line in mpileup_fh:
        line = line.rstrip('\n')
        cols = line.split('\t')
        chrom = cols[0]
        pos_1based = cols[1]

        pos_key = ':'.join([chrom, pos_1based])

        # If this is not at a position of interest, skip it
        if not (pos_key in poi):
            continue

        ref_allele = cols[2].upper()
        depth = int(cols[3])

        if depth > 0:
            base_str = cols[4]
            #base_quals = cols[5]
            base_pos_str = cols[6]

            read_counts, base_pos_counts = _get_counts(base_str, base_pos_str, ref_allele)
            unstranded_read_counts = _stranded2unstranded(read_counts)
            unstranded_base_pos_counts = _stranded2unstranded(base_pos_counts)
            alt_allele = _get_alt_allele(unstranded_read_counts, ref_allele)

            max_internal_base_pos = _get_max_internal_base_pos(unstranded_base_pos_counts, args.read_length)

            _print_line(chrom, pos_1based, ref_allele, alt_allele, read_counts, unstranded_read_counts, max_internal_base_pos, args.sample)
        else: # There is no coverage here

            _print_line_no_cov(chrom, pos_1based, ref_allele, args.sample)

    mpileup_fh.close()

def _load_poi(bed_fn):

    poi = {} # key = chrom:1-based position; value = 1 (dummy)
    with open(bed_fn, 'rU') as bed_fh:
        for line in bed_fh:
            line = line.rstrip('\n')
            cols = line.split('\t')
            chrom = cols[0]
            # start = cols[1]
            end = cols[2]

            pos_key = ':'.join([chrom, end])
            poi[pos_key] = 1

    return(poi)

def _parse_cmd_line():
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter,
        description=__doc__,
    )

    parser.add_argument(
        '--mpileup',
        '-m',
        help='path to mpileup of variants file.  If unspecfied, STDIN will be used.'
    )

    parser.add_argument(
        '--variant-bed',
        '-b',
        help='path to bed of variants.  Readcounts will be filtered for just these regions.'
    )
    parser.add_argument(
        '--sample',
        '-s',
        required=True,
        help='Sample ID, e.g., GTEX-ZZPU-0326-SM-5N9BJ.  Sample ID is output as a column in the STDOUT'
    )

    parser.add_argument(
        '--header',
        '-hd',
        action='store_true',
        help='if specified, a header will be printed.  Useful when not concatenating readcount files and when you want to know what those darn columns are.'
    )

    parser.add_argument(
        '--read-length',
        '-r',
        type=int,
        default=76,
        help='Length of sequencing reads.  Value is used to determine if alt alleles are near the ends of sequencing reads.'
    )

    args = parser.parse_args()

    return args
def _get_max_internal_base_pos(base_pos_counts, read_length):
    """Distance is relative to the closer end of the read"""
    max_internal_base_pos = {}
    read_midpoint = (read_length - 1)/ 2 + 1
    for base, base_positions in base_pos_counts.items():
        if len(base_positions):
            m = -1 # m will store the maxmimum internal base position.  Initialize it to a value that will never be a distance so can make sure it was set.
            for base_pos in base_positions:
                pos_from_closer_end = min(base_pos - 1, read_length - base_pos) # -1 so when base_pos = 1 (1-based), the distance = 0; similarily, when base_pos = read_length, disatnce = 0.  Note.  The distances are a little screwy if there is soft clipping, but we've already removed those heavily softclipped reads so it shouldn't affect the distance too much.
                if pos_from_closer_end < 0:
                    print('ERROR: read length must not be correct.  There are read positions (%d) greater than the specified read length (%d)' % (base_pos, read_length), file=sys.stderr)
                    sys.exit(1)
                if pos_from_closer_end >= m:
                    m = pos_from_closer_end
    
                    if abs(m - read_midpoint) <= 0.5: # if m = the midpoint (when read length is odd) or if m is within 1/2 of midpoint (when the read length is even), there can be no closer reads to the middle
                        break
            if m == -1:
                print('ERROR: could not find max internal base position for base pos counts \'%s\'' % (base_pos_counts), file=sys.stderr)
                sys.exit(1)
        else:
            m = 'NA' # If there is no coverage for this base, statistic is undefined
        max_internal_base_pos[base] = m

    return(max_internal_base_pos)

def _get_max_internal_base_pos_old(base_pos_counts, read_midpoint):
    """Distance is releative to the midpoint of the read"""

    max_internal_base_pos = {}
    for base, base_positions in base_pos_counts.items():
        if len(base_positions):
            m = float('inf') # Maxmimum internal base position
            for base_pos in base_positions:
                pos_from_middle = abs(base_pos - read_midpoint)
                if pos_from_middle <= m:
                    m = pos_from_middle
    
                    if abs(m - read_midpoint) <= 0.5: # if m = the midpoint (when read length is odd) or if m is within 1/2 of midpoint (when the read length is even), there can be no closer reads to the middle
                        break
    
            if m == float('inf'):
                print('ERROR: could not find max internal base position for base pos counts \'%s\'' % (base_pos_counts), file=sys.stderr)
                sys.exit(1)
        else:
            m = 'NA' # If there is no coverage for this base, statistic is undefined
        max_internal_base_pos[base] = m

    return(max_internal_base_pos)

def _get_alt_allele(unstranded_read_counts, ref_allele):
    # Alt allele is defined as the allele with the highest non-reference coverage. 
    # If there are ties, the first allele, when sorted alphabetically, is used.

    base2index = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    index2base = {i:b for b,i in base2index.items()}

    ref_index = base2index[ref_allele]

    unstranded_read_counts_arr = [unstranded_read_counts[i] for i in 'ACGT']

    sort_order = np.argsort(-1*np.asarray(unstranded_read_counts_arr)) # Sort from largest to smallest
    if sort_order[0] == ref_index: # If ref has highest coverage, alt has 2nd highest coverage
        alt_index = sort_order[1]
    else: # If ref does not have highest coverage, alt has highest coverage
        alt_index = sort_order[0]

    alt_allele = index2base[alt_index]
    
    return(alt_allele)

def _stranded2unstranded(x):

    unstranded = {} # key = base (unstranded); value = sum of base and BASE.  conveniently works well when values are numeric or lists
    for i in 'ACGTN':
        unstranded[i] =  x[i] + x[i.lower()]

    return(unstranded)

def _get_counts(base_str, base_pos_str, ref_allele):

    base_pos = list(map(int, base_pos_str.split(',')))
    i_base = 0
    j_pos = 0

    read_counts, base_pos_counts = _initialize_counts(ref_allele)

    while i_base < len(base_str):

        base = base_str[i_base]
        if base in read_counts:
            read_counts[base] += 1
            base_pos_counts[base].append(base_pos[j_pos])
            i_base += 1
            j_pos += 1
        elif base == "^": # Start of read
            i_base += 2 # Skip this and the next (alignment quality)
        elif base == '$': # End of read
            i_base += 1
        elif base == '>' or base == '<': # reference skip
            i_base += 1
            j_pos += 1
        elif base == '-' or base == '+': # Indel
            # -[0-9]+[ACGTNacgtn], example: -2AT
            m = re.match(r'[+-]([0-9]+)', base_str[i_base:])
            if m:
                indel_length = m.group(1)
                i_base += 1 + len(indel_length) + int(indel_length) # skip the strand (+1), indel length number and the bases in the indel
            else:
                print('ERROR: could not parse indel length starting at base string \'%s\'' % base_str[i_base:], file=sys.stderr)
                sys.exit(1)
        elif base == '*': # There was a deletion upstream
            i_base += 1
            j_pos += 1
        else:
            print('ERROR: got unknown character in base string \'%s\'', base, file=sys.stderr)
            sys.exit(1)

    # Check that we're at the end of the base pos
    if j_pos != len(base_pos):
        print('ERROR: did not parse the read base positions correctly', file=sys.stderr)
        sys.exit(1)

    
    read_counts = _finialize_counts(read_counts, ref_allele)
    base_pos_counts = _finialize_counts(base_pos_counts, ref_allele)

    return(read_counts, base_pos_counts)

def _print_line_no_cov(chrom, pos_1based, ref_allele, sample):

    alt_allele = 'NA' 
    read_counts = {i:0 for i in "ACGTNacgtn"}

    line = ([chrom, pos_1based, ref_allele, alt_allele, 0, 0, 0, 0, 0, 'NA'] + 
           [read_counts[i] for i in 'AaCcGgTtNn'] +
           [sample])

    print('\t'.join(list(map(str, line))))
    return None  
def _print_line(chrom, pos_1based, ref_allele, alt_allele, read_counts, unstranded_read_counts, max_internal_base_pos, sample):

    lc_alt_allele = alt_allele.lower()
    if (read_counts[alt_allele] > 0 and read_counts[lc_alt_allele] > 0):
        alt_on_both_strands = 1
    else:
        alt_on_both_strands = 0    

    line = ([chrom, 
             pos_1based, 
             ref_allele, 
             alt_allele, 
             unstranded_read_counts[ref_allele],
             unstranded_read_counts[alt_allele], 
             read_counts[alt_allele],
             read_counts[lc_alt_allele],
             alt_on_both_strands,
             max_internal_base_pos[alt_allele]] + 
           [read_counts[i] for i in 'AaCcGgTtNn'] +
           [sample])

    print('\t'.join(list(map(str, line))))
    return None  

def _print_header():

    print('\t'.join(['chr', 'pos', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'alt_pos_count', 'alt_neg_count', 'alt_on_both_strands', 'max_alt_pos'] + ['%s_count\t%s_count' % (i, i.lower()) for i in 'ACGTN'] + ['sample']))

    return None

def _finialize_counts(x, ref_allele):
    # Convert . and , back to uc ref and lc ref
    lc_ref_allele = ref_allele.lower()

    x[ref_allele] = x['.']
    x[lc_ref_allele] = x[',']

    del x['.']
    del x[',']

    return(x)

def _initialize_counts(ref_allele):
    read_counts = {i:0 for i in "ACGTNacgtn.,"} # key = base; value = count of base
    base_pos_counts = {i:[] for i in "ACGTNacgtn.,"} # key = base; value = list of base positions
    # Change ref to . and ,
    lc_ref_allele = ref_allele.lower()

    del read_counts[ref_allele]
    del read_counts[lc_ref_allele]

    del base_pos_counts[ref_allele]
    del base_pos_counts[lc_ref_allele]

    return(read_counts, base_pos_counts)


if __name__ == '__main__':
    main()
