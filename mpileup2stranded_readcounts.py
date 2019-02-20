#!/usr/bin/python3
"""Convert mpileup format to stranded readcounts format.  Stranded readcounts are printed to STDOUT.

OUTPUT COLUMNS:
If --output-format == 'txt': the first columns are:
   1  chr = chromosome
   2  pos = position (1-based)

If --output-format == 'bed': the first columns are:
   1  chr = chromosome
   2  start = start (0-based, INclusive)
   3  end = end (1-based, EXclusive)

The remaining columns are (indexes are relative to a txt output file):
   3  ref_allele = reference allele
   4  alt_allele = alternative allele (non-reference allele with highest coverage; if ties, the first allele when sorted alphabetically is used)
   5  ref_count = total ref coverage
   6  alt_count = total alt coverage
   7  alt_pos_count = alt coverage on + strand
   8  alt_neg_count = alt coverage on - strand
   9  alt_on_both_strands = 1 if alt was observed on + and - strand; 0 otherwise
  10  max_alt_pos = If position has > 0 alt, most internal position of the alt allele on the reads, e.g., 2 means the alt allele was observed 2bp from the end.  No other read containing the alt read had an alt base closer to the center of the read.  Else, = NA.  Note: distance is relative to the *closer end* of the read.
  11  avg_alt_pos = If position has > 0 alt, average position of the alt allele on the reads.  Else, = NA.  Note: distance is relative to the *start* of the read.
  11  std_alt_pos = If position has > 0 alt, standard deviation of the alt allele on the reads.  Else, = NA.  Note: distance is relative to the *start* of the read.
  12  A_count = A coverage on + strand
  13  a_count = A coverage on - strand
  14  C_count = C coverage on + strand
  15  c_count = C coverage on - strand
  16  G_count = G coverage on + strand
  17  g_count = G coverage on - strand
  18  T_count = T coverage on + strand
  19  t_count = T coverage on - strand
  20  N_count = N coverage on + strand
  21  n_count = N coverage on - strand
  22  sample = Sample ID (useful when concatenating multiple readcount files)

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

    if args.region_bed is None:
        poi = {}
        _process_pos_ref = _process_pos_dummy
    else:
        poi = _load_poi(args.region_bed)
        _process_pos_ref = _process_pos_poi

    if args.output_format == 'txt':
        _get_position_format_ref = _get_txt_position_format
    elif args.output_format == 'bed':
        _get_position_format_ref = _get_bed_position_format
    else:
        print('ERROR: unknown output format \'%s\'.' % (args.output_format))
        sys.exit(1)

    if args.header:
        _print_header(args.output_format)

    if args.mpileup is None:
        mpileup_fh = sys.stdin
    else:
        mpileup_fh = open(args.mpileup, 'rU')

    num_ambiguous_ref = 0
    num_poi = 0
    for line in mpileup_fh:
        line = line.rstrip('\n')
        cols = line.split('\t')
        chrom = cols[0]
        pos_1based = cols[1]

        pos_key = ':'.join([chrom, pos_1based])

        # Only process if this is a position of interest
        if _process_pos_ref(pos_key, poi):

            ref_allele = cols[2].upper()

            # Only process if the reference is unambiguous
            if ref_allele in 'ACGT':
                depth = int(cols[3])
        
                if depth > 0:
                    base_str = cols[4]
                    #base_quals = cols[5]
                    base_pos_str = cols[6]
        
                    read_counts, base_pos_counts, total_count = _get_counts(base_str, base_pos_str, ref_allele)
    
                    if total_count > 0:
                        unstranded_read_counts = _stranded2unstranded(read_counts)
                        unstranded_base_pos_counts = _stranded2unstranded(base_pos_counts)
                        alt_allele = _get_alt_allele(unstranded_read_counts, ref_allele)
            
                        # Calculate the max internal, average position, and standard deviation of each base
                        max_internal_base_pos, avg_base_pos, std_base_pos = _get_pos_statistics(unstranded_base_pos_counts, args.read_length, line)
            
                        _print_line(chrom, pos_1based, ref_allele, alt_allele, read_counts, unstranded_read_counts, max_internal_base_pos, avg_base_pos, std_base_pos, args.sample, _get_position_format_ref)
                    elif args.print_sites_w_no_coverage: # This can happen, e.g., only spliced reads overlap the position
                        _print_line_no_cov(chrom, pos_1based, ref_allele, args.sample, _get_position_format_ref)
                elif args.print_sites_w_no_coverage: # There is no coverage here
        
                    _print_line_no_cov(chrom, pos_1based, ref_allele, args.sample, _get_position_format_ref)
            else:
                num_ambiguous_ref += 1
            num_poi += 1

    mpileup_fh.close()
    print("INFO: number of POI with ambiguous reference alleles=%.1E%% (%d/%d)" % (num_ambiguous_ref/num_poi*100 , num_ambiguous_ref, num_poi), file=sys.stderr)

def _load_poi(bed_fn):

    poi = {} # key = chrom:1-based position; value = 1 (dummy).  (Use 1-based since mpileup is 1-based)
    with open(bed_fn, 'rU') as bed_fh:
        for line in bed_fh:
            line = line.rstrip('\n')
            cols = line.split('\t')

            # Add each position in the region to the poi (this was an easier hack than merging and comparing an mpileup position to regions)
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])

            for i in range(start,end):
                pos_key = ':'.join([chrom, str(i + 1)]) # Use 1-based.  See note above.
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
        '--region-bed',
        '-b',
        help='path to bed of regions/variants.  Readcounts will be filtered for just these regions.',
        required=False
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
        '--print-sites-w-no-coverage',
        '-p',
        action='store_true',
        help='if specified, print sites with no coverage.  These are sites that only have spliced read coverage and/or are sites in the regions of interest that had 0X coverage.'
    )

    parser.add_argument(
        '--read-length',
        '-r',
        type=int,
        default=76,
        help='Length of sequencing reads.  Value is used to determine if alt alleles are near the ends of sequencing reads.'
    )

    parser.add_argument(
        '--output-format',
        '-o',
        default='txt',
        choices=['bed','txt'],
        help='output format.  Bed outputs a stranded readcount in bed format with chrom/start/end.  Txt outputs a stranded readcount with just chrom/1-based pos.'
    )

    args = parser.parse_args()

    return args

def _get_txt_position_format(chrom, pos_1based):

    return([chrom, pos_1based])

def _get_bed_position_format(chrom, pos_1based):

    return([chrom, int(pos_1based) - 1, pos_1based])

def _get_pos_statistics(base_pos_counts, read_length, line):
    """For each base, calculate the 
        - maximum internal position (note: distance is relative to the closer end of the read)
        - average position (note: distance is relative to the start of the read)
        - standard deviation position (note: distance is relative to the start of the read)"""
    max_internal_base_pos = {}
    min_internal_base_pos = {} # Used for checking if submitted read length may not be correct
    avg_base_pos = {}
    std_base_pos = {}
    read_midpoint = (read_length - 1)/ 2 + 1
    for base, base_positions in base_pos_counts.items():
        if len(base_positions):

            pos_from_closer_end_list = [min(base_pos - 1, read_length - base_pos) for base_pos in base_positions] # -1 so when base_pos = 1 (1-based), the distance = 0; similarily, when base_pos = read_length, disatnce = 0.  Note.  The distances are a little screwy if there is soft clipping, but we've already removed those heavily softclipped reads so it shouldn't affect the distance too much.
            
            # If there are negative positions, the read length must not have been correct
            max_internal_base_pos[base] = np.max(pos_from_closer_end_list)
            min_internal_base_pos[base] = np.min(pos_from_closer_end_list)
            if min_internal_base_pos[base] < 0:
                print('ERROR: read length must not be correct.  There are read positions (%d) greater than the specified read length (%d).  Line was:\n    %s' % (base_pos, read_length, line), file=sys.stderr)
                sys.exit(1)

            avg_base_pos[base] = np.average(base_positions)
            std_base_pos[base] = np.std(base_positions)
        else:
            # If there is no coverage for this base, statistic is undefined
            max_internal_base_pos[base] = 'NA'
            avg_base_pos[base] = 'NA'
            std_base_pos[base] = 'NA'

    return(max_internal_base_pos, avg_base_pos, std_base_pos)


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
    total_count = sum(read_counts.values())

    return(read_counts, base_pos_counts, total_count)

def _print_line_no_cov(chrom, pos_1based, ref_allele, sample, _get_position_format_ref):

    alt_allele = 'NA' 
    read_counts = {i:0 for i in "ACGTNacgtn"}

    position_cols = _get_position_format_ref(chrom, pos_1based)
    line = (position_cols + 
           [ref_allele, alt_allele, 0, 0, 0, 0, 0, 'NA'] + 
           [read_counts[i] for i in 'AaCcGgTtNn'] +
           [sample])

    print('\t'.join(list(map(str, line))))
    return None  
def _print_line(chrom, pos_1based, ref_allele, alt_allele, read_counts, unstranded_read_counts, max_internal_base_pos, avg_base_pos, std_base_pos, sample, _get_position_format_ref):

    lc_alt_allele = alt_allele.lower()
    if (read_counts[alt_allele] > 0 and read_counts[lc_alt_allele] > 0):
        alt_on_both_strands = 1
    else:
        alt_on_both_strands = 0    

    position_cols = _get_position_format_ref(chrom, pos_1based)

    line = (position_cols +
           [ref_allele, 
            alt_allele, 
            unstranded_read_counts[ref_allele],
            unstranded_read_counts[alt_allele], 
            read_counts[alt_allele],
            read_counts[lc_alt_allele],
            alt_on_both_strands,
            max_internal_base_pos[alt_allele],
            avg_base_pos[alt_allele],
            std_base_pos[alt_allele]] + 
           [read_counts[i] for i in 'AaCcGgTtNn'] +
           [sample])

    print('\t'.join(list(map(str, line))))
    return None  

def _print_header(output_format):

    if output_format == 'txt':
        print('\t'.join(['chr', 'pos', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'alt_pos_count', 'alt_neg_count', 'alt_on_both_strands', 'max_internal_alt_pos', 'avg_alt_pos', 'std_alt_pos'] + ['%s_count\t%s_count' % (i, i.lower()) for i in 'ACGTN'] + ['sample']))
    elif output_format == 'bed':
        print('\t'.join(['chr', 'start', 'end', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'alt_pos_count', 'alt_neg_count', 'alt_on_both_strands', 'max_internal_alt_pos', 'avg_alt_pos', 'std_alt_pos'] + ['%s_count\t%s_count' % (i, i.lower()) for i in 'ACGTN'] + ['sample']))
    else:
        print('ERROR: unknown output format \'%s\'.' % (args.output_format))
        sys.exit(1)

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

def _process_pos_dummy(pos_key, poi):
    """Dummy function to always return True.  Needed if no POI were submitted"""

    process = True

    return(process)

def _process_pos_poi(pos_key, poi):

    process = False
    if pos_key in poi:
        process = True

    return(process)

if __name__ == '__main__':
    main()
