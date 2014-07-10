#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import collections
import copy
import math
import os
import sys
import time


PROG_NAME = "DEPTH_COV"


def setup_options():
    parser = OptionParser()
    parser.add_option("-a", "--abundance-file", dest="abundance_file", help="Assembler provided abundance file", metavar="FILE")
    parser.add_option("-m", "--mpileup-file", dest="mp_file_loc", help="mpileup output file", metavar="FILE")
    parser.add_option("-p", "--pileup-file", dest="pileup_filename", help="pileup file name", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_location", help="output location for depth of coverage tool", default=None, metavar="FILE")
    parser.add_option("-w", "--window-size", dest="window_size", help="window size to sweep over base pairs", type="int", default = 1)
    parser.add_option("-g", "--gff", dest="gff_format", default=False, action='store_true')
    parser.add_option("-e", "--empirical", dest="use_empirical", default=False, action='store_true')
    parser.add_option("-i", "--ignore", dest="ignore_ends", help="Ignore the first/last i bps of the read", type="int", default=0)

    (options,args) = parser.parse_args()

    should_quit = False
    if options.mp_file_loc == None:
        parser.error("You failed to provide the mpileup file")
    if should_quit:
        parser.help()
        exit(-1)

    return (options,args)


def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)


def debug(*objs):
    print("DEBUG: ", *objs, file=sys.stderr)


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d) and d:
        print("true")
        os.makedirs(d)


def main():
    (options, args) = setup_options()
    abundance_file = options.abundance_file
    mpile_file = options.mp_file_loc
    
    if options.output_location:
        ensure_dir(os.path.dirname(options.output_location))
    bad_cvg_file = open(options.output_location,'w') if options.output_location else sys.stdout
    
    window_size = options.window_size

    abundance_dict = {}
    
    if options.abundance_file:
        read_abundances(abundance_file, abundance_dict)
    else:
        # Calculate the coverage and window from the data itself.
        calculate_coverages(mpile_file, abundance_dict)

    coverage_window = collections.deque(maxlen = window_size)
    flagged_regions = []

    prev_contig = None
    lower_hinge = 0
    upper_hinge = 0
    end_pos = 0
    region_index = -1

    for line in open(mpile_file, 'r'):
        fields = line.split()

        # contig00001     1       A       1       ^~,     I
        if prev_contig is None:
            prev_contig = fields[0]
            lower_hinge = abundance_dict[prev_contig][1]
            upper_hinge = abundance_dict[prev_contig][2]
            end_pos = 0

        if prev_contig != fields[0]:
            # Output previous results and clear deque.
            coverage_window.clear()
            prev_contig = fields[0]
            lower_hinge = abundance_dict[prev_contig][1]
            upper_hinge = abundance_dict[prev_contig][2]
            end_pos = 0

        # Append the bp coverage to the window.
        coverage_window.append(int(fields[3]))
        end_pos += 1

        # If the coverage window is full, check and compare with median.
        if len(coverage_window) >= window_size:
            # TODO: Deepcopy is inefficient.
            copy_window = sorted(copy.deepcopy(coverage_window))
            median = copy_window[len(copy_window) / 2]
            if not len(copy_window) % 2:
                median = (copy_window[len(copy_window) / 2] + copy_window[len(copy_window) / 2 - 1]) / 2.0

            if not in_range(median, lower_hinge, upper_hinge):
                cov_type = "Low_coverage"
                color = "#7800ef"
                if median > upper_hinge:
                    cov_type = "High_coverage"
                    color = "#0077ee"

                # Extend previous window?
                if len(flagged_regions) > 0 and \
                        flagged_regions[region_index][0] == prev_contig and \
                        flagged_regions[region_index][2] == cov_type and \
                        int(flagged_regions[region_index][4]) >= end_pos - len(coverage_window) + 1 and \
                        int(flagged_regions[region_index][4]) <= end_pos:
                    flagged_regions[region_index][4] = end_pos

                else:
                    flagged_regions.append(\
                        [prev_contig, PROG_NAME, cov_type, end_pos - len(coverage_window) + 1, end_pos, median, lower_hinge, upper_hinge, color])
                    region_index += 1
 
    for region in flagged_regions:
        bad_cvg_file.write("%s\t%s\t%s\t%d\t%d\t%f\t.\t.\tlow=%f;high=%f;color=%s\n" % (region[0], region[1], \
                    region[2], region[3], region[4], region[5], region[6], region[7], region[8]))


def in_range(num, low, high):
    return num >= low and num <= high


def read_abundances(fp, a_dict):
    fp = open(fp,'r')
    for line in fp:
        (contig, avg_cov) = line.split()
        a_dict[contig] = float(avg_cov)
    fp.close()


def tukey_summary(array):
    """ Given an array of integers, return median, and lower/upper hinge tuple. """

    array.sort()

    median = array[len(array) / 2]
    if not len(array) % 2:
        median = (array[len(array) / 2] + array[len(array) / 2 - 1]) / 2.0

    first_quantile = array[len(array) / 4]
    third_quantile = array[3 * len(array) / 4]

    iqr = third_quantile - first_quantile

    lower_hinge = first_quantile - 1.5 * iqr
    upper_hinge = third_quantile + 1.5 * iqr

    return (median, lower_hinge, upper_hinge)


def calculate_coverages(mpileup_file, abundance_dict):
    """ For each contig, calculate coverage and lower/upper hinges. """

    prev_contig = None
    length = 0
    curr_coverages = []

    for record in open(mpileup_file, 'r'):
        fields = record.strip().split()

        if prev_contig != fields[0]:
            if prev_contig:
                # Calculate median coverage, and lower/upper hinges.
                abundance_dict[prev_contig] = tukey_summary(curr_coverages)
                
            prev_contig = fields[0]
            length = 0
            curr_coverages = []

        curr_coverages.append(int(fields[3]))
        length += 1

    if prev_contig:
        abundance_dict[prev_contig] = tukey_summary(curr_coverages)


def get_average_coverage(a_dict):
    l = a_dict.values()
    mean =  reduce(lambda x, y: x + y, l) / len(l)
    variances = map(lambda x: (x - avg)**2, l)
    var = reduce(lambda x, y: x + y, variances)/ len(variances)
    std = math.sqrt(var)
    return (mean, var, std)


if __name__ == '__main__':
    main()
