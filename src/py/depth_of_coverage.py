#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import math
import os
import time
import sys

PROG_NAME = "DEPTH_COV"

def setup_options():
    parser = OptionParser()
    parser.add_option("-a", "--abundance-file", dest="abundance_file", help="Assembler provided abundance file", metavar="FILE")
    parser.add_option("-m", "--mpileup-file", dest="mp_file_loc", help="mpileup output file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_location", help="output location for depth of coverage tool", default="data/output/coverage/bad_bases.csv", metavar="FILE")
    parser.add_option("-w", "--window-size", dest="window_size", help="window size to sweep over base pairs", type="int", default = 1)
    parser.add_option("-g", "--gff", dest="gff_format", default=False, action='store_true')
    parser.add_option("-e", "--empirical", dest="use_empirical", default=False, action='store_true')
    parser.add_option("-i", "--ignore", dest="ignore_ends", help="Ignore the first/last i bps of the read", type="int", default=0)

    (options,args) = parser.parse_args()

    should_quit = False
    #if options.abundance_file == None:
    #    parser.error("You failed to provide an abundance file")
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
    bad_cvg_dir = os.path.dirname(options.output_location)
    bad_cvg_file = options.output_location
    ensure_dir(bad_cvg_dir)

    window_size = options.window_size

    abundance_dict = {}
    
    if options.abundance_file:
        read_abundances(abundance_file, abundance_dict)
    else:
        # Calculate the coverage and window from the data itself.
        calculate_coverages(mpile_file, abundance_dict)

    
    with open(bad_cvg_file, 'w') as fp:
        for grouping in read_mpileup(mpile_file, window_size):
            check_grouping(grouping, abundance_dict, fp, 3, options.gff_format, options.use_empirical)


def in_range(num, low, high):
    return num >= low and num <= high

def check_grouping(grouping, a_dict, fp, n_devs, gff_format, use_empirical):
    #print(grouping)
    cvg = 0.0
    window_start = int(grouping[0][1])
    window_end = int(grouping[-1][1])
    win_size = 0.0
    for group in grouping:
        cvg += float(group[3])
        win_size+=1.0

    avg_contig_cvg = 0
    if use_empirical:
        avg_contig_cvg = float(a_dict[grouping[0][0]][0])
    else:
        avg_contig_cvg = float(a_dict[grouping[0][0]])

    a_cvg = cvg / win_size

    std_dev = avg_contig_cvg / 10.0
    low = avg_contig_cvg - (std_dev * n_devs)
    high = avg_contig_cvg + (std_dev * n_devs)

    if use_empirical:
        # Grab the lower and upper hinges.
        low = a_dict[grouping[0][0]][1]
        high = a_dict[grouping[0][0]][2]

    if not in_range(a_cvg, low , high) :
        if gff_format:
            cov_type = "Low_coverage"
            color = "#7800ef"
            if a_cvg > high:
                cov_type = "High_coverage"
                color = "#0077ee"
            fp.write("%s\t%s\t%s\t%d\t%d\t%f\t.\t.\tlow=%f;high=%f;color=%s\n" %(str(grouping[0][0]), PROG_NAME, \
                cov_type, window_start, window_end, a_cvg, low, high, color))
            
        else:
            fp.write("%s,%d,%d,%f,%f,%f\n" %(str(grouping[0][0]), window_start, \
                window_end, a_cvg, low, high))

        

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

def read_mpileup(fp, ws=1):
    fp = open(fp, 'r')
    last_contig = None
    save_line = None
    while True:
        line_bundle = []
        ret_bundle = []
        last_contig = None

        start_range = 0
        if save_line:
            start_range = 1
            line_bundle.append(save_line)
            save_line = None

        for i in range(start_range,ws):
            line_bundle.append(fp.readline())
        for line in line_bundle:
            if line:
                split_line = line.split()
                if last_contig and split_line[0] != last_contig:
                    save_line = line
                else:
                    last_contig = split_line[0]
                ret_bundle.append(split_line)
        if len(ret_bundle) == 0:
            break
        yield ret_bundle
    

def factorial(n):
    """ Return factorial using Stirling's approximation. """

    return math.sqrt(2 * math.pi * n) * (n / math.e) ** n


if __name__ == '__main__':
    main()
