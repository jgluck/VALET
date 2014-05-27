#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import math
import os
import time
import sys

def setup_options():
    parser = OptionParser()
    parser.add_option("-a", "--abundance-file", dest="abundance_file", help="Assembler provided abundance file", metavar="FILE")
    parser.add_option("-m", "--mpileup-file", dest="mp_file_loc", help="mpileup output file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_location", help="output location for depth of coverage tool", default="data/output/coverage/bad_bases.csv", metavar="FILE")
    parser.add_option("-w", "--window-size", dest="window_size", help="window size to sweep over base pairs", type="int", default = 1)
 
    (options,args) = parser.parse_args()

    should_quit = False
    if options.abundance_file == None:
        parser.error("You failed to provide an abundance file")
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
    if not os.path.exists(d):
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
    read_abundances(abundance_file, abundance_dict)
    
    
    with open(bad_cvg_file, 'w') as fp:
        for grouping in read_mpileup(mpile_file, window_size):
            check_grouping(grouping, abundance_dict, fp, 3)


def in_range(num, low, high):
    return num >= low and num <= high

def check_grouping(grouping, a_dict, fp, n_devs):
    cvg = 0.0
    window_start = int(grouping[0][1])
    window_end = int(grouping[-1][1])
    win_size = 0.0
    for group in grouping:
        cvg += float(group[3])
        win_size+=1.0
    
    avg_contig_cvg = float(a_dict[grouping[0][0]])
    a_cvg = cvg / win_size

    std_dev = avg_contig_cvg / 10.0
    low = avg_contig_cvg - (std_dev * n_devs)
    high = avg_contig_cvg + (std_dev * n_devs)
    if not in_range(a_cvg, low , high) :
        fp.write("%d, %d, %d, %f, %f, %f\n" %(int(grouping[0][0]), window_start, \
            window_end, a_cvg, low, high))

def read_abundances(fp, a_dict):
    fp = open(fp,'r')
    for line in fp:
        (contig, avg_cov) = line.split()
        a_dict[contig] = float(avg_cov)
    fp.close()

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
    

if __name__ == '__main__':
    main()
