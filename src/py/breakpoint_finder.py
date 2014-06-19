#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
from subprocess import call
import os
import sys
import math

class BreakpointFinder:

    def __init__(self):
        self.options = None
        self.assembly_file = None
        self.singleton_halves = []
        self.read_lengths = {}
        self.getOptions()
        self.set_locations()
    
    def set_locations(self):
        self.bowtie_dir = "/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/"
        self.bowtie_loc = self.bowtie_dir + "bowtie2"
        self.bowtie_build_loc = self.bowtie_dir + "bowtie2-build"
        self.breakpoint_dir = self.options.output_dir

        self.bowtie_index =  self.breakpoint_dir + "bowtie-index/"
        self.index_prefix = self.bowtie_index+ "breakpoint"
        ensure_dir(self.bowtie_index)
        
        self.sam_output_dir = self.breakpoint_dir + "sam/"
        self.sam_output_location = self.sam_output_dir + "breakpoints.sam"
        ensure_dir(self.sam_output_dir)

        self.singleton_dir = self.breakpoint_dir + "singleton/"
        self.singleton_loc = self.singleton_dir + "singletons.csv"

        ensure_dir(self.singleton_dir)

        self.conc_dir = self.breakpoint_dir + "concordance/"
        self.conc_loc = self.conc_dir + "concordants.csv"
        ensure_dir(self.conc_dir)

        self.breakpoint_file = self.breakpoint_dir + "breakpoints.csv"
        self.sorted_breakpoint_file = self.breakpoint_dir + "sorted_breakpoints.csv"
        self.binned_breakpoint_file = self.breakpoint_dir + "binned_breakpoints.csv"
        self.collapsed_breakpoint_file = self.breakpoint_dir\
                + "collapsed_breakpoints.csv"
        self.bins_of_interest_file = self.breakpoint_dir + "interesting_bins.csv"

    def run_bowtie_index(self):
        FNULL = open('/dev/null', 'w')
        call_arr = [self.bowtie_build_loc, self.assembly_file, self.index_prefix]
        #out_cmd(call_arr)
        call(call_arr, stdout = FNULL, stderr = FNULL)

    def run_bowtie_2(self):
        for file_name in os.listdir(self.reads_dir):
            if "reads" in file_name:
                call_arr = [self.bowtie_loc, '-x', self.index_prefix, '-U',\
                         self.reads_dir + file_name,\
                         '-S', self.sam_output_dir + file_name + '.sam',\
                         '--un', self.singleton_dir + file_name + '.singletons',\
                         '--al', self.conc_dir + file_name + '.reads',
                         '-q', '-I 50', '-X 800', '-p 10']
                #out_cmd(call_arr)
                call(call_arr)
            else:
                warning("Skipping potential read file: " + file_name)

    def sort_breakpoints(self):
        call_arr = ["sort", "-t\t", '-k 1n,1', '-k 2n,2', '-k 5n,5', self.breakpoint_file]
        out_cmd(call_arr)
        out_file = open(self.sorted_breakpoint_file, 'w')
        warning("This call outputs to file: ", self.sorted_breakpoint_file)
        call(call_arr, stdout=out_file)

    def bin_breakpoints(self):
        with open(self.sorted_breakpoint_file,'r') as breakpoints,\
                open(self.binned_breakpoint_file,'w') as out_file:
            for contig_bundle in self.read_contig(breakpoints):
                self.w_s = 1
                self.w_e = self.w_s + self.bin_size
                for match in contig_bundle:
                    match_split = match.split('\t')
                    match_start = int(match_split[1])
                    match_end = int(match_split[4]) + match_start
                    if match_start == 0:
                        continue
                    if match_start >= self.w_s and match_start <= self.w_e:
                        out_file.write(match.strip()+('\t%s\n'%(self.w_s)))
                    else:
                        while match_start < self.w_s or match_start > self.w_e:
                            self.w_s = self.w_e
                            self.w_e = self.w_s + self.bin_size
                        out_file.write(match.strip()+('\t%s\n'%(self.w_s)))

    def collapse_bins(self):
        call_str = "awk '{printf \"%%s\\t%%s\\n\",$1,$6}' %s | sort | uniq -c | tee %s" %(self.binned_breakpoint_file, self.collapsed_breakpoint_file)
        out_cmd(call_str)
        call(call_str,shell=True)

    def trim_bins(self):
        avg_bin_size = 0
        num_bins = 0
        summed_var = 0
        std_dev = 0
        with open(self.collapsed_breakpoint_file,'r') as pass_1:
            for line in pass_1:
                avg_bin_size += int(line.split()[0])
                num_bins += 1
            avg_bin_size = float(avg_bin_size)/num_bins
        with open(self.collapsed_breakpoint_file,'r') as pass_2:
            for line in pass_2:
                summed_var += (float(line.split()[0]) - avg_bin_size)**2
            summed_var = summed_var / num_bins
            std_dev = math.sqrt(summed_var)
        with open(self.collapsed_breakpoint_file,'r') as pass_3,\
                open(self.bins_of_interest_file,'w') as out_file:
                    cutoff = 3 * std_dev
                    for line in pass_3:
                        split_line = line.split()
                        if (float(split_line[0]) - avg_bin_size) > cutoff:
                            out_file.write("%s\tBreakpoint_finder\tExcessive_alignment\
                                    \t%d\t%d\t%d\t.\t.\tavg_bin_size=%f;std_dev=%f\n"\
                                    %(split_line[1],split_line[2],\
                                    split_line[2]+self.bin_size,split_line[0],\
                                    avg_bin_size, std_dev))



    def read_contig(self, fp):
        run_flag = True
        contig_bundle = []
        line_buffer = ""
        while run_flag:
            current_contig = ""
            if line_buffer != "":
                current_contig = line_buffer.split('\t')[0]
                contig_bundle.append(line_buffer)
                line_buffer = ""
            while True:
                line = fp.readline()
                if line == '':
                    run_flag = False
                    break
                elif line.split('\t')[0] != current_contig:
                    line_buffer = line
                    break
                else:
                    contig_bundle.append(line)
            ret_bundle = contig_bundle
            contig_bundle = []
            yield ret_bundle
 
    def read_in_lengths(self):
        log_file = open(self.breakpoint_dir + "log.log", 'w')
        for file_name in os.listdir(self.conc_dir):
            print("Checking name %s\n" %(file_name))
            if '.reads' in file_name:
                print("Reading file: %s\n" %(file_name))
                with open(self.conc_dir + file_name, 'r') as reads_file:
                    for (read,length) in self.read_read(reads_file):
                        #print ("Read: %s Length: %s" %(read,length))
                        self.read_lengths[read] = length
    
    def read_read(self, fp):
        run_flag = True
        while run_flag:
            line_bundle = []
            for i in range(4):
                line_bundle.append(fp.readline())
            if line_bundle[0] == '':
                run_flag = False
            else:
                read = line_bundle[0][1:]
                length = len(line_bundle[1].strip())
                yield (read,length)

    '''
    Outputs breakpoints to file self.breakpoint_file
    Position Sequence Name Flag Length
    Tab sep
    '''
    def detect_breakpoints(self):
        with open(self.breakpoint_file, 'w') as out_file:
            for file_name in os.listdir(self.sam_output_dir):
                if "sam" in file_name:
                    with open(self.sam_output_dir+file_name, 'r') as in_file:
                        for line in in_file:
                            if line[0] != '@':
                                line_components = line.split('\t')
                                if line_components[3] != '0':
                                    out_file.write("%s\t%s\t%s\t%s\t%d\n" % (line_components[2],\
                                            line_components[3],\
                                            line_components[0], line_components[1],\
                                            len(line_components[9])))

    def go(self):
        self.run_bowtie_index()
        self.run_bowtie_2()
        self.detect_breakpoints()
        self.sort_breakpoints()
        self.bin_breakpoints()
        self.collapse_bins()
        self.trim_bins()
    
    def getOptions(self):
        parser = OptionParser()
        parser.add_option("-a", "--assembly-file", dest="assembly_file",\
                help="Assembly File to search", metavar="FILE")
        parser.add_option("-r", "--reads-dir", dest="reads", \
                help="Directory of Reads", metavar="PATH")
        parser.add_option("-b", "--bin-size", dest="bin_size", \
                help="Bin size", metavar="SIZE", type="int")
        parser.add_option("-o", "--output", dest="output_dir", \
                help="Output directory", metavar="DIR",\
                default="data/output/breakpoint/")

        (options, args) = parser.parse_args()
        self.options = options
        if options.reads:
            self.reads_dir = options.reads
        if options.bin_size:
            self.bin_size = options.bin_size
        else:
            self.bin_size = 500
        #if options.reads_two:
        #    self.singleton_halves.append(options.reads_two)
        if not options.assembly_file:
            warning("Did not provide assembly file")
            parser.print_help()
            exit(1)
        else:
            self.assembly_file = options.assembly_file
        warning(self.singleton_halves)

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def warning(*objs):
    print("\tWARNING: ",*objs, file=sys.stderr)

def out_cmd(*objs):
    print("="*75, file=sys.stderr)
    print("About to exec: ", *objs, file=sys.stderr)

def main():
    '''
    splits read files for breakpoint
    '''
    finder = BreakpointFinder()
    finder.go()

if __name__=='__main__':
    main()

