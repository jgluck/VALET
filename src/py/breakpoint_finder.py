#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
from subprocess import call
import os
import sys


class BreakpointFinder:

    def __init__(self):
        self.options = None
        self.assembly_file = None
        self.singleton_halves = []
        self.lengths = {}
        self.set_locations()
        self.getOptions()

    def set_locations(self):
        self.bowtie_dir = "/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/"
        self.bowtie_loc = self.bowtie_dir + "bowtie2"
        self.bowtie_build_loc = self.bowtie_dir + "bowtie2-build"
        self.breakpoint_dir = "data/output/breakpoint/"

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

    def run_bowtie_index(self):
        call_arr = [self.bowtie_build_loc, self.assembly_file, self.index_prefix]
        out_cmd(call_arr)
        call(call_arr)

    def run_bowtie_2_archived(self): #For mate pairs, not what we wanted
        call_arr = [self.bowtie_loc, "-x", self.index_prefix, "-1", \
                self.singleton_halves[0], "-2", self.singleton_halves[1],\
                "-S", self.sam_output_location, "--un-conc", self.singleton_loc,\
                "--al-conc", self.conc_loc, "-q", "-I 50", "-X", "1500", "-p 10"]
        out_cmd(call_arr)
        call(call_arr)

    def run_bowtie_2(self):
        for file_name in os.listdir(self.reads_dir):
            if "csv" in file_name:
                call_arr = [self.bowtie_loc, '-x', self.index_prefix, '-U',\
                         self.reads_dir + file_name,\
                         '-S', self.sam_output_dir + file_name + '.sam',\
                         '--un', self.singleton_dir + file_name + '.singletons',\
                         '--al', self.conc_dir + file_name + '.reads',
                         '-q', '-I 50', '-X 800', '-p 10']
                out_cmd(call_arr)
                call(call_arr)
            else:
                warning("Skipping potential read file: " + file_name)

    def read_lengths(self):
        for file_name in os.listdir(self.conc_dir):
            if '.reads' in file_name:
                with open(self.conc_dir + file_name, 'r') as reads_file:
                    for read in read_read(reads_file):
                        print(read)
 
    def read_read(self, fp):
        while True:
            line_bundle = []
            for i in range(4):
                line_bundle.append(fp.readline())
            if len(line_bundle) == 0:
                break
            yield line_bundle

    def detect_breakpoints(self):
        with open(self.breakpoint_file, 'w') as out_file:
            for file_name in os.listdir(self.sam_output_dir):
                if "sam" in file_name:
                    with open(self.sam_output_dir+file_name, 'r') as in_file:
                        for line in in_file:
                            if line[0] != '@':
                                line_components = line.split('\t')
                                if line_components[3] != '0':
                                    out_file.write("%s\t%s\t%s\n" % (line_components[3],\
                                            line_components[0], line_components[1]))

    def go(self):
        self.run_bowtie_index()
        self.run_bowtie_2()
        self.detect_breakpoints()

    def getOptions(self):
        parser = OptionParser()
        parser.add_option("-a", "--assembly-file", dest="assembly_file",\
                help="Assembly File to search", metavar="FILE")
        parser.add_option("-r", "--reads-dir", dest="reads", \
                help="Directory of Reads", metavar="PATH")
        #parser.add_option("-l", "--reads-2", dest="reads_two", \
        #        help="Second read file to split", metavar="FILE")

        (options, args) = parser.parse_args()
        self.options = options
        if options.reads:
            self.reads_dir = options.reads
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

