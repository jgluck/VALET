#!/usr/bin/env python
from __future__ import print_function
from subprocess import call
from optparse import OptionParser
import os
import time
import sys

def main():
    #Future echo implementation
    #command_echo_file = open("data/output/commands.sh",'w')
    #command_echo_file.write("#!/bin/bash")
    
    (options, args) = get_options()

    start_time = time.time()
    fastaFile = options.fasta_file #fastaFile = "data/input/pairs/soap.trimmed.fna"
    outputPrefixDir = "data/bowtieIndex/"
    outputPrefix = outputPrefixDir + "soap"
    
    #reads_untrimmed_location = ["data/input/pairs/SRS011061.denovo_duplicates_marked\
    #        .trimmed.1.fastq","data/input/pairs/SRS011061.denovo_duplicates_marked.\
    #        trimmed.2.fastq"]
    reads_untrimmed_location = [options.read_file_1, options.read_file_2]
    reads_trimmed_location = []
    
    samOutputLocationDir = "data/output/sam/"
    samOutputLocation = samOutputLocationDir + "soap.sam"
    
    singletonOutputDir = "data/output/singleton/"
    singletonOutputLocation = singletonOutputDir + "singletons.csv"

    ensure_dir(outputPrefixDir)
    ensure_dir(singletonOutputDir)
    

    #TODO: Need to fix Chris bowtie2 alias. Install locally
    call_string = "/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/bowtie2-build %s %s"%(fastaFile, outputPrefix)
    call_arr = ["/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/bowtie2-build", fastaFile, outputPrefix]
    out_cmd(call_arr)
    call(call_arr) 

    ensure_dir(samOutputLocationDir)

    should_align = True #Should be true but not right now
    min_read_length = 300


    if len(reads_untrimmed_location)>1 and should_align:
        for readFile in reads_untrimmed_location:
            reads_trimmed_location.append(readFile) #Added for interim, until trimming actually works
        warning("About to run bowtie2")
        call_arr = ["/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/bowtie2", "-x", outputPrefix, "-1", reads_trimmed_location[0], "-2", reads_trimmed_location[1], "-S", samOutputLocation, "--un", singletonOutputLocation, "-q", "-I 100", "-X 500", "-p 8"]
        out_cmd(call_arr)
        call(call_arr)
    elif should_align:
        reads_trimmed_location.append(reads_untrimmed_location[0])
        call_arr = ["/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/bowtie2", "-x", outputPrefix, "-U", reads_trimmed_location[0], "-S",  samOutputLocation, "--un", singletonOutputLocation, "-p 8"]
        out_cmd(call_arr)
        call(call_arr)

    outputProbsDir = "data/output/lap/"
    ensure_dir(outputProbsDir)
    outputProbsLocation = outputProbsDir + "output.prob"
    fp = open(outputProbsLocation, "w")

    ##lap
    call_arr = ["bin/lap/aligner/calc_prob.py", "-a", fastaFile,  "-s", samOutputLocation,  "-q", "-i",  reads_trimmed_location[0]]
    out_cmd(call_arr)
    warning("That command outputs to: ", outputProbsLocation)
    call(call_arr, stdout=fp)
    call(["bin/lap/aligner/sum_prob.py", "-i", outputProbsLocation])


    outputBreakpointDir = "data/output/breakpoint/"
    ouputBreakpointLocation = outputBreakpointDir + "errorsDetected.csv"
    ensure_dir(outputBreakpointDir)

    ##breakpoint --- takes too long with current implementation
    ##alpha is length that must match
    #call(["bin/assembly-testing/breakpoint-detection/breakpoint_indices.py", "-a", fastaFile, "-u", singletonOutputLocation, "-o", ouputBreakpointLocation ,  "--alpha", "20","--algorithm","naive"])


    reapr_command = "bin/Reapr_1.0.17/reapr"

    warning("About to run facheck")
    call_arr = [reapr_command, "facheck", fastaFile ]
    out_cmd(call_arr)
    call(call_arr)

    bamDir = "data/output/bam/"
    ensure_dir(bamDir)
    bamLocation = bamDir + "soap.bam"
    sorted_bam_location = bamDir+"sorted_soap"
    bam_fp = open(bamLocation, 'w+')

    warning("About to run samtools view to create bam")
    call_arr = ["samtools", "view", "-bS", samOutputLocation]
    out_cmd(call_arr)
    warning("That command outputs to file: ", bamLocation)
    call(call_arr,stdout=bam_fp)

    warning("About to attempt to sort bam")
    call_arr = ["samtools", "sort", bamLocation,sorted_bam_location]
    out_cmd(call_arr)
    call(call_arr)

    coverage_file_dir = "data/output/coverage/"
    ensure_dir(coverage_file_dir)
    pileup_file = coverage_file_dir + "mpileup_output.out"
    p_fp = open(pileup_file, 'w')
    call_arr = ["samtools", "mpileup", "-f", fastaFile, sorted_bam_location]
    out_cmd(call_arr)
    warning("That command outputs to file: ", pileup_file)
    call(call_arr,stdout=p_fp)

    abundance_file = options.coverage_file #"data/input/pairs/soap.abundance.cvg"
    call_arr = ["src/py/depth_of_coverage.py", "-a", abundance_file, "-m", pileup_file]
    out_cmd(call_arr)
    call(call_arr)

    reapr_output_dir = "data/output/reapr"
    reapr_perfect_prefix = "data/output/r_perfect_prefix"
    
    warning("About to run reapr pipeline")
    call_arr = [reapr_command,"pipeline",fastaFile,sorted_bam_location+".bam",reapr_output_dir]
    out_cmd(call_arr)
    call(call_arr)
    
    if options.email:
        notify_complete(options.email,time.time()-start_time)


def get_options():
    parser = OptionParser()
    parser.add_option("-a", "--assembly-fasta", dest="fasta_file", \
            help="Candidate assembly file", metavar="FILE")
    parser.add_option("-r", "--reads-1", dest="read_file_1", \
            help="First Read File", metavar="FILE")
    parser.add_option("-d", "--reads-2", dest="read_file_2", \
            help="Second Read File", metavar="FILE")
    parser.add_option("-c", "--coverage-file", dest="coverage_file", \
            help="Assembly created per-contig coverage file")
    parser.add_option("-e", "--email", dest="email", \
            help="Email to notify when job completes")

    (options, args) = parser.parse_args()

    should_err = False
    if not options.fasta_file:
        warning("You need to provide a fasta file with -a")
        should_err = True
    if not options.read_file_1:
        warning("You need to provide the first read file with -r")
        should_err = True
    if not options.read_file_2:
        warning("You need to provide the second read file with -d")
        should_err = True
    if not options.coverage_file:
        warning("You need to provide the coverage file with -c")
        should_err = True
    
    if should_err:
        parser.print_help()
        exit(-1)

    return (options,args)

def ran_command(st, fp):
    fp.write(st)

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def notify_complete(target_email,t):
    call(['echo "Completed in %d" | mail -s "Job Completed" %s' % (t, target_email) ],shell=True)

def line(x):
    print ("-"*x)

def out_cmd(*objs):
    line(75)
    print("About to exec: ", *objs, file=sys.stderr)

def warning(*objs):
    print("\tWARNING: ",*objs, file=sys.stderr)

if __name__ == '__main__':
    main()
