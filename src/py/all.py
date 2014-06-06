#!/usr/bin/env python
from __future__ import print_function
from subprocess import call
from optparse import OptionParser
from tempfile import mkstemp
import os
import random
import shlex
import subprocess
import sys
import time

FNULL = open('/dev/null', 'w')

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def main():
    #Future echo implementation
    #command_echo_file = open(options.output_dir + "/commands.sh",'w')
    #command_echo_file.write("#!/bin/bash")
    
    (options, args) = get_options()

    start_time = time.time()
    fasta_file = options.fasta_file #fasta_file = "data/input/pairs/soap.trimmed.fna"

    error_files = []
    
    #reads_untrimmed_location = ["data/input/pairs/SRS011061.denovo_duplicates_marked\
    #        .trimmed.1.fastq","data/input/pairs/SRS011061.denovo_duplicates_marked.\
    #        trimmed.2.fastq"]
    reads_untrimmed_location = [options.first_mates, options.second_mates]
    reads_trimmed_location = []
    
    sam_output_location_dir = options.output_dir + "/sam/"
    sam_output_location = sam_output_location_dir + "library.sam"
    
    singleton_output_dir = options.output_dir + "/singleton/"
    singleton_output_location = singleton_output_dir + "singletons.csv"

    ensure_dir(sam_output_location_dir)
    ensure_dir(singleton_output_dir)

    step("ALIGNING READS")
    run_bowtie2(options, sam_output_location)

    step("CALCULATING ASSEMBLY PROBABILITY")
    run_lap(options, sam_output_location, reads_trimmed_location)

    outputBreakpointDir = options.output_dir + "/breakpoint/"
    ouputBreakpointLocation = outputBreakpointDir + "errorsDetected.csv"
    ensure_dir(outputBreakpointDir)

    ##breakpoint --- takes too long with current implementation
    ##alpha is length that must match
    #call(["bin/assembly-testing/breakpoint-detection/breakpoint_indices.py", "-a", fasta_file, "-u", singleton_output_location, "-o", ouputBreakpointLocation ,  "--alpha", "20","--algorithm","naive"])

    step("RUNNING SAMTOOLS")
    bam_location, sorted_bam_location, pileup_file = run_samtools(options, sam_output_location)

    step("DEPTH OF COVERAGE")
    error_files.append(run_depth_of_coverage(options, pileup_file))

    step("MATE-PAIR HAPPINESS")
    error_files.append(run_reapr(options, sorted_bam_location))

    step("SUMMARY")
    summary_file = open(options.output_dir + "/summary.gff", 'w')
    misassemblies = []
    for error_file in error_files:
        for line in open(error_file, 'r'):
            misassemblies.append(line.strip().split())

    # Sort misassemblies by start site.
    for misassembly in sorted(misassemblies, key = lambda misassembly: int(misassembly[3])):
        summary_file.write('\t'.join(misassembly) + '\n')

    summary_file.close()
    results(options.output_dir + "/summary.gff")

    if options.email:
        notify_complete(options.email,time.time()-start_time)


def get_options():
    parser = OptionParser()
    parser.add_option("-a", "--assembly-fasta", dest="fasta_file", \
            help="Candidate assembly file", metavar="FILE")
    parser.add_option("-r", "--reads", dest="reads_filenames", \
            help="First Read File", metavar="FILE")
    parser.add_option("-1", "--1", dest="first_mates", \
            help="Fastq filenames separated by commas that contain the first mates.")
    parser.add_option("-2", "--2", dest="second_mates", \
            help="Fastq filenames separated by commas that contain the second mates.")
    parser.add_option("-c", "--coverage-file", dest="coverage_file", \
            help="Assembly created per-contig coverage file")
    parser.add_option("-o", "--output-dir", dest="output_dir", \
            help = "Output directory", default="data/output/")
    parser.add_option("-w", "--window-size", dest="window_size", \
            help = "Sliding window size when determining misassemblies.", default = "100")
    parser.add_option("-q", "--fastq", dest="fastq_file", \
            default=False, action='store_true', \
            help="if set, input reads are fastq format (fasta by default).")    
    parser.add_option("-p", "--threads", dest="threads", \
            help = "Number of threads", default="8")
    parser.add_option("-I", "--minins", dest="min_insert_size", \
            help="Min insert sizes for mate pairs separated by commas.", default="0")
    parser.add_option("-X", "--maxins", dest="max_insert_size", \
            help="Max insert sizes for mate pairs separated by commas.", default="500")
    parser.add_option("-n", "--orientation", dest="orientation", default="fr", \
            help="Orientation of the mates.")
    parser.add_option("-m", "--mu" , dest="mu", default = "180", \
            help="average mate pair insert sizes.")
    parser.add_option("-t", "--sigma" , dest="sigma", default = "18", \
            help="standard deviation of mate pair insert sizes.")
    parser.add_option("-x", "--max_alignments", dest="max_alignments", default = "10000", \
            help="bowtie2 parameter to set the max number of alignments.")
    parser.add_option("-e", "--email", dest="email", \
            help="Email to notify when job completes")

    (options, args) = parser.parse_args()

    should_err = False
    if not options.fasta_file:
        warning("You need to provide a fasta file with -a")
        should_err = True
    #if not options.read_file_1:
    #    warning("You need to provide the first read file with -r")
    #    should_err = True
    #if not options.read_file_2:
    #    warning("You need to provide the second read file with -d")
    #    should_err = True
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

def step(*objs):
    line(75)
    print(bcolors.HEADER + "STEP:\t" + bcolors.ENDC, *objs, file=sys.stderr)

def out_cmd(*objs):
    #line(75)
    print(bcolors.OKBLUE + "COMMAND:\t" + bcolors.ENDC, ' '.join(*objs), file=sys.stderr)

def results(*objs):
    print(bcolors.WARNING + "RESULTS:\t" + bcolors.ENDC,*objs, file=sys.stderr)

def warning(*objs):
    print("INFO:\t",*objs, file=sys.stderr)


def build_bowtie2_index(index_name, reads_file):
    """
    Build a Bowtie2 index.
    """
    command = "bin/bowtie2-2.2.2/bowtie2-build " + os.path.abspath(reads_file) + " " + os.path.abspath(index_name)

    # Bad workaround.
    out_cmd([command])

    bowtie2_build_proc = subprocess.Popen(command, shell = True, stdout = FNULL, stderr = FNULL)
    bowtie_output, err = bowtie2_build_proc.communicate()

    return index_name


def run_bowtie2(options = None, output_sam = 'temp.sam'):
    """
    Run Bowtie2 with the given options and save the SAM file.
    """

    # Using bowtie2.
    # Create the bowtie2 index if it wasn't given as input.
    #if not assembly_index:
    if not os.path.exists(os.path.abspath(options.output_dir) +'/indexes'):
        os.makedirs(os.path.abspath(options.output_dir)+'/indexes')
    fd, index_path = mkstemp(prefix='temp_', dir=(os.path.abspath(options.output_dir)+'/indexes/'))
    try:
        os.mkdir(os.path.dirname(index_path))
    except:
        pass
    
    build_bowtie2_index(os.path.abspath(index_path), os.path.abspath(options.fasta_file))
    assembly_index = os.path.abspath(index_path)

    #input_sam_file = output_sam_file
    read_type = " -f "
    if options.fastq_file:
        read_type = " -q "
    
    bowtie2_args = ""
    if options.first_mates:
        bowtie2_args = "-a -x " + assembly_index + " -1 " + options.first_mates + " -2 " + options.second_mates + \
                " -p " + options.threads + " --very-sensitive -a " + " --reorder --" + options.orientation + " -I " + options.min_insert_size + \
                " -X " + options.max_insert_size
    else:
        bowtie2_args = "-a -x " + assembly_index + read_type + " -U " + options.reads_filenames + \
                " --very-sensitive -a " + " --reorder -p " + options.threads

    if not options:
        sys.stderr.write("[ERROR] No Bowtie2 options specified" + '\n')
        return
    
    # Using bowtie 2.
    command = "bin/bowtie2-2.2.2/bowtie2 " + bowtie2_args + " -S " + output_sam
    
    out_cmd([command])

    ignore = open('/dev/null', 'w')
    args = shlex.split(command) 
    bowtie_proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
            stderr=ignore)
    bowtie_output, err = bowtie_proc.communicate()


def run_lap(options, sam_output_location, reads_trimmed_location):
    """ Calculate the LAP using the previously computed SAM file. """

    output_probs_dir = options.output_dir + "/lap/"
    ensure_dir(output_probs_dir)
    output_probs_location = output_probs_dir + "output.prob"
    fp = open(output_probs_location, "w")

    reads = [options.reads_filenames]
    if options.first_mates:
        reads = [options.first_mates, options.second_mates]

    call_arr = ["bin/lap/aligner/calc_prob.py", "-a", options.fasta_file,  "-s", sam_output_location,  "-q", "-i",  ','.join(reads)]
    out_cmd(call_arr)
    #warning("That command outputs to: ", output_probs_location)
    results(output_probs_location)

    call(call_arr, stdout=fp)
    output_sum_probs_location = output_probs_dir + "output.sum"

    call_arr = ["bin/lap/aligner/sum_prob.py", "-i", output_probs_location]
    out_cmd(call_arr)
    call(call_arr, stdout=open(output_sum_probs_location,'w'))
    results(output_sum_probs_location)


def run_samtools(options, sam_output_location):
    """ Takes a sam file and runs samtools to create bam, sorted bam, and mpileup. """

    bam_dir = options.output_dir + "/bam/"
    ensure_dir(bam_dir)
    bam_location = bam_dir + "library.bam"
    sorted_bam_location = bam_dir + "sorted_library"
    bam_fp = open(bam_location, 'w+')

    #warning("About to run samtools view to create bam")
    call_arr = ["bin/Reapr_1.0.17/src/samtools", "view", "-bS", sam_output_location]
    out_cmd(call_arr)
    #warning("That command outputs to file: ", bam_location)
    call(call_arr, stdout = bam_fp, stderr = FNULL)

    #warning("About to attempt to sort bam")
    call_arr = ["bin/Reapr_1.0.17/src/samtools", "sort", bam_location, sorted_bam_location]
    out_cmd(call_arr)
    call(call_arr, stderr = FNULL)

    coverage_file_dir = options.output_dir + "/coverage/"
    ensure_dir(coverage_file_dir)
    pileup_file = coverage_file_dir + "mpileup_output.out"
    p_fp = open(pileup_file, 'w')
    call_arr = ["bin/Reapr_1.0.17/src/samtools", "mpileup", "-A", "-f", options.fasta_file, sorted_bam_location + ".bam"]
    out_cmd(call_arr)
    results(pileup_file)
    #warning("That command outputs to file: ", pileup_file)
    call(call_arr, stdout = p_fp, stderr = FNULL)

    return (bam_location, sorted_bam_location, pileup_file)


def run_depth_of_coverage(options, pileup_file):
    """ Run depth of coverage. """

    dp_fp = options.output_dir + "/coverage/errors_cov.gff"
    abundance_file = options.coverage_file
    call_arr = ["src/py/depth_of_coverage.py", "-a", abundance_file, "-m", pileup_file, "-w", options.window_size, "-o", dp_fp, "-g"]
    out_cmd(call_arr)
    call(call_arr)
    results(dp_fp)

    return dp_fp


def run_reapr(options, sorted_bam_location):
    """ Run REAPR. """

    reapr_command = "bin/Reapr_1.0.17/reapr"

    #warning("About to run facheck")
    call_arr = [reapr_command, "facheck", options.fasta_file ]
    out_cmd(call_arr)
    call(call_arr)

    reapr_output_dir = options.output_dir + "/reapr"
    reapr_perfect_prefix = options.output_dir + "/r_perfect_prefix"
    
    #warning("About to run reapr pipeline")
    call_arr = [reapr_command, "pipeline", options.fasta_file, sorted_bam_location + ".bam", reapr_output_dir]
    out_cmd(call_arr)
    call(call_arr, stdout=FNULL)

    call_arr = ["gunzip", reapr_output_dir + "/03.score.errors.gff"]
    out_cmd(call_arr)
    call(call_arr, stdout=FNULL)

    return reapr_output_dir + "/03.score.errors.gff"


if __name__ == '__main__':
    main()
