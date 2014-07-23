#!/usr/bin/env python
from __future__ import print_function
from collections import defaultdict
from subprocess import call
from optparse import OptionParser
from tempfile import mkstemp
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import time

FNULL = open('/dev/null', 'w')
base_path = os.path.dirname(sys.argv[0])[:-len('src/py/')]

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def main():
   
    (options, args) = get_options()

    start_time = time.time()
    fasta_file = options.fasta_file

    error_files = []
    
    reads_untrimmed_location = [options.first_mates, options.second_mates]
    reads_trimmed_location = []
    
    sam_output_location_dir = options.output_dir + "/sam/"
    sam_output_location = sam_output_location_dir + "library.sam"
    
    singleton_output_dir = options.output_dir + "/singleton/"
    singleton_output_location = singleton_output_dir + "singletons.csv"

    ensure_dir(sam_output_location_dir)
    ensure_dir(singleton_output_dir)
    
    bins_dir = options.output_dir + "/bins/"
    ensure_dir(bins_dir)
 
    input_fasta_saved = options.fasta_file
    output_dir_saved = options.output_dir
    
    all_contig_lengths = {}

    if options.min_contig_length > 0:
        step("FILTERING ASSEMBLY CONTIGS LESS THAN " + str(options.min_contig_length) + ' BPs')
        all_contig_lengths = filter_short_contigs(options)
        results(options.fasta_file)
        fasta_file = options.fasta_file
        input_fasta_saved = options.fasta_file

    step("ALIGNING READS")
    unaligned_dir = run_bowtie2(options, sam_output_location)

    contig_lengths = get_contig_lengths(sam_output_location)

    step("RUNNING SAMTOOLS")
    bam_location, sorted_bam_location, pileup_file = \
            run_samtools(options, sam_output_location, index=True)

    if options.coverage_file is None:
        step("CALCULATING CONTIG COVERAGE")
        options.coverage_file = calculate_contig_coverage(options, pileup_file)
        results(options.coverage_file)

    step("CALCULATING ASSEMBLY PROBABILITY")
    run_lap(options, sam_output_location, reads_trimmed_location)

    step("DEPTH OF COVERAGE")
    error_files.append(run_depth_of_coverage(options, pileup_file))

    contig_to_bin_map, bin_dir_dict = bin_coverage(options,bins_dir)
    split_sam_by_bin(sam_output_location, contig_to_bin_map, bin_dir_dict)


    outputBreakpointDir = options.output_dir + "/breakpoint/"
    ouputBreakpointLocation = outputBreakpointDir + "errorsDetected.csv"
    ensure_dir(outputBreakpointDir)
    
    step("BREAKPOINT")
    error_files.append(run_breakpoint_finder(options,\
            unaligned_dir, outputBreakpointDir))

    for bin_dir in os.listdir(bins_dir):
        #if 'bin' in bin_dir:
        coverages = bin_dir
        options.fasta_file = os.path.abspath(output_dir_saved) + '/bins/'\
                + bin_dir + '/' + os.path.basename(input_fasta_saved)

        options.output_dir = os.path.abspath(output_dir_saved) + '/bins/'\
                + bin_dir

        bin_dir_infix = '/bins/' + bin_dir + '/'
        bin_dir = os.path.abspath(options.output_dir) + '/bins/' + bin_dir + '/'
        #warning("Bin dir is: %s" % bin_dir)
        sam_output_location_dir = options.output_dir + '/sam/'
        sam_output_location = sam_output_location_dir + 'library.sam'
        
        step("RUNNING SAMTOOLS ON COVERAGE BIN " + coverages)
        bam_location, sorted_bam_location, pileup_file = \
                run_samtools(options, sam_output_location, with_pileup=False)

        #step("DEPTH OF COVERAGE")
        #error_files.append(run_depth_of_coverage(options, pileup_file))
        
        step("MATE-PAIR HAPPINESS ON COVERAGE BIN " + coverages)
        try:
            error_files.append(run_reapr(options, sorted_bam_location))
        except:
            e = sys.exc_info()[0]
            error("Reapr failed to run with: %s" %  str(e))
    
    options.output_dir = output_dir_saved
    options.fasta_file = input_fasta_saved

    step("SUMMARY")
    summary_file = open(options.output_dir + "/summary.gff", 'w')
    suspicious_file = open(options.output_dir + "/suspicious.gff", 'w')
    summary_table_file = open(options.output_dir + "/summary.tsv", 'w')
    #suspicious_table_file = open(options.output_dir + "/suspicious.tsv", 'w')

    misassemblies = []
    for error_file in error_files:
        if error_file:
            for line in open(error_file, 'r'):
                misassemblies.append(line.strip().split('\t'))

    # Sort misassemblies by start site.
    misassemblies.sort(key = lambda misassembly: (misassembly[0], int(misassembly[3])))
    final_misassemblies = []
    for misassembly in misassemblies:

        # Truncate starting/ending region if it is near the end of the contigs.
        if int(misassembly[3]) <= options.ignore_end_distances and \
            int(misassembly[4]) > options.ignore_end_distances:
          misassembly[3] = str(options.ignore_end_distances + 1)

        if int(misassembly[4]) >= (contig_lengths[misassembly[0]] - options.ignore_end_distances) and \
            int(misassembly[3]) < (contig_lengths[misassembly[0]] - options.ignore_end_distances):
          misassembly[4] = str(contig_lengths[misassembly[0]] - options.ignore_end_distances - 1)

        # Don't print a flagged region if it occurs near the ends of the contig.
        if int(misassembly[3]) > options.ignore_end_distances and \
                int(misassembly[4]) < (contig_lengths[misassembly[0]] - options.ignore_end_distances):
            summary_file.write('\t'.join(misassembly) + '\n')

            final_misassemblies.append(misassembly)
    
    summary_file.close()

    results(options.output_dir + "/summary.gff")

    # Find regions with multiple misassembly signatures.
    suspicious_regions = find_suspicious_regions(misassemblies, options.min_suspicious_regions)
    final_suspicious_misassemblies = []
    for region in suspicious_regions:
        if int(region[3]) > options.ignore_end_distances and \
                int(region[4]) <= (contig_lengths[region[0]] - options.ignore_end_distances):
            suspicious_file.write('\t'.join(region) + '\n')
            final_suspicious_misassemblies.append(region)

    results(options.output_dir + "/suspicious.gff")

    # Output summary table.
    generate_summary_table(options.output_dir + "/summary.tsv", all_contig_lengths, \
        contig_lengths, final_misassemblies)

    results(options.output_dir + "/summary.tsv")

    # Output suspicious table.
    #generate_summary_table(options.output_dir + "/suspicious.tsv", all_contig_lengths, \
    #    contig_lengths, final_suspicious_misassemblies)

    #results(options.output_dir + "/suspicious.tsv")

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
            help = "Sliding window size when determining misassemblies.", default = "200")
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
    parser.add_option("-x", "--max-alignments", dest="max_alignments", default = "10000", \
            help="bowtie2 parameter to set the max number of alignments.")
    parser.add_option("-e", "--email", dest="email", \
            help="Email to notify when job completes")
    parser.add_option("-g", "--min-coverage", dest="min_coverage", type="int", default=0, \
            help="Minimum average coverage to run misassembly detection.")
    parser.add_option("-l", "--coverage-multiplier", dest="coverage_multiplier", type=float, default=0.0, \
            help="When binning by coverage, the new high = high + high * multiplier")
    parser.add_option("-s", "--min-suspicious", dest="min_suspicious_regions", default=2, type=int, \
            help="Minimum number of overlapping flagged miassemblies to mark region as suspicious.")
    parser.add_option('-z', "--min-contig-length", dest="min_contig_length", default=1000, type=int, \
            help="Ignore contigs smaller than this length.")
    parser.add_option('-b', "--ignore-ends", dest="ignore_end_distances", default=0, type=int, \
            help="Ignore flagged regions within b bps from the ends of the contigs.") 
    parser.add_option('-k', "--breakpoint-bin", dest="breakpoints_bin", default="50", type=str, \
            help="Bin sized used to find breakpoints.")   

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
        warning("Coverage file not provided, will create one.")
        
        #should_err = True
    
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
    print ("-"*x, file=sys.stderr)

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


def error(*objs):
    print(bcolors.WARNING + "ERROR:\t" + bcolors.ENDC, *objs, file=sys.stderr)


def filter_short_contigs(options):
    """
    Filter out contigs less than a certain length.
    """

    filtered_fasta_filename = options.output_dir + '/filtered_assembly.fasta'
    filtered_assembly_file = open(filtered_fasta_filename, 'w')

    all_contig_lengths = {}

    curr_length = 0
    with open(options.fasta_file,'r') as assembly:
        for contig in contig_reader(assembly):
            curr_length = len(''.join(contig['sequence']))
            
            if curr_length >= options.min_contig_length:
                filtered_assembly_file.write(contig['name'])
                filtered_assembly_file.writelines(contig['sequence'])
                filtered_assembly_file.write('\n')

            all_contig_lengths[contig['name'].strip()[1:]] = curr_length

    filtered_assembly_file.close()
    options.fasta_file = filtered_fasta_filename

    return all_contig_lengths


def get_contig_lengths(sam_filename):
    """
    Return a dictionary of contig names => contig lengths from a SAM file.
    """

    sam_file = open(sam_filename, 'r')

    # Build dictionary of contig lengths.
    contig_lengths = {}
    pattern = re.compile('SN:(?P<contig>[\w_\|\.]+)\s*LN:(?P<length>\d+)')
    line = sam_file.readline()
    while line.startswith("@"):

        if line.startswith("@SQ"):
            matches = pattern.search(line)
            
            if len(matches.groups()) == 2:
                contig_lengths[matches.group('contig')] = int(matches.group('length'))

        line = sam_file.readline()

    return contig_lengths


def find_suspicious_regions(misassemblies, min_cutoff = 2):
    """
    Given a list of miassemblies in gff format
    """

    regions =[]

    for misassembly in misassemblies:
        regions.append([misassembly[0], misassembly[3], 'START', misassembly[2]])
        regions.append([misassembly[0], misassembly[4], 'END', misassembly[2]])

    regions.sort(key = lambda region: (region[0], int(region[1])))

    """
    Example:

    relocref        36601   START   Breakpoint_finder                                    
    relocref        36801   END     Breakpoint_finder                                    
    relocref        67968   START   REAPR
    relocref        68054   START   REAPR
    relocref        69866   END     REAPR
    relocref        69867   START   REAPR
    relocref        71833   END     REAPR
    relocref        73001   START   Breakpoint_finder                                    
    relocref        73201   END     Breakpoint_finder   
    """

    curr_contig = None
    curr_index = 0
    curr_length = -1

    start_indexes = []
    start_region = 0
    end_index = 0
    signatures = []
    recording = False

    signature_starts = defaultdict(list)

    curr_coverage = 0
    suspicious_regions = []

    for region in regions:

        if curr_contig is None:
            curr_contig = region[0]
            recording = False
            signature_starts = defaultdict(list)

        # We have found a new contig, process the previous contig results.
        if region[0] != curr_contig:

            curr_contig = region[0]
            recording = False

        if region[2] == 'START':
            curr_coverage += 1
            if region[3] not in signatures: signatures.append(region[3])
            signature_starts[region[3]].append(region[1])

            # Record start point.
            if curr_coverage == min_cutoff:
                start_region = region[1]
                recording == True

            start_indexes.append(region[1])

        else:
            
            curr_coverage -= 1

            end_index = region[1]
            if region[3] in signatures: signatures.remove(region[3]) 

            # If we were recording, and min signatures drop belows threshold,
            # then we need to output our results
            if curr_coverage < min_cutoff and recording:
                min_start = None

                suspicious_regions.append([region[0], '.', 'SUSPICIOUS', str(start_region), str(end_index), '.', '.', '.', 'color=#181009;' + ','.join(signatures)])
                signatures = []
                recording = False

        if curr_coverage >= min_cutoff:
            recording = True

    # Hack to correct for overlapping suspicious regions.
    compressed_suspicious_regions = []

    prev_region = None
    for region in suspicious_regions:

        if prev_region is None:
            prev_region = region
        else:
            if prev_region[0] == region[0] and int(prev_region[4]) >= int(region[3]):
                prev_region[4] = region[4]
            else:
                compressed_suspicious_regions.append(prev_region)
                prev_region = region


    if prev_region:
        compressed_suspicious_regions.append(prev_region)

    return compressed_suspicious_regions


def generate_summary_table(table_filename, all_contig_lengths, filtered_contig_lengths, misassemblies):
    """
    Output the misassemblies in a table format:

    contig_name  contig_length  low_cov  low_cov_bps  high_cov  high_cov_bps ...
    CONTIG1 12000   1   100 0   0 ...
    CONTIG2 100 NA  NA  NA ...  
    """

    table_file = open(table_filename, 'w')
    table_file.write("contig_name\tcontig_length\tlow_cov\tlow_cov_bps\thigh_cov\thigh_cov_bps\treapr\treapr_bps\tbreakpoints\tbreakpoints_bps\n")

    prev_contig = None
    curr_contig = None

    # Misassembly signatures
    low_coverage = 0
    low_coverage_bps = 0
    high_coverage = 0   
    high_coverage_bps = 0
    reapr = 0
    reapr_bps = 0
    breakpoints = 0
    breakpoints_bps = 0

    processed_contigs = set()

    for misassembly in misassemblies:
        """
        contig00001     REAPR   Read_orientation        88920   97033   .       .       .       Note=Warning: Bad read orientation;colour=1
        contig00001     REAPR   FCD     89074   90927   0.546142        .       .       Note=Error: FCD failure;colour=17
        contig00001     DEPTH_COV       low_coverage    90818   95238   29.500000       .       .       low=30.000000;high=70.000000;color=#7800ef
        """

        curr_contig = misassembly[0]

        if prev_contig is None:
            prev_contig = curr_contig

        if curr_contig != prev_contig:
            # Output previous contig stats.
            table_file.write(prev_contig + '\t' + str(filtered_contig_lengths[prev_contig]) + '\t' + \
                str(low_coverage) + '\t' + str(low_coverage_bps) + '\t' + str(high_coverage) + '\t' + \
                str(high_coverage_bps) + '\t' + str(reapr) + '\t' + str(reapr_bps) + '\t' + str(breakpoints) + '\t' + \
                str(breakpoints_bps) + '\n')

            processed_contigs.add(prev_contig)

            # Reset misassembly signature counts.
            low_coverage = 0
            low_coverage_bps = 0
            high_coverage = 0   
            high_coverage_bps = 0
            reapr = 0
            reapr_bps = 0
            breakpoints = 0
            breakpoints_bps = 0

            prev_contig = curr_contig

        # Process the current contig misassembly.
        if misassembly[1] == 'REAPR':
            reapr += 1
            reapr_bps += (int(misassembly[4]) - int(misassembly[3]) + 1)

        elif misassembly[1] == 'DEPTH_COV':
            if misassembly[2] == 'Low_coverage':
                low_coverage += 1
                low_coverage_bps += (int(misassembly[4]) - int(misassembly[3]) + 1)
            else:
                high_coverage += 1
                high_coverage_bps += (int(misassembly[4]) - int(misassembly[3]) + 1)

        elif misassembly[1] == 'Breakpoint_finder':
            breakpoints += 1
            breakpoints_bps += (int(misassembly[4]) - int(misassembly[3]) + 1)

        else:
            print("Unhandled error: " + misassembly[1])

    if prev_contig:
        # Output previous contig stats.
        table_file.write(prev_contig + '\t' + str(filtered_contig_lengths[prev_contig]) + '\t' + \
            str(low_coverage) + '\t' + str(low_coverage_bps) + '\t' + str(high_coverage) + '\t' + \
            str(high_coverage_bps) + '\t' + str(reapr) + '\t' + str(reapr_bps) + '\t' + str(breakpoints) + '\t' + \
            str(breakpoints_bps) + '\n')

        processed_contigs.add(prev_contig)

    # We need to add the remaining, error-free contigs.
    for contig in filtered_contig_lengths:
        if contig not in processed_contigs:
            table_file.write(contig + '\t' + str(filtered_contig_lengths[contig]) + '\t' + \
                '0\t0\t0\t0\t0\t0\t0\t0\n')
            processed_contigs.add(contig)


    # Finally, add the contigs that were filtered out prior to evaluation.
    for contig in all_contig_lengths:
        if contig not in processed_contigs:
            table_file.write(contig + '\t' + str(all_contig_lengths[contig]) + '\t' + \
                'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n')
            processed_contigs.add(contig)


def calculate_contig_coverage(options, pileup_file):
    """
    Calculate contig coverage.  The coverage of a contig is the mean per-bp coverage.
    """

    coverage_filename = options.output_dir + '/coverage/temp.cvg'
    coverage_file = open(coverage_filename, 'w')

    prev_contig = None
    curr_contig = None

    length = 0
    curr_coverage = 0

    for record in open(pileup_file, 'r'):
        fields = record.strip().split()

        if prev_contig != fields[0]:
            if prev_contig:
                coverage_file.write(prev_contig + '\t' + str(float(curr_coverage) / length) + '\n')

            prev_contig = fields[0]
            length = 0
            curr_coverage = 0

        curr_coverage += int(fields[3])
        length += 1
    if prev_contig:
        coverage_file.write(prev_contig + '\t' + str(float(curr_coverage) / length) + '\n')
    coverage_file.close()

    return coverage_filename


def build_bowtie2_index(index_name, reads_file):
    """
    Build a Bowtie2 index.
    """
    command = os.path.join(base_path, "bin/bowtie2-2.2.2/bowtie2-build ") + os.path.abspath(reads_file) + " " + os.path.abspath(index_name)

    # Bad workaround.
    out_cmd([command])

    bowtie2_build_proc = subprocess.Popen(command, shell = True, stdout = FNULL, stderr = FNULL)
    bowtie_output, err = bowtie2_build_proc.communicate()
    bowtie2_build_proc.wait()

    return index_name


def run_bowtie2(options = None, output_sam = 'temp.sam'):
    """
    Run Bowtie2 with the given options and save the SAM file.
    """

    # Using bowtie2.
    # Create the bowtie2 index if it wasn't given as input.
    #if not assembly_index:
    if not os.path.exists(os.path.abspath(options.output_dir) + '/indexes'):
        os.makedirs(os.path.abspath(options.output_dir) + '/indexes')
    fd, index_path = mkstemp(prefix='temp_',\
            dir=(os.path.abspath(options.output_dir)   + '/indexes/'))
    try:
        os.mkdirs(os.path.dirname(index_path))
    except:
        pass
    
    fasta_file = options.fasta_file

    build_bowtie2_index(os.path.abspath(index_path), os.path.abspath(fasta_file))
    assembly_index = os.path.abspath(index_path)

    unaligned_dir = os.path.abspath(options.output_dir) + '/unaligned_reads/'
    ensure_dir(unaligned_dir)
    unaligned_file = unaligned_dir + 'unaligned.reads'

    #input_sam_file = output_sam_file
    read_type = " -f "
    if options.fastq_file:
        read_type = " -q "
    
    bowtie2_args = ""
    bowtie2_unaligned_check_args = ""
    if options.first_mates:
        bowtie2_args = "-a -x " + assembly_index + " -1 " + options.first_mates\
                + " -2 " + options.second_mates + " -p " + options.threads\
                + " --very-sensitive -a " + " --reorder --"\
                + options.orientation + " -I " + options.min_insert_size\
                + " -X " + options.max_insert_size #+ " --un-conc "\
                #+ unaligned_file
        
        bowtie2_unaligned_check_args = "-a -x " + assembly_index + read_type + " -U "\
                + options.first_mates + "," + options.second_mates + " --very-sensitive -a "\
                + " --reorder -p " + options.threads + " --un " + unaligned_file
 
    else:
        bowtie2_args = "-a -x " + assembly_index + read_type + " -U "\
                + options.reads_filenames + " --very-sensitive -a "\
                + " --reorder -p " + options.threads + " --un " + unaligned_file

    if not options:
        sys.stderr.write("[ERROR] No Bowtie2 options specified" + '\n')
        return
    
    # Using bowtie 2.
    command = os.path.join(base_path, "bin/bowtie2-2.2.2/bowtie2 ") + bowtie2_args + " -S " + output_sam
    
    out_cmd([command])


    ignore = open('/dev/null', 'w')
    #call(command.split())
    args = shlex.split(command) 
    bowtie_proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=ignore)
    bowtie_output, err = bowtie_proc.communicate()
 


    if bowtie2_unaligned_check_args != "":
        command = os.path.join(base_path, "bin/bowtie2-2.2.2/bowtie2 ") + bowtie2_unaligned_check_args + " -S " + output_sam + "_2.sam"
        out_cmd([command])
        args = shlex.split(command)
        bowtie_proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=ignore)
        bowtie_output, err = bowtie_proc.communicate()

    return unaligned_dir


def run_breakpoint_finder(options,unaligned,breakpoint_dir):
    '''
    attempts to find breakpoints
    '''
    std_err_file = open(breakpoint_dir + 'splitter_std_err.log', 'w')
    call_arr = [os.path.join(base_path,'src/py/breakpoint_splitter.py'),\
            '-u', unaligned,\
            '-o', breakpoint_dir + 'split_reads/']

    out_cmd(call_arr)
    call(call_arr, stderr=std_err_file)
    std_err_file.close()
    
    std_err_file = open(breakpoint_dir + 'std_err.log','w')
    call_arr = [os.path.join(base_path, 'src/py/breakpoint_finder.py'),\
            '-a', options.fasta_file,\
            '-r', breakpoint_dir + 'split_reads/',\
            '-b', options.breakpoints_bin, '-o', breakpoint_dir,\
            '-c', options.coverage_file]
    out_cmd(call_arr)
    call(call_arr,stderr=std_err_file)
    results(breakpoint_dir + 'interesting_bins.gff')
    return breakpoint_dir + 'interesting_bins.gff'


def split_sam_by_bin(sam_output_location, contig_to_bin_map, bin_dir_dict):
    common_header = ""
    output_bin = {}
    output_fp = {}
    for bin in set(contig_to_bin_map.values()):
        output_bin[bin] = ""
        bin_dir = bin_dir_dict[bin]
        if os.path.exists(bin_dir):
            ensure_dir(bin_dir + "sam/")
            output_fp[bin]  = open(bin_dir + "sam/"\
                    + os.path.basename(sam_output_location), 'w')

    with open(sam_output_location, 'r') as sam_file:
        for line in sam_file:
            if line.split()[0] == "@HD" or line.split()[0] == "@PG"\
                    or line.split()[0] == "@CO" or line.split()[0] == "@RG":
                        for fp in output_fp.values():
                            fp.write(line)
            elif line.split()[0] == "@SQ":
                # TODO: Clean up.
                if line.split()[1].split(':')[1] in contig_to_bin_map:
                    bin = contig_to_bin_map[line.split()[1].split(':')[1]]
                    output_fp[bin].write(line)
            else:
                line_split = line.split('\t')
                if line_split[2] == '*':
                    pass
                else:
                    # TODO: Clean up.
                    if line_split[2] in contig_to_bin_map:
                        bin = contig_to_bin_map[line_split[2]]
                        output_fp[bin].write(line)
    
        
def increment_coverage_window(options, low, high):
    """ Find new low/high boundaries for coverage bins. """

    low = high
    prev_high = high
    high = int(high + high * options.coverage_multiplier)
    if high == prev_high:
        high = high + 1

    return low, high


def bin_coverage(options, bin_dir):
    contig_to_coverage_map = {}
    contig_to_bin_map = {}
    bin_to_name_map = {}

    with open(options.coverage_file,'r') as coverage_file:
        for line in coverage_file:
            split_line = line.split()
            if float(split_line[1]) >= options.min_coverage:
                # Only store contigs who are above minimum avg coverage.
                contig_to_coverage_map[split_line[0]] = float(split_line[1])

    max_cvg = max(contig_to_coverage_map.values())

    high = int(options.min_coverage + options.min_coverage * .1)
    if high <= options.min_coverage:
        high = high + 1
    low = options.min_coverage

    curr_bin = 0
    bins = []
    while len(contig_to_bin_map.keys()) < len(contig_to_coverage_map.keys()):
        slice_dict = {k: v for k,v in contig_to_coverage_map.iteritems() if low<=v and high>v}
        for contig in slice_dict.keys():
            contig_to_bin_map[contig] = curr_bin

        bin_to_name_map[curr_bin] = (low, high)

        low, high = increment_coverage_window(options, low, high)

        curr_bin += 1

    bin_set = set(contig_to_bin_map.values())
    fp_dict = {}
    bin_dir_dict = {}
    for bin in bin_set:
        #a_new_bin = bin_dir + "bin" + str(bin) + "/"
        a_new_bin = bin_dir + str(bin_to_name_map[bin][0]) + "x-" + str(bin_to_name_map[bin][1]) + "x/"
        bin_dir_dict[bin] = a_new_bin
        ensure_dir(a_new_bin)
        shutil.copy(options.coverage_file, a_new_bin +\
                os.path.basename(options.coverage_file))
        fp_dict[bin] = open(a_new_bin + os.path.basename(options.fasta_file),'w')

    with open(options.fasta_file,'r') as assembly:
        for contig in contig_reader(assembly):
            # TODO: Clean up.
            if contig['name'][1:].strip() in contig_to_bin_map:
                bin = contig_to_bin_map[contig['name'][1:].strip()]
                fp_dict[bin].write(contig['name'])
                fp_dict[bin].writelines(contig['sequence'])

    for fp in fp_dict.values():
        fp.close()

    for fp in fp_dict.values():
        name = fp.name
        if os.stat(fp.name).st_size <= 10:
            shutil.rmtree(os.path.dirname(fp.name))

    return contig_to_bin_map,bin_dir_dict


def contig_reader(fasta_file):
    save_line = ""
    contig = {}
    in_contig = False
    for line in fasta_file:
        if line[0] == '>' and in_contig:
            save_line = line
            ret_contig = contig
            contig = {}
            contig['sequence'] = []
            contig['name'] = line
            yield ret_contig
        elif line[0] == '>':
            contig['name'] = line
            contig['sequence'] = []
            in_contig = True
        else:
            contig['sequence'].append(line.strip())
    yield contig



def run_lap(options, sam_output_location, reads_trimmed_location):
    """ Calculate the LAP using the previously computed SAM file. """
    output_probs_dir = options.output_dir + "/lap/"
    
    ensure_dir(output_probs_dir)
    output_probs_location = output_probs_dir + "output.prob"
    fp = open(output_probs_location, "w")

    reads = [options.reads_filenames]
    if options.first_mates:
        reads = [options.first_mates, options.second_mates]

    call_arr = [os.path.join(base_path, "bin/lap/aligner/calc_prob.py"), "-a", options.fasta_file,  "-s", sam_output_location,  "-q", "-i",  ','.join(reads), "-n", options.coverage_file]
    out_cmd(call_arr)
    #warning("That command outputs to: ", output_probs_location)
    results(output_probs_location)

    call(call_arr, stdout=fp)
    output_sum_probs_location = output_probs_dir + "output.sum"

    call_arr = [os.path.join(base_path, "bin/lap/aligner/sum_prob.py"), "-i", output_probs_location]
    out_cmd(call_arr)
    call(call_arr, stdout=open(output_sum_probs_location,'w'))
    results(output_sum_probs_location)


def run_samtools(options, sam_output_location, with_pileup = True, index=False):
    """ Takes a sam file and runs samtools to create bam, sorted bam, and mpileup. """

    bam_dir = options.output_dir + "/bam/"
    ensure_dir(bam_dir)
    bam_location = bam_dir + "library.bam"
    sorted_bam_location = bam_dir + "sorted_library"
    bam_fp = open(bam_location, 'w+')
    error_file_location = bam_dir + "error.log"
    error_fp = open(error_file_location, 'w+')

    #warning("About to run samtools view to create bam")
    call_arr = [os.path.join(base_path, "bin/Reapr_1.0.17/src/samtools"), "view", "-bS", sam_output_location]
    out_cmd(call_arr)
    #warning("That command outputs to file: ", bam_location)
    call(call_arr, stdout = bam_fp, stderr = error_fp)

    #warning("About to attempt to sort bam")
    call_arr = [os.path.join(base_path, "bin/Reapr_1.0.17/src/samtools"), "sort", bam_location, sorted_bam_location]
    out_cmd(call_arr)
    call(call_arr, stderr = FNULL)

    coverage_file_dir = options.output_dir + "/coverage/"
    ensure_dir(coverage_file_dir)
    pileup_file = coverage_file_dir + "mpileup_output.out"
    p_fp = open(pileup_file, 'w')

    if with_pileup:
        call_arr = [os.path.join(base_path, "bin/Reapr_1.0.17/src/samtools"), "mpileup", "-A", "-f", options.fasta_file, sorted_bam_location + ".bam"]
        out_cmd(call_arr)
        results(pileup_file)
        #warning("That command outputs to file: ", pileup_file)
        call(call_arr, stdout = p_fp, stderr = FNULL)

    if index:
        call_arr = [os.path.join(base_path, "bin/Reapr_1.0.17/src/samtools"), "index", sorted_bam_location + ".bam"]
        out_cmd(call_arr)
        call(call_arr, stdout = FNULL, stderr = FNULL)

    return (bam_location, sorted_bam_location, pileup_file)


def run_depth_of_coverage(options, pileup_file):
    """ Run depth of coverage. """

    dp_fp = options.output_dir + "/coverage/errors_cov.gff"
    abundance_file = options.coverage_file
    #call_arr = ["src/py/depth_of_coverage.py", "-a", abundance_file, "-m", pileup_file, "-w", options.window_size, "-o", dp_fp, "-g", "-e"]
    call_arr = [os.path.join(base_path, "src/py/depth_of_coverage.py"), "-m", pileup_file, "-w", options.window_size, "-o", dp_fp, "-g", "-e"]
    out_cmd(call_arr)
    call(call_arr)
    results(dp_fp)

    return dp_fp    


def run_reapr(options, sorted_bam_location):
    """ Run REAPR. """

    reapr_command = os.path.join(base_path, "bin/Reapr_1.0.17/reapr")

    #warning("About to run facheck")
    call_arr = [reapr_command, "facheck", options.fasta_file ]
    out_cmd(call_arr)
    call(call_arr)

    reapr_output_dir = options.output_dir + "/reapr"
    reapr_perfect_prefix = options.output_dir + "/r_perfect_prefix"
    
    #warning("About to run reapr pipeline")
    call_arr = [reapr_command, "pipeline", options.fasta_file,\
            sorted_bam_location + ".bam", reapr_output_dir]
    out_cmd(call_arr)
    call(call_arr, stdout=FNULL, stderr=FNULL)

    call_arr = ["gunzip", reapr_output_dir + "/03.score.errors.gff"]
    out_cmd(call_arr)
    call(call_arr, stdout=FNULL, stderr=FNULL)

    if os.path.exists(reapr_output_dir + "/03.score.errors.gff"):
        return reapr_output_dir + "/03.score.errors.gff"
    else:
        return None



if __name__ == '__main__':
    main()
