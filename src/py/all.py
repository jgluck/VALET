#!/usr/bin/env python
from subprocess import call
import os
import time

def main():
    start_time = time.time()
    fastaFile = "data/input/pairs/soap.contig"
    outputPrefixDir = "data/bowtieIndex/"
    outputPrefix = outputPrefixDir + "soap"
    
    reads_untrimmed_location = ["data/input/pairs/SRS011061.denovo_duplicates_marked.trimmed.1.fastq","data/input/pairs/SRS011061.denovo_duplicates_marked.trimmed.2.fastq"]
    reads_trimmed_location = []
    
    samOutputLocationDir = "data/output/sam/"
    samOutputLocation = samOutputLocationDir + "soap.sam"
    
    singletonOutputDir = "data/output/singleton/"
    singletonOutputLocation = singletonOutputDir + "singletons.csv"

    ensure_dir(outputPrefixDir)
    ensure_dir(singletonOutputDir)

    #call(["/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/bowtie2-build", fastaFile, outputPrefix]) 


    ensure_dir(samOutputLocationDir)

    should_align = True
    min_read_length = 300


    if len(reads_untrimmed_location)>1 and should_align:
        for readFile in reads_untrimmed_location:
            call(["fastx_clipper","-Q33","-l %d" %(min_read_length), "-i%s" %(readFile), "-o%s" %(readFile+"trimmed.fastq")])
            reads_trimmed_location.append(readFile+"trimmed.fastq")
        call(["/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/bowtie2", "-x", outputPrefix, "-1", reads_trimmed_location[0], "-2", reads_trimmed_location[1], "-S", samOutputLocation, "--un", singletonOutputLocation, "-q", "-I 100", "-X 500", "-p 8"])
    elif should_align:
        call(["fastx_clipper","-Q33","-l %d" %(min_read_length), "-i%s" %(reads_untrimmed_location[0]), "-o%s" %(reads_untrimmed_location[0]+"trimmed.fastq")])
        reads_trimmed_location.append(readFile+"trimmed.fastq")
        call(["/cbcb/project-scratch/cmhill/tools/bowtie2-2.2.2/bowtie2", "-x", outputPrefix, "-U", reads_trimmed_location[0], "-S",  samOutputLocation, "--un", singletonOutputLocation, "-p 8"])

    outputProbsDir = "data/lap/"
    ensure_dir(outputProbsDir)
    outputProbsLocation = outputProbsDir + "output.prob"
    fp = open(outputProbsLocation, "w")

    ##lap
    #call(["bin/lap/aligner/calc_prob.py", "-a", fastaFile,  "-s", samOutputLocation,  "-q", "-i",  reads_trimmed_Location[0]], stdout=fp)
    #call(["bin/lap/aligner/sum_prob.py", "-i", outputProbsLocation])


    outputBreakpointDir = "data/output/breakpoint/"
    ouputBreakpointLocation = outputBreakpointDir + "errorsDetected.csv"
    ensure_dir(outputBreakpointDir)

    ##breakpoint --- takes too long with current implementation
    ##alpha is length that must match
    #call(["bin/assembly-testing/breakpoint-detection/breakpoint_indices.py", "-a", fastaFile, "-u", singletonOutputLocation, "-o", ouputBreakpointLocation ,  "--alpha", "20","--algorithm","naive"])


    reapr_command = "bin/Reapr_1.0.17/reapr"

    call([reapr_command, "facheck", "data/input/pairs/soap.contig"])

    bamDir = "data/output/bam/"
    ensure_dir(bamDir)
    bamLocation = bamDir + "soap.bam"
    sorted_bam_location = bamDir+"sorted_soap"
    bam_fp = open(bamLocation, "rw")

    call(["samtools", "view", "-bS", samOutputLocation],stdout=bam_fp)

    call(["samtools", "sort", bamLocation,sorted_bam_location])

    reapr_output_dir = "data/output/reapr"
    reapr_perfect_prefix = "data/output/r_perfect_prefix"
    ##does not want this
    #ensure_dir(reapr_output_dir)
    call([reapr_command,"perfectmap",fastaFile,reads_trimmed_location[0],reads_trimmed_location[1], "100", reapr_perfect_prefix])
    call([reapr_command,"pipeline",fastaFile,sorted_bam_location+".bam",reapr_output_dir, reapr_perfect_prefix])

    notify_complete("jonathangluck08854@gmail.com",time.time()-start_time)


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def cleanup(cleanOpts):
    for (field,val) in cleanOpts.iterItems():
        if field == "sam":
            #remove sam
            print ""
        elif field == "fastx":
            print ""


def notify_complete(target_email,t):
    call(['echo "Completed in %d" | mail -s "Job Completed" %s' % (t, target_email) ],shell=True)
if __name__ == '__main__':
    main()
