#VALET
Pipeline for evaulating metagenomic assemblies.

## Installing VALET
Once the repository has been cloned, to install the required tools run the command:
```
git clone https://github.com/jgluck/VALET.git
cd VALET/
./setup.sh

...

# Let your shell know where to find the VALET pipeline.
export VALET=`pwd`/src/py/
```

Test the installation by running the following command:

```
$VALET/pipeline.py -a test/carsonella_asm.fna -c test/carsonella_asm.cvg -q -1 test/lib1.1.fastq -2 test/lib1.2.fastq -o test_validate
```
```
---------------------------------------------------------------------------
STEP:    FILTERING ASSEMBLY CONTIGS LESS THAN 1000 BPs
RESULTS:         test_validate/filtered_assembly.fasta
---------------------------------------------------------------------------
STEP:    ALIGNING READS

...
...

---------------------------------------------------------------------------
STEP:    SUMMARY
RESULTS:         test_validate/summary.gff
RESULTS:         test_validate/suspicious.gff
RESULTS:         test_validate/summary.tsv

```

The flagged regions (potential misassemblies) are stored in two files **[OUTPUT_DIR]/summary.gff** and **[OUTPUT_DIR]/suspicious.gff**.
The flagged regions are stored in [GFF format](http://www.sanger.ac.uk/resources/software/gff/spec.html).  If multiple misassembly signatures overlap, their intersection is written to **suspicious.gff**.

```
contig00001     REAPR   Read_orientation        88920   97033   .       .       .       Note=Warning: Bad read orientation;colour=1
contig00001     REAPR   FCD     89074   90927   0.546142        .       .       Note=Error: FCD failure;colour=17
contig00001     DEPTH_COV       Low_coverage    90818   95238   29.500000       .       .       low=30.000000;high=70.000000;color=#7800ef
contig00001     REAPR   Frag_cov        90847   95221   0.000228571     .       .       Note=Error: Fragment coverage too low;color=15
contig00001     REAPR   FCD     95142   96991   0.534348        .       .       Note=Error: FCD failure;colour=17
contig00001     REAPR   Read_orientation        132895  134216  .       .       .       Note=Warning: Bad read orientation;colour=1
contig00001     REAPR   FCD     132942  136743  0.519051        .       .       Note=Error: FCD failure;colour=17
contig00001     REAPR   Read_orientation        134218  135390  .       .       .       Note=Warning: Bad read orientation;colour=1
contig00001     REAPR   Frag_cov        134757  134882  0.00793651      .       .       Note=Error: Fragment coverage too low;color=15
contig00001     REAPR   Read_orientation        135392  136654  .       .       .       Note=Warning: Bad read orientation;colour=1
```

In addition, a breakdown of each contig's number of misassemblies is available in the **[OUTPUT_DIR]/summary.tsv** file:

```
contig_name     contig_length   low_cov low_cov_bps     high_cov        high_cov_bps    reapr   reapr_bps       breakpoints     breakpoints_bps
contig00001     160502  1       4421    0       0       9       23879   3       103
```

## Example usages

* I have a pair-end FASTQ library.  I want to ignore misassemblies found in contigs smaller than 1000bp, and contigs with less than 10x coverage.

```
$VALET/pipeline.py -a [ASSEMBLY_FASTA] -q -1 [FIRST_MATES] -2 [SECOND_MATES] \
    --min-contig-length 1000 \
    --min-coverage 10 \
    -o [OUTPUT_DIR]
```


## Tutorial: Finding misassemblies in the Human Microbiome Project

Here we show how VALET can be used to find misassemblies in the Human Microbiome Project (http://www.hmpdacc.org/).

Before you continue, make sure the following tools are installed:
* **VALET**
* **git**
* GNU tools:
    * **wget** (can also use **curl**)
    * **tar**
* **Integrative Genomics Viewer** (http://www.broadinstitute.org/igv/home).  Any other genome browser that accepts BAM/SAM alignments and can overlay GFF files are acceptable.

### Installing VALET
```
git clone https://github.com/jgluck/VALET.git
cd VALET/
./setup.sh

...

# Let your shell know where to find the VALET pipeline.
export VALET=`pwd`/src/py/
```

### Downloading HMP sample SRS014465
```
mkdir samples
cd samples

# Download reads
wget ftp://public-ftp.hmpdacc.org/Illumina/vaginal_introitus/SRS014465.tar.bz2
tar xvjf SRS014465.tar.bz2
cd ..

# Download assembly
wget ftp://public-ftp.hmpdacc.org/HMASM/PGAs/vaginal_introitus/SRS014465.scaffolds.fa.bz2
tar xvjf SRS014465.scaffolds.fa.bz2

# Export sample directory to a path variable
export HMP_SAMPLE=`pwd`/SRS014465/

cd ..
```

### Running VALET

```
$VALET/pipeline.py -a $HMP_SAMPLE/SRS014465.scaffolds.fa \
    -q -1 $HMP_SAMPLE/SRS014465.denovo_duplicates_marked.trimmed.1.fastq \
    -2 $HMP_SAMPLE/SRS014465.denovo_duplicates_marked.trimmed.2.fastq \
    -o SRS014465_valet \
    --window-size 100 --min-coverage 10 --threads 32 \
    --ignore-ends 100 --min-contig-length 1000
```

Lets take a closer look at the parameters we've selected:
* **-window-size 100** sets the sliding window size to be 100bp.
* **--min-coverage 10** sets the minimum contig coverage to 10x.  Any contig with less than a median contig coverage of 10x will not be flagged for misassemblies.
* **--threads 32** sets the threads to 32. 
* **--ignore-ends 100** any flagged region within 100bp of the ends of a contig will be ignored.
* **--min-contig-length 1000** will have VALET ignore any flagged regions on contigs smaller than 1000bp.

```
INFO:    Coverage file not provided, will create one.
---------------------------------------------------------------------------
STEP:    FILTERING ASSEMBLY CONTIGS LESS THAN 1000 BPs
RESULTS:         SRS014465_valet/filtered_assembly.fasta
---------------------------------------------------------------------------
STEP:    ALIGNING READS
...
```

### Investigating potential misassemblies using IGV
Any genomics viewer that supports FASTA, BAM, and GFF formats should be able to display the flagged misassemblies.  For this example, we will use Broad's IGV(http://www.broadinstitute.org/software/igv/download#binary).

```
# Download and install IGV
wget http://www.broadinstitute.org/igv/projects/downloads/IGV_2.3.34.zip
unzip IGV_2.3.34.zip
export IGV_PATH=`pwd`/IGV_2.3.34
```

Next, create a file containing the IGV batch parameters below (*this will be automated in a future releases!*):

```
echo "new
genome filtered_assembly.fasta
load summary.gff
load suspicious.gff
load bam/sorted_library.bam" > SRS014465_valet/IGV.batch
```

In order to view the assembly and flagged regions, change the directory to the SRS014465 VALET directory and run:
```
cd SRS014465_valet/
$IGV_PATH/igv.sh -b IGV.batch
```

Now you are free to explore the flagged regions!

## Options
```
$VALET/pipeline.py -h

Options:
  -h, --help            show this help message and exit
  -a FILE, --assembly-fasta=FILE
                        Candidate assembly file
  -r FILE, --reads=FILE
                        First Read File
  -1 FIRST_MATES, --1=FIRST_MATES
                        Fastq filenames separated by commas that contain the
                        first mates.
  -2 SECOND_MATES, --2=SECOND_MATES
                        Fastq filenames separated by commas that contain the
                        second mates.
  -c COVERAGE_FILE, --coverage-file=COVERAGE_FILE
                        Assembly created per-contig coverage file
  -o OUTPUT_DIR, --output-dir=OUTPUT_DIR
                        Output directory
  -w WINDOW_SIZE, --window-size=WINDOW_SIZE
                        Sliding window size when determining misassemblies.
  -q, --fastq           if set, input reads are fastq format (fasta by
                        default).
  -p THREADS, --threads=THREADS
                        Number of threads
  -I MIN_INSERT_SIZE, --minins=MIN_INSERT_SIZE
                        Min insert sizes for mate pairs separated by commas.
  -X MAX_INSERT_SIZE, --maxins=MAX_INSERT_SIZE
                        Max insert sizes for mate pairs separated by commas.
  -n ORIENTATION, --orientation=ORIENTATION
                        Orientation of the mates.
  -m MU, --mu=MU        average mate pair insert sizes.
  -t SIGMA, --sigma=SIGMA
                        standard deviation of mate pair insert sizes.
  -x MAX_ALIGNMENTS, --max-alignments=MAX_ALIGNMENTS
                        bowtie2 parameter to set the max number of alignments.
  -e EMAIL, --email=EMAIL
                        Email to notify when job completes
  -g MIN_COVERAGE, --min-coverage=MIN_COVERAGE
                        Minimum average coverage to run misassembly detection.
  -l COVERAGE_MULTIPLIER, --coverage-multiplier=COVERAGE_MULTIPLIER
                        When binning by coverage, the new high = high + high *
                        multiplier
  -s MIN_SUSPICIOUS_REGIONS, --min-suspicious=MIN_SUSPICIOUS_REGIONS
                        Minimum number of overlapping flagged miassemblies to
                        mark region as suspicious.
  -z MIN_CONTIG_LENGTH, --min-contig-length=MIN_CONTIG_LENGTH
                        Ignore contigs smaller than this length.
  -b IGNORE_END_DISTANCES, --ignore-ends=IGNORE_END_DISTANCES
                        Ignore flagged regions within b bps from the ends of
                        the contigs.
```


