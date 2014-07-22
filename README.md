#VALET

Pipeline for evaulating metagenomic assemblies.

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

```

The summary of misassemblies are stored in [GFF format](http://www.sanger.ac.uk/resources/software/gff/spec.html):

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
    -1 $HMP_SAMPLE/SRS014465.denovo_duplicates_marked.trimmed.1.fastq \
    -2 $HMP_SAMPLE/SRS014465.denovo_duplicates_marked.trimmed.2.fastq \
    -o SRS014465_valet \
    --window-size 100 --min-coverage 10 --coverage-multiplier 0.0 --threads 32 \
    --ignore-ends 100 --min-contig-length 1000
```

Lets take a closer look at the parameters we've selected:
* **-window-size 100** sets the sliding window size to be 100bp.
* **--min-coverage 10** sets the minimum contig coverage to 10x.  Any contig with less than a median contig coverage of 10x will not be flagged for misassemblies.
* **--coverage-multiplier 0.0** contigs are binned by exponentially increasing the bin size. At each new iteration, low = previous high, and high = high * coverage_multiplier. If coverage-multiplier is 0, each contig is binned with contigs sharing the same whole number.
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



