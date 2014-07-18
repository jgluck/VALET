#Valet

Pipeline for evaulating metagenomic assemblies.

Once the repository has been cloned, to install the required tools run the command:
```
./setup.sh
```

Test the installation by running the following command:

```
src/py/pipeline.py -a test/carsonella_asm.fna -c test/carsonella_asm.cvg -q -1 test/lib1.1.fastq -2 test/lib1.2.fastq -o test_validate
```

The summary of misassemblies are stored in [GFF format](http://www.sanger.ac.uk/resources/software/gff/spec.html):

```
contig00001     DEPTH_COV       Low_coverage    1       100     27.910000       .       .       low=35.007000;high=65.013000
contig00001     DEPTH_COV       High_coverage   23501   23600   68.330000       .       .       low=35.007000;high=65.013000
contig00001     DEPTH_COV       High_coverage   52301   52400   66.290000       .       .       low=35.007000;high=65.013000
contig00001     DEPTH_COV       High_coverage   57801   57900   65.610000       .       .       low=35.007000;high=65.013000
contig00001     DEPTH_COV       Low_coverage    77901   78000   34.440000       .       .       low=35.007000;high=65.013000
contig00001     DEPTH_COV       Low_coverage    78001   78100   34.550000       .       .       low=35.007000;high=65.013000
contig00001     REAPR   Read_orientation        88920   97033   .       .       .       Note=Warning:   Bad     read    orientation;colour=1
contig00001     REAPR   FCD     89074   90927   0.546142        .       .       Note=Error:     FCD     failure;colour=17
contig00001     REAPR   Frag_cov        90847   95221   0.000228571     .       .       Note=Error:     Fragment        coverage        too     low;color=15
```
