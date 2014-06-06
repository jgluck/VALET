#! /bin/bash
set -o verbose

cd bin/

# Install LAP
unzip lap_release_1.1.zip
mv lap_release_1.1 lap

# Install REAPR
tar -xzf Reapr_1.0.17.tar.gz
cd Reapr_1.0.17
./install.sh
cd .. 

# Install bowtie2
wget http://cbcb.umd.edu/~cmhill/files/bowtie2-2.2.2-linux-x86_64.zip
unzip bowtie2-2.2.2-linux-x86_64.zip
cd ..