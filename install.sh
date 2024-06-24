#!/bin/bash

# fastp installation
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# diamond installation
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
chmod +x diamond
rm diamond-linux64.tar.gz

# prodigal installation
wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux
mv prodigal.linux prodigal
chmod +x prodigal

# themisto installation
wget https://github.com/algbio/themisto/releases/download/v3.2.2/themisto_linux-v3.2.2.tar.gz
tar xzf themisto_linux-v3.2.2.tar.gz
mv ./themisto_linux-v3.2.2/themisto ./
chmod +x themisto
rm themisto_linux-v3.2.2.tar.gz
rm -r ./themisto_linux-v3.2.2
