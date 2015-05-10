*****************************************
DCC: detect circRNAs from chimeric reads
*****************************************
DCC is a python package intended to detect and quantify circRNAs with high specificity. DCC works with the STAR (Dobin et al., 2013) chimeric.out.junction 
files which contains chimerically aligned reads including circRNA junction spanning reads. DCC detect and quantify circRNA junction 
spanning reads with high reliability. 

=============
Installation
=============

DCC depences on pysam, pybedtools, numpy.
The installation process of DCC with automatically check for the dependencies, if any dependence missing from path, it will automatically
install.

1) From git clone

.. code-block:: bash

  $ git clone git@github.com:dieterich-lab/DCC.git
  
  $ cd DCC
  
  $ python setup.py install
  
2) Download from github, on the write site of https://github.com/dieterich-lab/DCC, click 'Download Zip'. Unzip the file, do the same as 1.


========
Usage
========
The detection of circRNAs from RNAseq data through DCC can be summarised by three steps:

1. Map RNAseq data from quality checked fastq files. For paired end data, it is recommended to map the with two pairs jointly, and **also** separately. 

2. Prepare files needed by DCC. In summary, only one file is mandatory: 'junction_files'. Three other files are recommended: 'Repetitive_regions.gtf', 'mate1' and 'mate2'.

3. Run DCC. 

========================
Step by step tutorial
========================

1) Map RNA-seq data with STAR (Dobin et al., 2013)

.. code-block:: bash

  $ STAR --runThreadN 10   --genomeDir <genome directory>  --outSAMtype BAM Unsorted  --genomeLoad LoadAndKeep   --readFilesIn SamplePairedRead_1.fastq.gz  SamplePairedRead_2.fastq.gz   --readFilesCommand zcat   --outFileNamePrefix <sample prefix> --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15

2) Detect circRNAs from chimeric.out.junction file

.. code-block:: bash

  $ DCC

If DCC is not in you path
