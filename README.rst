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

1) Map RNA-seq data with STAR (Dobin et al., 2013). What is important is --alignSJoverhangMin and --chimJunctionOverhangMin showed use the same criteria, to make the circRNA expression and linear gene expression level comparable. 

a. Do pairs joined mapping first if your data is paired end, and then do mates separate mapping. If the data is single end, only one mapping step is needed.

.. code-block:: bash

  $ mkdir Sample1
  $ cd Sample1
  $ STAR --runThreadN 10   --genomeDir <genome directory>  --outSAMtype BAM Unsorted  --genomeLoad LoadAndKeep   --readFilesIn Sample1_1.fastq.gz  SamplePairedRead_2.fastq.gz   --readFilesCommand zcat   --outFileNamePrefix <sample prefix> --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15


b. (Skip when you have single end data). Mates separate mapping. Be careful that, what you define as first mate (mate1) should also appears the first in the joined mapping. In this case, SamplePairedRead_1.fastq.gz is the first mate which come first.

.. code-block:: bash

  # Make a directory for mate1
  $ mkdir mate1
  $ STAR --runThreadN 10   --genomeDir <genome directory>  --outSAMtype None --genomeLoad LoadAndKeep   --readFilesIn Sample1_1.fastq.gz  --readFilesCommand zcat   --outFileNamePrefix <sample prefix> --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15


.. code-block:: bash

  # Go back to the sample directory
  $ cd ..
  # Make a directory for mate2
  $ mkdir mate2
  # Do mapping for mate2, you can chose not to output linear aligned SAM/BAM files to same space.
  $ STAR --runThreadN 10   --genomeDir <genome directory>  --outSAMtype None --genomeLoad LoadAndKeep   --readFilesIn Sample1_2.fastq.gz  --readFilesCommand zcat   --outFileNamePrefix <sample prefix> --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15


2) Detect circRNAs from chimeric.out.junction file

a. It is strongly recommended to prepare a repetitive region file for filtering. You can get this file through:

.. code-block:: bash

  $ DCC @samples -mt1 @mate1 -mt2 @mate2 -D -S -R /data/genomes/drosophila_melanogaster/BDGP5_75/BDGP5.75.RepetitiveRegions.gtf -an /data/genomes/drosophila_melanogaster/BDGP5_75/BDGP5.75.gtf -Pi -F -M -Nr 10 5 20 6 -fg -G -A /data/genomes/drosophila_melanogaster/BDGP5_75/BDGP5.75.dna.toplevel.fa

If DCC is not in you path, you can run DCC by:

.. code-block:: bash

  $ python <DCC directory>/scripts/DCC <options>


