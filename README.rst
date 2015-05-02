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
  
2) Download from github, on the write site of https://github.com/dieterich-lab/DCC, click Download Zip. Unzip the file, do the
same as 1.

=====
Usage
=====
1) Map RNA-seq data with STAR (Dobin et al., 2013)

.. code-block:: bash

STAR --runThreadN 10   --genomeDir /data/Indices/STAR/Rattus_norvegicus/rn6  --outSAMtype BAM Unsorted  --genomeLoad LoadAndKeep   --readFilesIn /data/projects/Rat/${GSM}/${GSM}.fastq.gz   --readFilesCommand zcat   --outFileNamePrefix ${GSM} --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15
