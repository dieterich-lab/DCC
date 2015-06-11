*****************************************
DCC: detect circRNAs from chimeric reads
*****************************************
DCC is a python package intended to detect and quantify circRNAs with high specificity. DCC works with the STAR (Dobin et al., 2013) chimeric.out.junction 
files which contains chimerically aligned reads including circRNA junction spanning reads. 

=============
Installation
=============

DCC dependes on pysam, pybedtools, numpy and HTSeq
The installation process of DCC will automatically check for the dependencies, if any dependence missing from path, it will be automatically
installed.

1) From git clone

.. code-block:: bash

  $ git clone git@github.com:dieterich-lab/DCC.git
  
  $ cd DCC
  
  $ python setup.py install
  
2) Download from github. 

On the right site of https://github.com/dieterich-lab/DCC, click 'Download Zip'. Unzip the file, do the same as 1.

3) To check for installation, do

.. code-block:: bash
  
  # Get help
  $ DCC -h
  # If somehow DCC is not included in your path, you also call DCC by:
  $ python <DCC directory>/scripts/DCC <Options>

========
Usage
========
The detection of circRNAs from RNAseq data through DCC can be summarised by three steps:

1. Map RNAseq data from quality checked fastq files. For paired end data, it is recommended to map the with two pairs jointly, and **also** separately. This is because STAR does not output reads or read pairs which have more more than one chimeric junctions. 

2. Prepare files needed by DCC. In summary, only one file is mandatory: 'samplesheet', which specifies where your Chimeric.out.junction files are stored, with one sample per line. Three other files are recommended: 1). 'Repetitive_regions.gtf', a GTF format annotation of repetitive regions, which is used to filter out circRNA candidates from repetitive regions. 2). For paired end sequencing, 'mate1' and 'mate2', which are Chimeric.out.junction files from mate separate mapping.

3. Run DCC. DCC can be run for different perpers with different models. In summary, 1) Run DCC to detect circRNAs and host gene expression (use -D and -G option ) 2). Run DCC only to detec circRNAs (use -D option only). 3) Run DCC to count host gene expression with a custom provided circRNA list in BED format (use -G option, provide custom circRNA with -C option). 

========================
Step by step tutorial
========================
In this tutorial, we use Westholm et al. data as an example. The data are paired end, stranded RiboMinus RNAseq data from Drosophila.melanogaster, consisting of samples 3 developmental stages (1days, 4days and 20days) from heads. You can download the data as fastq files with NCBI SRA accession number: SRP001696. 

1. Map RNA-seq data with STAR (Dobin et al., 2013). Note that --alignSJoverhangMin and --chimJunctionOverhangMin should use the same value, to make the circRNA expression and linear gene expression level comparable. 

* Do pairs joined mapping first if your data is paired end, and then do mates separate mapping. If the data is single end, only one mapping step is needed.

.. code-block:: bash

  $ mkdir Sample1
  $ cd Sample1
  $ STAR --runThreadN 10   --genomeDir [genome]  --outSAMtype BAM Unsorted --readFilesIn Sample1_1.fastq.gz  Sample1_2.fastq.gz   --readFilesCommand zcat  --outFileNamePrefix [sample prefix] --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15


* (Skip when you have single end data). Mates separate mapping. Be careful that, what you define as first mate (mate1) should also appears the first in the joined mapping. In this case, SamplePairedRead_1.fastq.gz is the first mate which come first.

.. code-block:: bash

  # Make a directory for mate1
  $ mkdir mate1
  $ STAR --runThreadN 10   --genomeDir [genome]  --outSAMtype None --readFilesIn Sample1_1.fastq.gz  --readFilesCommand zcat   --outFileNamePrefix [sample prefix] --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15

  $ cd ..
  $ mkdir mate2
  # Do the same mapping as mate1 for mate2

2. Detect circRNAs from chimeric.out.junction file

* It is strongly recommended to specify a repetitive region file in GTF format for filtering. You can obtain this file through UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables. Select your genome, select group as "Repeats" or "Variation and Repeats". For the track, I recommend chose all possible repeats and combine the results. **NOTE**: the output file needs to comply with GTF format specification.

* Prepare path files to specify where is your chimeric.junction.out files are. 

  First, "samplesheet" file, in which you specify your chimeric.out.junction file's absolute paths (mates joined mapping chimeric.out.junction files, for paired end data), one line per sample. 

  Second (only if you have paired end sequencing data), "mate1" and "mate2" files. As with the "samplesheet" file, you specify where your mate1 and mate2 separately mapped chimeric.junction.out files are.

  You can find a example of this files for Westholm et al. data at:
  
.. code-block:: bash

  $ <DCC directory>/DCC/data/chimeric_junctions # Mates jointly mapped chimeric.junction.out files
  $ <DCC directory>/DCC/data/mate1 # Mate1 independently mapped chimeric.junction.out files
  $ <DCC directory>/DCC/data/mate1 # Mate2 independently mapped chimeric.junction.out files

* After all the preparation steps, you can now run DCC for circRNA detection. 


  # Call DCC to detect circRNAs. 

  $ DCC @samplesheet -mt1 @mate1 -mt2 @mate2 -D -S -R [Repeats].gtf -an [Annotation].gtf -Pi -F -M -Nr 10 5 20 6 -fg -G -A [Reference].fa

NOTE: -F flag is mandatory, if you want to filter on the results. All filtering steps are not mandatory, but strongly recommended.

**Finished!!!**


The output of DCC include: CircRNACount, CircCoordinates, LinearCount (this three in DCC working directories) and [header]Circ_Skip_Count (in the sample directories).


**CircRNACount:** a table containing read counts for circRNAs detected. First thre columns are chr, circRNA start, circRNA end. From fourth column one are the circRNA read counts, one sample per column, shown in the order given in your samplesheet.
**CircCoordinates:** CircRNA annotation in BED format. The columns are chr, start, end, genename, junctiontype (come from STAR, 1 for AG/GT, 2 for TC/CA, 0 for non-canonical junction), strand, circRNA region (startregion-endregion), overall regions (the genomic features circRNA coordinates interval covers).

**LinearCount:** host gene expression count table, same setup with CircRNACount file.

**[header]Circ_Skip_Count:** CircSkip junctions.

If you only want to detect circRNA without counting host gene expression, you can do

.. code-block:: bash

  $ DCC @samplesheet -mt1 @mate1 -mt2 @mate2 -D -S -R [Repeats].gtf -an [Annotation].gtf -Pi -F -M -Nr 10 5 20 6 -fg

If you have your own list of circRNAs in BED format, you can cant host gene expression for your list of circRNAs using DCC by:

.. code-block:: bash

  $ DCC @samplesheet -C [your list] -G -A [Reference].fa


========================================================================
Test for host-independently regulated circRNAs with CircTest package
========================================================================

1) Install CircTest package as described: https://github.com/dieterich-lab/CircTest

2) Read and load DCC output into R

.. code-block:: R

  CircRNACount <- read.delim('CircRNACount',header=T)
  LinearCount <- read.delim('LinearCount',header=T)
  CircCoordinates <- read.delim('CircCoordinates',header=T)

  CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount, Nreplicates = 6, filter.sample = 6, filter.count = 5, percentage = 0.1)
  CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]
  LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]

Alternatively, load the processed Westholm et al. data from CircTest package.

.. code-block:: R

  CircRNACount_filtered <- data(CircRNACount)
  CircCoordinates_filtered <- data(CircCoordinates)
  LinearCount_filtered <- data(LinearCount)

3) Test for host-independently regulated circRNAs

.. code-block:: R 

 test=Circ.test(CircRNACount_filtered,LinearCount_filtered,CircCoordinates_filtered,group=c(rep(1,6),rep(2,6),rep(3,6)))
 # Significant result show in a summary table
 View(test$summary_table)

4) Visuallize the significantly host-independently regulated circRNAs

.. code-block:: R

 for (i in rownames(test$summary_table))  {
  Circ.ratioplot( CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered, plotrow=i, 
                  groupindicator1=c(rep('1days',6),rep('4days',6),rep('20days',6)), 
                  lab_legend='Ages' )
 }









