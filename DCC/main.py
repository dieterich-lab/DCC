#! /usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import functools
import logging
import multiprocessing
import os
import random
import string
import sys
import time

import pysam

import CombineCounts as Cc
import circAnnotate as Ca
import circFilter as Ft
import findcircRNA as Fc
import genecount as Gc
from fix2chimera import Fix2Chimera


def main():
    version = "0.4.6"

    parser = argparse.ArgumentParser(prog="DCC", formatter_class=argparse.RawDescriptionHelpFormatter,
                                     fromfile_prefix_chars="@",
                                     description="Contact: tobias.jakobi@med.uni-heidelberg.de || s6juncheng@gmail.com")

    parser.add_argument("--version", action="version", version=version)
    parser.add_argument("Input", metavar="Input", nargs="+",
                        help="Input of the Chimeric.out.junction file from STAR. Alternatively, a sample sheet "
                             "specifying where your chimeric.out.junction files are, each sample per line, "
                             "provide with @ prefix (e.g. @samplesheet)")
    parser.add_argument("-k", "--keep-temp", dest="temp", action="store_true", default=False,
                        help="Temporary files will not be deleted [default: False]")
    parser.add_argument("-T", "--threads", dest="cpu_threads", type=int, default=2,
                        help="Number of CPU threads used for computation [default: 2]")
    parser.add_argument("-O", "--output", dest="out_dir", default="./",
                        help="DCC output directory [default: .]")
    parser.add_argument("-t", "--temp", dest="tmp_dir", default="_tmp_DCC/",
                        help="DCC temporary directory [default: _tmp_DCC/]")

    group = parser.add_argument_group("Find circRNA Options", "Options to find circRNAs from STAR output")
    group.add_argument("-D", "--detect", action="store_true", dest="detect", default=False,
                       help="Enable circRNA detection from Chimeric.out.junction files [default: False]")
    group.add_argument("-ss", action="store_true", dest="secondstrand", default=False,
                       help="Must be enabled for stranded libraries, aka 'fr-secondstrand' [default: False]")
    group.add_argument("-N", "--nonstrand", action="store_false", dest="strand", default=True,
                       help="The library is non-stranded [default stranded]")
    group.add_argument("-E", "--endTol", dest="endTol", type=int, default=5, choices=range(0, 10),
                       help="Maximum base pair tolerance of reads extending over junction sites [default: 5]")
    group.add_argument("-m", "--maximum", dest="max", type=int, default=1000000,
                       help="The maximum length of candidate circRNAs (including introns) [default: 1000000]")
    group.add_argument("-n", "--minimum", dest="min", type=int, default=30,
                       help="The minimum length of candidate circRNAs (including introns) [default 30]")
    group.add_argument("-an", "--annotation", dest="annotate",
                       help="Gene annotation file in GTF/GFF3 format, to annotate "
                            "circRNAs by their host gene name/identifier")

    group.add_argument("-Pi", "--PE-independent", action="store_true", dest="pairedendindependent", default=False,
                       help="Has to be specified if the paired end mates have also been mapped separately."
                            "If specified, -mt1 and -mt2 must also be provided [default: False]")
    group.add_argument("-mt1", "--mate1", dest="mate1", nargs="+",
                       help="For paired end data, Chimeric.out.junction files from mate1 independent mapping result")
    group.add_argument("-mt2", "--mate2", dest="mate2", nargs="+",
                       help="For paired end data, Chimeric.out.junction files from mate2 independent mapping result")
    parser.add_argument_group(group)

    group = parser.add_argument_group("Filtering Options", "Options to filter the circRNA candidates")
    group.add_argument("-F", "--filter", action="store_true", dest="filter", default=False,
                       help="If specified, the program will perform a recommended filter step on the detection results")
    group.add_argument("-f", "--filter-only", dest="filteronly", nargs=2,
                       help="If specified, the program will only filter based on two files provided: "
                            "1) a coordinates file [BED6 format] and 2) a count file. E.g.: -f example.bed counts.txt")
    group.add_argument("-M", "--chrM", action="store_true", dest="chrM", default=False,
                       help="If specified, circRNA candidates located on the mitochondrial chromosome will be removed")
    group.add_argument("-R", "--rep_file", dest="rep_file",
                       help="Custom repetitive region file in GTF format to filter out "
                            "circRNA candidates in repetitive regions")
    group.add_argument("-L", "--Ln", dest="length", type=int, default=50,
                       help="Minimum length in base pairs to check for repetitive regions [default 50]")
    group.add_argument("-Nr", nargs=2, type=int, metavar=("countthreshold", "replicatethreshold"), default=[2, 5],
                       help="countthreshold replicatethreshold [default: 2,5]")
    group.add_argument("-fg", "--filterbygene", action="store_true", dest="filterbygene", default=False,
                       help="If specified, filter also by gene annotation (candidates are not allowed to span"
                            " more than one gene) default: False")
    parser.add_argument_group(group)

    group = parser.add_argument_group("Host gene count Options", "Options to count host gene expression")
    group.add_argument("-G", "--gene", action="store_true", dest="gene", default=False,
                       help="If specified, the program will count host gene expression given circRNA coordinates "
                            "[default: False]")
    group.add_argument("-C", "--circ", dest="circ",
                       help="User specified circRNA coordinates, any tab delimited file with first three "
                            "columns as circRNA coordinates: chr\tstart\tend, which DCC will use to count "
                            "host gene expression")
    group.add_argument("-B", "--bam", dest="bam", nargs="+",
                       help="A file specifying the mapped BAM files from which host gene expression is computed; "
                            "must have the same order as input chimeric junction files")
    group.add_argument("-A", "--refseq", dest="refseq",
                       help="Reference sequence FASTA file")

    parser.add_argument_group(group)

    options = parser.parse_args()

    timestr = time.strftime("%Y-%m-%d_%H%M")

    if not os.path.isdir(options.out_dir):
        try:
            os.makedirs(options.out_dir)
        except OSError:
            print "Could not create output folder %s" % options.out_dir
            logging.info("Could not create output folder %s" % options.out_dir)

            exit(-1)
    else:
        print "Output folder %s already exists, reusing" % options.out_dir

    # create temporary directory if not existing

    if not os.path.isdir(options.tmp_dir):
        try:
            os.makedirs(options.tmp_dir)
        except OSError:
            print "Could not create temporary folder %s" % options.tmp_dir
            exit(-1)
    else:
        print "Temporary folder %s already exists, reusing" % options.tmp_dir

    logging.basicConfig(filename=os.path.join(options.out_dir, "DCC-" + timestr + ".log"),
                        filemode="w", level=logging.DEBUG,
                        format="%(asctime)s %(message)s")

    logging.info("DCC %s started" % version)
    print "DCC %s started" % version
    logging.info('DCC command line: ' + ' '.join(sys.argv))

    # Get input file names
    filenames = [os.path.basename(name) for name in options.Input]
    samplelist = "\t".join(filenames)

    # make sure the user supplied path variables have a trailing /
    options.tmp_dir = os.path.join(options.tmp_dir, "")
    options.out_dir = os.path.join(options.out_dir, "")

    # set output file names
    output_coordinates = options.out_dir + "CircCoordinates"
    output_circ_counts = options.out_dir + "CircRNACount"
    output_linear_counts = options.out_dir + "LinearCount"
    output_skips = options.out_dir + "CircSkipJunctions"

    circfiles = []  # A list for .circRNA file names

    # check whether the junction file names have duplicates
    same = False
    if len(set(filenames)) != len(options.Input):
        logging.info(
            "Input file names have duplicates, add number suffix in input order to output files for distinction")
        print ("Input file names have duplicates, add number suffix in input order to output files for distinction")
        same = True

    cpu_count = multiprocessing.cpu_count()

    if options.cpu_threads <= cpu_count:
        print "%s CPU cores available, using %s" % (cpu_count, options.cpu_threads)
    else:
        print "Only %s CPU cores available while %s requested, falling back to %s" % \
              (cpu_count, options.cpu_threads, cpu_count)
        options.cpu_threads = cpu_count

    pool = multiprocessing.Pool(processes=options.cpu_threads)

    if options.annotate is not None:
        checkfile(options.annotate, True)
    if options.rep_file is not None:
        checkfile(options.rep_file, True)

    # Make instance
    cm = Cc.Combine(options.tmp_dir)
    circann = Ca.CircAnnotate(tmp_dir=options.tmp_dir, strand=options.strand)

    if (not options.mate1 or not options.mate1) and options.pairedendindependent:
        logging.info('-Pi (Paired independent mode) specified but -mt1, -mt2, or both are missing.  '
                     'Will not use mate information.')
        print('-Pi (Paired independent mode) specified but -mt1, -mt2, or both are missing.  '
              'Will not use mate information.')

        options.pairedendindependent = False

    if checkjunctionfiles(options.Input, options.mate1, options.mate2, options.pairedendindependent):
        logging.info("circRNA detection skipped due to empty junction files")
        print("circRNA detection skipped due to empty junction files")

        options.detect = False

    if options.bam and options.mate1 and (len(options.bam) != len(options.mate1)):
        logging.info("BAM file list is shorter than mate list. Maybe you forgot the @ (@file.list)?")
        print("BAM file list is shorter than mate list. Maybe you forgot the @ (@file.list)?")
        exit(-1)

    if options.detect:
        logging.info("Starting to detect circRNAs")
        if options.strand:
            logging.info("Stranded data mode")
        else:
            logging.info("Non-stranded data, the strand of circRNAs guessed from the strand of host genes")
            print "WARNING: non-stranded data, the strand of circRNAs guessed from the strand of host genes"

        # Start de novo circular RNA detection model
        # Create instances
        f = Fc.Findcirc(endTol=options.endTol, maxL=options.max, minL=options.min)
        sort = Fc.Sort()

        if options.pairedendindependent:
            print "Please make sure that the read pairs have been mapped both, combined and on a per mate basis"
            logging.info("Please make sure that the read pairs have been mapped both, combined and on a per mate basis")

            # Fix2chimera problem by STAR
            print ("Collecting chimera information from mates-separate mapping")
            logging.info("Collecting chimera information from mates-separate mapping")
            Input = fixall(options.Input, options.mate1, options.mate2, options.out_dir, options.tmp_dir)
        else:
            Input = options.Input

        if options.strand:

            if options.pairedendindependent:
                circfiles = pool.map(
                    functools.partial(wrapfindcirc, tmp_dir=options.tmp_dir, endTol=options.endTol, maxL=options.max,
                                      minL=options.min, strand=True, pairdendindependent=True, same=same), Input)
            else:
                circfiles = pool.map(
                    functools.partial(wrapfindcirc, tmp_dir=options.tmp_dir, endTol=options.endTol, maxL=options.max,
                                      minL=options.min, strand=True, pairdendindependent=False, same=same), Input)

        else:
            if options.pairedendindependent:
                circfiles = pool.map(
                    functools.partial(wrapfindcirc, tmp_dir=options.tmp_dir, endTol=options.endTol, maxL=options.max,
                                      minL=options.min, strand=False, pairdendindependent=True, same=same), Input)
            else:
                circfiles = pool.map(
                    functools.partial(wrapfindcirc, tmp_dir=options.tmp_dir, endTol=options.endTol, maxL=options.max,
                                      minL=options.min, strand=False, pairdendindependent=False, same=same), Input)

        # Combine the individual count files
        # Create a list of ".circRNA" file names
        print("Combining individual circRNA read counts")
        logging.info("Combining individual circRNA read counts")

        cm.comb_coor(circfiles, strand=options.strand)
        cm.map(options.tmp_dir + "tmp_coordinates", circfiles, strand=options.strand)

        res = cm.combine([files + "mapped" for files in circfiles], col=7, circ=True)

        # swap strand if the sequences are sense strand
        if (options.secondstrand and options.strand):
            logging.info("Swapping strand information")
            strand_swap = {}
            strand_swap["+\n"] = "-\n"
            strand_swap["-\n"] = "+\n"
            toswap = open(options.tmp_dir + "tmp_coordinates").readlines()
            swaped = open(options.tmp_dir + "tmp_coordinatesswaped", "w")
            for lin in toswap:
                lin_split = lin.split("\t")
                lin_split[5] = strand_swap[lin_split[5]]
                swaped.write("\t".join(lin_split))
            swaped.close()
            os.remove(options.tmp_dir + "tmp_coordinates")
            os.rename(options.tmp_dir + "tmp_coordinatesswaped", options.tmp_dir + "tmp_coordinates")

        if options.filter:
            cm.writeouput(options.tmp_dir + "tmp_circCount", res)
            if options.annotate:
                logging.info("Write in annotation")
                logging.info("Select gene features in Annotation file")
                annotation_tree = circann.selectGeneGtf(options.annotate)
                circann.annotate(options.tmp_dir + "tmp_coordinates", annotation_tree,
                                 options.tmp_dir + "tmp_coordinates_annotated")
                os.remove(options.tmp_dir + "tmp_coordinates")
                os.rename(options.tmp_dir + "tmp_coordinates_annotated", options.tmp_dir + "tmp_coordinates")
        else:
            cm.writeouput(output_circ_counts, res, samplelist, header=True)
            if options.annotate:
                logging.info("Write in annotation")
                logging.info("Select gene features in Annotation file")
                annotation_tree = circann.selectGeneGtf(options.annotate)
                circann.annotate(options.tmp_dir + "tmp_coordinates", annotation_tree,
                                 options.tmp_dir + "tmp_coordinates_annotated")
                circann.annotateregions(options.tmp_dir + "tmp_coordinates_annotated", annotation_tree,
                                        output_coordinates)
            else:

                os.rename(options.tmp_dir + "tmp_coordinates", output_coordinates)

    # Filtering
    if options.filter:
        logging.info("Filtering started")

        filt = Ft.Circfilter(length=options.length,
                             countthreshold=options.Nr[0],
                             replicatethreshold=options.Nr[1],
                             tmp_dir=options.tmp_dir)

        if not options.detect and not options.gene and options.filteronly:
            try:
                file2filter = options.filteronly[0]
                coorfile = options.filteronly[1]
                logging.info("Using files %s and %s for filtering" % (options.filteronly[0], options.filteronly[1]))
                print "Using files %s and %s for filtering" % (options.filteronly[0], options.filteronly[1])

            except IndexError:
                logging.error("Program exit because input error. Please check the input. If only use the program "
                              "for filtering, a coordinate file in bed6 format and a count file is needed")

                sys.exit("Please check the input. If only use the program for filtering, a coordinate file in "
                         "bed6 format and a count file is needed")

        elif not options.detect:

            sys.exit("Filter mode for detected circRNAs enabled without detection module.\nCombine with -f or -D.")

        elif options.detect:

            file2filter = options.tmp_dir + "tmp_circCount"
            coorfile = options.tmp_dir + "tmp_coordinates"
            logging.info("Using files _tmp_DCC/tmp_circCount and _tmp_DCC/tmp_coordinates for filtering")
            print "Using files _tmp_DCC/tmp_circCount and _tmp_DCC/tmp_coordinates for filtering"

        if options.rep_file:
            rep_file = options.rep_file
        else:
            # from pkg_resources import resource_filename
            # rep_file = resource_filename("DCC", "data/DCC.Repeats")
            rep_file = None
        count, indx = filt.readcirc(file2filter, coorfile)
        logging.info("Filtering by read counts")
        count0, indx0 = filt.filtercount(count, indx)  # result of first filtering by read counts

        # filt.makeregion(indx0)
        # nonrep_left,nonrep_right = filt.nonrep_filter("_tmp_DCC/tmp_left","_tmp_DCC/tmp_right",rep_file)
        # filt.intersectLeftandRightRegions(nonrep_left,nonrep_right,indx0,count0)

        if rep_file is not None:
            logging.info("Filter by non repetitive region")
            filt.filter_nonrep(rep_file, indx0, count0)
        else:
            filt.dummy_filter(indx0, count0)
            
        if not options.chrM and not options.filterbygene:
            filt.sortOutput(options.tmp_dir + "tmp_unsortedWithChrM", output_circ_counts,
                            output_coordinates, samplelist)

        # Filter chrM, if no further filtering, return "CircRNACount" and "CircCoordinates",
        # else return "_tmp_DCC/tmp_unsortedNoChrM"
        if options.chrM:
            logging.info("Deleting circRNA candidates from mitochondrial chromosome")
            filt.removeChrM(options.tmp_dir + "tmp_unsortedWithChrM")
            if not options.filterbygene:
                filt.sortOutput(options.tmp_dir + "tmp_unsortedNoChrM", output_circ_counts,
                                output_coordinates, samplelist)
        else:
            os.rename(options.tmp_dir + "tmp_unsortedWithChrM",
                      options.tmp_dir + "tmp_unsortedNoChrM")
            # Note in this case "_tmp_DCC/tmp_unsortedNoChrM" actually has chrM

        # Filter by gene annotation, require one circRNA could not from more than one gene.
        # return final "CircRNACount" and "CircCoordinates"
        if options.filterbygene:
            if options.annotate:
                logging.info("Filtering by gene annotation. "
                             "CircRNA candidates from more than one genes are deleted")
                circann.filtbygene(options.tmp_dir + "tmp_unsortedNoChrM",
                                   options.tmp_dir + "tmp_unsortedfilterbygene")
                filt.sortOutput(options.tmp_dir + "tmp_unsortedfilterbygene", output_circ_counts,
                                output_coordinates, samplelist)
            else:
                logging.warning(
                    "To filter by gene annotation, a annotation file in GTF/GFF format needed, "
                    "skiped filter by gene annotation")
                filt.sortOutput(options.tmp_dir + "tmp_unsortedNoChrM", output_circ_counts,
                                output_coordinates, samplelist)

        # Add annotation of regions
        if options.annotate:
            circann.annotateregions(output_coordinates, annotation_tree,
                                    output_coordinates)

        logging.info("Filtering finished")

    if options.gene:
        # import the list of bamfile names as a file
        if not options.bam:
            print "No BAM files provided (-B) trying to automatically guess BAM file names"
            logging.info("No BAM files provided (-B) trying to automatically guess BAM file names")
            bamfiles = convertjunctionfile2bamfile(options.Input)
            if not bamfiles:
                print "Could not guess BAM file names, please provides them manually via -B"
                logging.info("Could not guess BAM file names, please provides them manually via -B")
        else:
            bamfiles = options.bam

        if not options.refseq:
            print "Please provide reference sequence, program will not count host gene expression"
            logging.warning("Please provide reference sequence, program will not count host gene expression")

        if options.refseq:
            linearfiles = []  # A list for .linear file names
            unsortedBAMS = checkBAMsorting(bamfiles)

            # check whether the number of bamfiles is equale to the number of chimeric.junction.out files
            if len(bamfiles) != len(options.Input):
                logging.error("The number of bam files does not match with chimeric junction files")
                sys.exit("The number of bam files does not match with chimeric junction files")

            elif len(unsortedBAMS) > 0:
                logging.error("The following BAM files seem to be not sorted by coordinate or are missing an index:")
                logging.error(', '.join(unsortedBAMS))
                print("The following BAM files seem to be not sorted by coordinate or are missing an index:")
                print(', '.join(unsortedBAMS))
                sys.exit("Error: not all BAM files are sorted by coordinate or are missing indices")
            else:
                # For each sample (each bamfile), do one host gene count, and then combine to a single table

                if options.circ:
                    linearfiles = pool.map(
                        functools.partial(wraphostgenecount, tmp_dir=options.tmp_dir, circ_coor=options.circ,
                                          ref=options.refseq,
                                          countlinearsplicedreads=False), bamfiles)
                else:
                    linearfiles = pool.map(
                        functools.partial(wraphostgenecount, tmp_dir=options.tmp_dir, circ_coor=output_circ_counts,
                                          ref=options.refseq,
                                          countlinearsplicedreads=False), bamfiles)

                logging.info("Finished linear gene expression counting, start to combine individual sample counts")

            # Combine all to a individual sample host gene count to a single table
            res = cm.combine(linearfiles, col=6, circ=False)
            cm.writeouput_linear(output_linear_counts, res, samplelist, header=True)
            logging.info("Finished combine individual linear gene expression counts")

            if not options.temp:
                deleted = cm.deletefile("", linearfiles)
                logdeleted(deleted)

    # CircSkip junction
    if options.annotate and options.detect and not options.circ:
        logging.info("Count CircSkip junctions")
        print("Count CircSkip junctions")
        SJ_out_tab = getSJ_out_tab(options.Input)
        CircSkipfiles = findCircSkipJunction(output_coordinates, options.tmp_dir,
                                             options.annotate, circfiles, SJ_out_tab,
                                             strand=options.strand, same=same)
        fin = open(output_coordinates, "r").readlines()[1:]
        with open(options.tmp_dir + "tmp_CircCoordinatesNoheader", "w") as fout:
            fout.writelines(fin)
        cm.map(options.tmp_dir + "tmp_CircCoordinatesNoheader", CircSkipfiles, strand=options.strand, col=4)
        CircSkipfilesmapped = [fname + "mapped" for fname in CircSkipfiles]
        res = cm.combine(CircSkipfilesmapped, col=9)
        cm.writeouput(output_skips, res, samplelist, header=True)
    else:
        logging.info("CircSkip junctions not counted, not running in detected mode (-D)")

    # Delete temporary files
    if not options.temp:
        deletion_regex = r"^tmp_\.*"
        deleted = cm.deletefile(options.tmp_dir, deletion_regex)
        logdeleted(deleted)
        try:
            os.rmdir(options.tmp_dir)
        except OSError:
            print "Could not delete temporary folder %s: not empty" % options.tmp_dir
            logging.info("Could not delete temporary folder %s: not empty" % options.tmp_dir)

        try:
            os.rmdir(options.out_dir)
        except OSError:
            print "Not deleting output folder %s: contains files" % options.out_dir
            logging.info("Not deleting output folder %s: contains files" % options.out_dir)

        print("Temporary files deleted")

    logging.info("DCC completed successfully")


def fixall(joinedfnames, mate1filenames, mate2filenames, out_dir, tmp_dir):
    # Fix all 2chimera in one read/read pair for all inputs
    # outputs as a list of fixed filenames, with .fixed end
    outputs = []
    fx = Fix2Chimera(tmp_dir)
    # check mate1 and mate2 input
    if len(mate1filenames) == len(mate2filenames) == len(joinedfnames):
        for i in range(len(joinedfnames)):
            extension = "." + id_generator()
            fx.fixchimerics(mate1filenames[i], mate2filenames[i], joinedfnames[i],
                            os.path.join(tmp_dir, os.path.basename("tmp_" + joinedfnames[i]) + extension))
            outputs.append(os.path.join(tmp_dir, os.path.basename("tmp_" + joinedfnames[i]) + extension))
    else:
        logging.error("The number of input mate1, mate2 and joined mapping files are different")
        sys.exit("The number of input mate1, mate2 and joined mapping files are different")

    return outputs


def checkfile(filename, previousstate):
    # check for file existence
    if not os.path.isfile(filename):
        sys.exit("ERROR: Required file " + str(filename) + " is missing, exiting")
    # check for file content
    elif os.stat(filename).st_size == 0:
        print ("WARNING: File " + str(filename) + " is empty!")
        return True
    return previousstate


def checkjunctionfiles(joinedfnames, mate1filenames, mate2filenames, pairedendindependent):
    # Check if the junctions files have actually any content
    # if no, skip circRNA detection (and return True)

    skipcirc = False

    mate1empty = False
    mate2empty = False
    joinedempty = False

    if pairedendindependent:

        # check input files
        if len(mate1filenames) == len(mate2filenames) == len(joinedfnames):

            for i in range(len(joinedfnames)):
                # check for mate 1 files
                mate1empty = checkfile(mate1filenames[i], mate1empty)

                # check for mate 2 files
                mate2empty = checkfile(mate2filenames[i], mate2empty)

                # check for combined files
                joinedempty = checkfile(joinedfnames[i], joinedempty)

            if mate1empty or mate2empty or joinedempty:
                skipcirc = True
                logging.warning('One of the input junctions files is empty.')
                print('One of the input junctions files is empty.')
        else:
            skipcirc = True

            logging.warning('Input file lists have different length (mate 1 %d, mate 2 %d, joined %d).' % (
            len(mate1filenames), len(mate2filenames), len(joinedfnames)))

            print('Input file lists have different length (mate 1 %d, mate 2 %d, joined %d).' % (
            len(mate1filenames), len(mate2filenames), len(joinedfnames)))

        if skipcirc:
            logging.warning('Junction files seem empty, skipping circRNA detection module.')
            print('Junction files seem empty, skipping circRNA detection module.')

        return skipcirc

    else:

        for i in range(len(joinedfnames)):
            # check for combined files
            joinedempty = checkfile(joinedfnames[i], joinedempty)

        if joinedempty:
            skipcirc = True

        if skipcirc:
            logging.warning('Junction files seem empty, skipping circRNA detection module.')
            print('Junction files seem empty, skipping circRNA detection module.')

        return skipcirc


def logdeleted(deleted):
    for itm in deleted:
        logging.info("Deleted temporary file " + itm)


def mergefiles(output, *fnames):
    import shutil
    destination = open(output, "wb")
    for fname in fnames:
        shutil.copyfileobj(open(fname, "rb"), destination)
    destination.close()


def convertjunctionfile2bamfile(junctionfilelist):
    # only works for STAR-like names: Aligned.noS.bam
    def getbamfname(junctionfname):
        import re
        import os
        # Get the stored directory
        dirt = "/".join((junctionfname.split("/")[:-1])) + "/"
        p = r".*Aligned\..*bam"
        bamfname = ""
        for fname in os.listdir(dirt):
            if re.match(p, fname):
                bamfname = dirt + re.findall(p, fname)[0]
        if bamfname:
            return bamfname

    bamfnames = []
    for fname in junctionfilelist:
        entry = getbamfname(fname)
        if entry:
            bamfnames.append(getbamfname(fname))
    return bamfnames


# CircSkip junctions
def findCircSkipJunction(CircCoordinates, tmp_dir, gtffile, circfiles, SJ_out_tab, strand=True, same=False):
    from Circ_nonCirc_Exon_Match import CircNonCircExon
    CircSkipfiles = []
    CCEM = CircNonCircExon(tmp_dir)
    # Modify gtf file
    if not os.path.isfile(tmp_dir + "tmp_" + os.path.basename(gtffile) + ".exon.sorted"):
        CCEM.select_exon(gtffile)
    if CCEM.modifyExon_id(tmp_dir + "tmp_" + os.path.basename(gtffile) + ".exon.sorted"):
        # Start and end coordinates
        start2end = CCEM.print_start_end_file(CircCoordinates)
        Iv2Custom_exon_id, Custom_exon_id2Iv, Custom_exon_id2Length = CCEM.readNonUniqgtf(
            tmp_dir + "tmp_" + os.path.basename(gtffile) + ".exon.sorted.modified")
        if strand:
            circStartExons = CCEM.intersectcirc(tmp_dir + "tmp_start.bed", tmp_dir + "tmp_" + os.path.basename(
                gtffile) + ".exon.sorted.modified", isStartBED=True)  # Circle start or end to corresponding exons
        else:
            circStartExons = CCEM.intersectcirc(tmp_dir + "tmp_start.bed",
                                                tmp_dir + "tmp_" + os.path.basename(gtffile) + ".exon.sorted.modified",
                                                strand=False, isStartBED=True)
        circStartAdjacentExons, circStartAdjacentExonsIv = CCEM.findcircAdjacent(circStartExons, Custom_exon_id2Iv,
                                                                                 Iv2Custom_exon_id, start=True)
        if strand:
            circEndExons = CCEM.intersectcirc(tmp_dir + "tmp_end.bed", tmp_dir + "tmp_" + os.path.basename(
                gtffile) + ".exon.sorted.modified", isStartBED=False)  # Circle start or end to corresponding exons
        else:
            circEndExons = CCEM.intersectcirc(tmp_dir + "tmp_end.bed",
                                              tmp_dir + "tmp_" + os.path.basename(gtffile) + ".exon.sorted.modified",
                                              strand=False, isStartBED=False)
        circEndAdjacentExons, circEndAdjacentExonsIv = CCEM.findcircAdjacent(circEndExons, Custom_exon_id2Iv,
                                                                             Iv2Custom_exon_id, start=False)
        exonskipjunctions = CCEM.exonskipjunction(circStartAdjacentExonsIv, circEndAdjacentExonsIv, start2end)
        for indx, fname in enumerate(SJ_out_tab):
            if same:
                path = tmp_dir + "tmp_" + os.path.basename(fname).replace("SJ.out.tab", str(indx))
            else:
                path = tmp_dir + "tmp_" + os.path.basename(fname).replace("SJ.out.tab", "")
            junctionReadCount = CCEM.readSJ_out_tab(fname)

            if len(junctionReadCount) == 0:
                logging.error("Do you have SJ.out.tab files in your sample folder? DCC cannot find it")
                logging.info("Cannot fine SJ.out.tab files, please check the path. circSkip will not be output")
                break
            else:
                skipJctCount = CCEM.getskipjunctionCount(exonskipjunctions, junctionReadCount)
                circCount = CCEM.readcircCount(circfiles[indx])
                CircSkipfile = CCEM.printCirc_Skip_Count(circCount, skipJctCount, path)
                CircSkipfiles.append(CircSkipfile)
        return CircSkipfiles


def getSJ_out_tab(chimeralist):
    SJ_out_tab = []
    for fname in chimeralist:
        SJ_out_tab.append(fname.replace("Chimeric.out.junction", "SJ.out.tab"))
    return SJ_out_tab


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return "".join(random.choice(chars) for _ in range(size))


def checkBAMsorting(bamfiles):
    unsortedBAMs = []

    for file in bamfiles:

        bamfile = pysam.AlignmentFile(file, "rb")

        try:
            bamfile.check_index()
        except ValueError:
            print "BAM file %s has no index (%s.bai is missing)" % (file, file)
            logging.info("BAM file %s has no index (%s.bai is missing)" % (file, file))
            unsortedBAMs.append(file)
            break

        # Checking fort sorted reads is disabled to to issue #38:
        # https://github.com/dieterich-lab/DCC/issues/38

        # readcount = 0
        # readstarts = []
        #
        # # until_eof=False makes sure that we get the actual ordering in the file
        # for read in bamfile.fetch(until_eof=False):
        #
        #     readcount += 1
        #     readstarts.append(read.reference_start)
        #
        #     # sample only the first 100 reads
        #     if readcount > 100:
        #         break
        #
        # # close the BAM file again
        # bamfile.close()
        #
        # # compare sorted with unsorted list
        # if sorted(readstarts) != readstarts:
        #     # BAM file is not not sorted
        #     unsortedBAMs.append(file)

    # we return a list of unsorted files
    # -> if empty everything is sorted
    # -> if not, we have a list of files to look at for the user
    return unsortedBAMs


def wraphostgenecount(bamfile, tmp_dir, circ_coor, ref, countlinearsplicedreads=True):
    # create the Genecount object
    gc = Gc.Genecount(tmp_dir)

    # generate a unique thread ID
    tid = id_generator()

    # create an (temporary) output file based on tid and file name
    output = tmp_dir + "tmp_" + os.path.basename(bamfile) + "_" + tid + "_junction.linear"

    print "Counting host gene expression based on " \
          "detected and filtered circRNA coordinates for %s" % bamfile

    # launch the gene counting
    gc.comb_gen_count(circ_coor, bamfile, ref, output, countlinearsplicedreads)

    # return this input file's output name
    return output


def wrapfindcirc(files, tmp_dir, endTol, maxL, minL, strand=True, pairdendindependent=True, same=False):
    # create local instance
    f = Fc.Findcirc(endTol=endTol, maxL=maxL, minL=minL)

    # Start de novo circular RNA detection model
    sort = Fc.Sort()
    indx = id_generator()
    logging.info("started circRNA detection from file %s" % files)
    print "started circRNA detection from file %s" % files

    if same:
        circfilename = files + indx + ".circRNA"
    else:
        circfilename = files + ".circRNA"
    if pairdendindependent:
        f.printcircline(files, tmp_dir + "tmp_printcirclines." + indx)

        print "\t=> separating duplicates [%s]" % files
        f.sepDuplicates(tmp_dir + "tmp_printcirclines." + indx, tmp_dir + "tmp_duplicates." + indx,
                        tmp_dir + "tmp_nonduplicates." + indx)

        # Find small circles
        print "\t=> locating small circRNAs [%s]" % files
        f.smallcirc(tmp_dir + "tmp_duplicates." + indx, tmp_dir + "tmp_smallcircs." + indx)

        if strand:
            # Find normal circles
            print "\t=> locating circRNAs (stranded mode) [%s]" % files
            f.findcirc(tmp_dir + "tmp_nonduplicates." + indx, tmp_dir + "tmp_normalcircs." + indx, strand=True)
        else:
            print "\t=> locating circRNAs (unstranded mode) [%s]" % files
            f.findcirc(tmp_dir + "tmp_nonduplicates." + indx, tmp_dir + "tmp_normalcircs." + indx, strand=False)

        # Merge small and normal circles
        print "\t=> merging circRNAs [%s]" % files
        mergefiles(tmp_dir + "tmp_findcirc." + indx, tmp_dir + "tmp_smallcircs." + indx,
                   tmp_dir + "tmp_normalcircs." + indx)
    else:
        if strand:
            print "\t=> locating circRNAs (stranded mode) [%s]" % files
            f.findcirc(files, tmp_dir + "tmp_findcirc." + indx, strand=True)
        else:
            print "\t=> locating circRNAs (unstranded mode) [%s]" % files
            f.findcirc(files, tmp_dir + "tmp_findcirc." + indx, strand=False)

    # Sort
    if strand:
        print "\t=> sorting circRNAs (stranded mode) [%s]" % files
        sort.sort_count(tmp_dir + "tmp_findcirc." + indx, circfilename, strand=True)
    else:
        print "\t=> sorting circRNAs (unstranded mode) [%s]" % files
        sort.sort_count(tmp_dir + "tmp_findcirc." + indx, circfilename, strand=False)

    logging.info("finished circRNA detection from file %s" % files)
    print "finished circRNA detection from file %s" % files

    return circfilename


if __name__ == "__main__":
    main()
