# given circRNA coordinates file and exon annotation gtf file, print out the adjacent exon coordinates and exon_id

import os
import re

import HTSeq

from IntervalTree import IntervalTree


class CircNonCircExon(object):
    def __init__(self, tmp_dir):
        self.tmp_dir = tmp_dir

    def print_start_end_file(self, circcoordinates):
        # Print start.bed and end.bed
        # Get start and end corresponding relationship
        start2end = {}
        circ = open(circcoordinates, 'r')
        start_bed = open(self.tmp_dir + 'tmp_start.bed', 'w')
        end_bed = open(self.tmp_dir + 'tmp_end.bed', 'w')
        header = True
        for lin in circ:
            lin = lin.rstrip()
            if header:
                header = False
            else:
                lin_split = lin.split('\t')
                start_bed.write(
                    lin_split[0] + '\t' + lin_split[1] + '\t' + lin_split[1] + '\t' + '.' + '\t' + '.' + '\t' +
                    lin_split[5] + '\n')
                end_bed.write(
                    lin_split[0] + '\t' + lin_split[2] + '\t' + lin_split[2] + '\t' + '.' + '\t' + '.' + '\t' +
                    lin_split[5] + '\n')
                # Here is only for strand data, because here store GenomicInterval object to
                # guide the matching of circStartAdjacentExon with circEndAdjacentExon,
                # If the data is non-strand, the matching will suffer from uncertain strand
                start2end.setdefault(
                    HTSeq.GenomicInterval(lin_split[0], int(lin_split[1]), int(lin_split[1]), str(lin_split[5])), []).append(
                    HTSeq.GenomicInterval(lin_split[0], int(lin_split[2]), int(lin_split[2]), str(lin_split[5])))
        circ.close()
        start_bed.close()
        end_bed.close()
        return start2end

    def select_exon(self, gtf_file):
        gtf = HTSeq.GFF_Reader(gtf_file, end_included=True)
        new_gtf = open(self.tmp_dir + 'tmp_' + os.path.basename(gtf_file) + '.exon.sorted', 'w')
        gtf_exon = []
        for feature in gtf:
            # Select only exon line
            if feature.type == 'exon':
                gtf_exon.append(feature.get_gff_line().split('\t'))
            else:
                pass
        gtf_exon_sorted = sorted(gtf_exon, key=lambda x: (x[0], int(x[3]), int(x[4])))
        gtf_exon_sorted = ['\t'.join(s) for s in gtf_exon_sorted]
        new_gtf.writelines(gtf_exon_sorted)
        new_gtf.close()

    def modifyExon_id(self, exon_gtf_file):
        rtrn = True
        # write custom_exon_id as transcript_id:exon_number
        gtf = HTSeq.GFF_Reader(exon_gtf_file, end_included=True)
        # gff = True
        new_gtf = open(self.tmp_dir + os.path.basename(exon_gtf_file) + '.modified', 'w')
        exon_number = {}
        for feature in gtf:
            # if gff:
            # 	try:
            # 		feature.attr['ID']
            # 	except KeyError:
            # 		gff = False
            # else:
            # 	pass
            # First look for transcript_id
            try:
                if feature.attr['transcript_id'] in exon_number:
                    # if feature.iv.strand == '+':
                    exon_number[feature.attr['transcript_id']] = exon_number[feature.attr['transcript_id']] + 1
                # else:
                # exon_number[feature.attr['transcript_id']] = exon_number[feature.attr['transcript_id']] - 1
                else:
                    # if feature.iv.strand == '+':
                    exon_number[feature.attr['transcript_id']] = 1
                # else:
                # exon_number[feature.attr['transcript_id']] = -1
                # if gff:
                # 	custom_exon_id = ';custom_exon_id'+'='+feature.attr['transcript_id']+':'+str(exon_number[feature.attr['transcript_id']])
                # else:
                custom_exon_id = '; custom_exon_id' + ' ' + '"' + feature.attr['transcript_id'] + ':' + str(
                    exon_number[feature.attr['transcript_id']]) + '"'
            except KeyError:
                # Try assume gff format
                try:
                    if feature.attr['Parent'] in exon_number:
                        # if feature.iv.strand == '+':
                        exon_number[feature.attr['Parent']] = exon_number[feature.attr['Parent']] + 1
                    # else:
                    #	exon_number[feature.attr['Parent']] = exon_number[feature.attr['Parent']] - 1
                    else:
                        # if feature.iv.strand == '+':
                        exon_number[feature.attr['Parent']] = 1
                    # else:
                    #	exon_number[feature.attr['Parent']] = -1
                    custom_exon_id = ';custom_exon_id' + '=' + feature.attr['Parent'] + ':' + str(
                        exon_number[feature.attr['Parent']])
                except (KeyError, TypeError):
                    print (
                        'DCC confused with the annotation, cannot determine CircSkip junctions. '
                        'If gtf file provided, one or two of the features cannot find: transcript_id. '
                        'If gff file provided, cannot determine Parent feature.')
                    rtrn = False
                    break
            new_gtf.write(feature.get_gff_line().strip('\n') + custom_exon_id + '\n')

        # Do not print non exon features
        #### MAKE sure only modify and interact with exons!!!!!!!!  FIrst get only exons!!!
        ####  for gff3 files, go for IDs!!!!!!		# Solved
        new_gtf.close()
        return rtrn

    def readModifiedgtf(self, modifiedgtf):
        # store a dictionary which custom_exon_id are keys, exon_id are values.
        gtf = HTSeq.GFF_Reader(modifiedgtf, end_included=True)
        custom_exon_id2exon_id = {}
        for feature in gtf:
            custom_exon_id2exon_id[feature.attr['custom_exon_id']] = feature.attr['exon_id']
        return custom_exon_id2exon_id

    def readuniqgtf(self, uniqgtf):
        # Requires exon_id are uniq in the file
        gtf = HTSeq.GFF_Reader(uniqgtf, end_included=True)
        exon_id2custom_exon_id = {}
        for feature in gtf:
            exon_id2custom_exon_id[feature.attr['exon_id']] = feature.attr['custom_exon_id']
        return exon_id2custom_exon_id

    def intersectcirc(self, circ_file, modified_gtf_file, strand=True, isStartBED=True):
        # input the result file of print_start_end_file
        input_bed_file = open(circ_file).readlines()
        exon_gtf_file = HTSeq.GFF_Reader(modified_gtf_file, end_included=True)
        gtf_exon_sorted = IntervalTree()
        for feature in exon_gtf_file:
            row = feature.attr
            current_bed_interval = feature.iv
            gtf_exon_sorted.insert(current_bed_interval, annotation=row)

        circ_exon_set = {}
        for bed_line in input_bed_file:
            bed_field = bed_line.split('\t')
            custom_exon_list = []

            # we add 1bp in order for intersect to work correctly
            # different case for start or end bed file
            if isStartBED:
                start = int(bed_field[1])
                end = int(bed_field[1]) + 1
            else:
                start = int(bed_field[1]) - 1
                end = int(bed_field[1])

            # in order for the intersect to work, we need at least 1bp frame size
            current_bed_interval = HTSeq.GenomicInterval(bed_field[0],
                                                         start,
                                                         end,
                                                         bed_field[5].strip()
                                                         )

            # for later processing however, we again need the "0" bp frame window
            insert_bed_interval = HTSeq.GenomicInterval(bed_field[0],
                                                        int(bed_field[1]),
                                                        int(bed_field[2]),
                                                        bed_field[5].strip()
                                                        )
            # extract all customs exons
            gtf_exon_sorted.intersect(current_bed_interval,
                                      lambda x: custom_exon_list.append(x.annotation['custom_exon_id'])
                                      )

            if custom_exon_list:  # if we found one or more custom exons
                for custom_exon in custom_exon_list:  # go through the list
                    circ_exon_set.setdefault(insert_bed_interval, set()).add(custom_exon)  # and add them to the set

        # return the filled set
        return circ_exon_set

    def printuniq(self, Infile):
        f = open(Infile, 'r').readlines()
        keys = []
        for lin in f:
            lin_split = lin.split('\t')
            keys.append(lin_split[0] + '\t' + lin_split[1] + '\t' + lin_split[2])
        for lin in f:
            lin_split = lin.split('\t')
            if keys.count(lin_split[0] + '\t' + lin_split[1] + '\t' + lin_split[2]) == 1:
                print lin.strip('\n')

    def readgtf(self, gtf_file):
        # store nonCircExons based on transcript_id and exon_number with all its annotations from different transcripts
        gtf = HTSeq.GFF_Reader(gtf_file, end_included=True)
        nonCircExons = {}  # A list of dictionaries
        for feature in gtf:
            nonCircExons.setdefault(feature.iv, []).append(feature.attr)
        return nonCircExons

    def getAdjacent(self, custom_exon_id, start=True, reverse=False):
        # Need to determine the oder. (Some exon ids increasing from first exon to last of the transcript,
        # irrelevant to strand. But some id increase with coordinates, in this case, for - strand, exon id will
        # ACTUALLY decrease from 5' to 3'. First determine whether exon id is reverse oder for - strand.
        #
        # Need a function to determine order. !!!
        # Solved by sort gtf and assign id to exon based on occuring order. Thus for - strand, order reverse.
        #
        if reverse:
            if start:
                exon_number = int(custom_exon_id.split(':')[1]) - 1
            else:
                exon_number = int(custom_exon_id.split(':')[1]) + 1
        else:
            if start:
                exon_number = int(custom_exon_id.split(':')[1]) - 1
            else:
                exon_number = int(custom_exon_id.split(':')[1]) + 1
        if exon_number == -1:
            return 'None'
        else:
            transcript_id = custom_exon_id.split(':')[0]
            return transcript_id + ':' + str(exon_number)

    def getgtforder(self, modified_gtf_file):
        gtf = HTSeq.GFF_Reader(modified_gtf_file, end_included=True)
        seen = {}
        for feature in gtf:
            if feature.type == 'exon' and feature.iv.strand == '-':
                try:
                    parent = feature.attr['transcript_id']
                    seen[parent] = feature  # GTF
                except KeyError:
                    parent = feature.attr['Parent']
                    seen[parent] = feature  # GFF3
                if parent in seen:
                    # found one transcript with multiple exon, try determine order now
                    if (seen[parent].iv.start - feature.iv.start) * (
                                self.getcustom_id_num(seen[parent].attr['custom_exon_id']) - self.getcustom_id_num(
                                feature.attr['custom_exon_id'])) > 0:
                        return True
                    else:
                        return False


    def getIDnum(self, ID):
        # Get the ID number of ID feature of GFF3 file
        res = re.findall(r'\d+', ID)
        if len(res) > 1:
            return None
            print('Cannot correctly determine ID.')
        else:
            return res[0]

    def getcustom_id_num(self, custom_id):
        return int(custom_id.split(':')[1])

    def readHTSeqCount(self, HTSeqCount, exon_id2custom_exon_id):
        Count = open(HTSeqCount, 'r').readlines()
        Count_custom_exon_id = {}
        for lin in Count:
            lin_split = lin.split('\t')
            if lin_split[0] == '__no_feature':
                break
            else:
                Count_custom_exon_id[exon_id2custom_exon_id[lin_split[0]]] = int(lin_split[1])
        return Count_custom_exon_id

    def findcircAdjacent(self, circExons, Custom_exon_id2Iv, Iv2Custom_exon_id, start=True):
        circAdjacentExons = {}
        circAdjacentExonsIv = {}
        for key in circExons.keys():
            for ids in circExons[key]:
                try:
                    interval = Custom_exon_id2Iv[self.getAdjacent(ids, start=start)]
                    interval2ids = Iv2Custom_exon_id[interval]
                except KeyError:  # CircExon is the start or end of that transcript
                    interval = 'None'
                    interval2ids = 'None'
                circAdjacentExons.setdefault(key, []).extend(interval2ids)
                circAdjacentExonsIv.setdefault(key, []).append(interval)
                # From custom_exon_id find out the interval, from interval find out all the custom_exon_ids for that interval
        return circAdjacentExons, circAdjacentExonsIv

    def printCounts(self, Exons, Count_custom_exon_id, Custom_exon_id2Length):
        # Print the counts of circexons and adjacentexons
        # Exons: dictionaries with intervals as key, custom_exon_id as values
        ExonCounts = {}
        for key in Exons.keys():
            counts = []
            for ids in Exons[key]:  # If for circAdjacentExons, ids here is a list
                try:
                    counts.append(ids + ' FPKM=' + str(float(Count_custom_exon_id[ids]) / Custom_exon_id2Length[ids]))
                except KeyError:
                    pass
            ExonCounts[key] = counts
        return ExonCounts

    # For circExons, no problem because all the circExons are selected from deduplicates, but those one region has more
    # than one count (one region has more than one distinct exon) shoul left out, because the reads counting process is ambigous.

    def readNonUniqgtf(self, NonUniqgtf):
        gtf = HTSeq.GFF_Reader(NonUniqgtf, end_included=True)
        Iv2Custom_exon_id = {}
        Custom_exon_id2Iv = {}
        Custom_exon_id2Length = {}
        for feature in gtf:
            Iv2Custom_exon_id.setdefault(feature.iv, []).append(feature.attr['custom_exon_id'])
            Custom_exon_id2Iv.setdefault(feature.attr['custom_exon_id'], feature.iv)
            Custom_exon_id2Length.setdefault(feature.attr['custom_exon_id'], feature.iv.length)
        return (Iv2Custom_exon_id, Custom_exon_id2Iv, Custom_exon_id2Length)

    def printresults(self, circCount, circAdjacentCount, circExons, circAdjacentExons, prefix):
        result = open(prefix + 'exonFPKM', 'w')
        result_clean = open(prefix + 'exonFPKM_clean', 'w')
        for key in circCount:
            # This step discard those start/end mapped to more than one overlaping exons
            if len(circCount[key]) > 1 or len(circAdjacentCount[key]) > 1:
                pass
            else:
                result.write(str(key) + ': ' + '\t' + str(circCount[key]) + '\t' + str(circAdjacentCount[key]) + '\n')
                try:
                    circ = str(circCount[key][0]).split('=')[1]
                except IndexError:
                    circ = 'NA'
                try:
                    circAdjacent = str(circAdjacentCount[key][0]).split('=')[1]
                except IndexError:
                    circAdjacent = 'NA'

                result_clean.write(circ + '\t' + circAdjacent + '\n')
        result.close()
        result_clean.close()

    def exonskipjunction(self, circStartAdjacentExonsIv, circEndAdjacentExonsIv, start2end, strand=True):
        # A list of CircRNA interval to exon skip (circ skip) junctions
        junctions = {}
        for key in circStartAdjacentExonsIv:
            start = set()
            for itv0 in circStartAdjacentExonsIv[key]:
                try:
                    start.add(str(itv0.end))  # Last base of left adjacent non-circEXON
                except AttributeError:
                    pass

            # Find the end interval
            endiv = start2end.get(key,
                                  None)  # A list, if more than one circ start from the same position, but have different ending.
            if endiv:
                for itv1 in endiv:
                    end = set()
                    # store the circ
                    circ = HTSeq.GenomicInterval(itv1.chrom, key.start, itv1.end, key.strand)
                    # but not every itv have adjacentexon
                    try:
                        for itv2 in circEndAdjacentExonsIv[itv1]:
                            try:
                                end.add(str(itv2.start))  # First base of right adjacent non-circEXON
                            except AttributeError:
                                pass
                    except KeyError:
                        pass

                    if len(start) > 0 and len(end) > 0:
                        for i in start:
                            for j in end:
                                junctions.setdefault(circ, []).append(itv1.chrom + '\t' +
                                                                      str(int(i) + 1) + '\t'
                                                                      + j + "\t" + key.strand)
        return junctions

    def readSJ_out_tab(self, SJ_out_tab):
        # read SJ.out.tab, store coordinates and read counts into a dictionary
        junctionReadCount = {}
        try:
            sj = open(SJ_out_tab, 'r')
            for lin in sj:

                lin_split = lin.split('\t')

                if lin_split[3] == "1":
                    strand = "+"
                elif lin_split[3] == "2":
                    strand = "-"
                else:
                    strand = "?"

                junctionReadCount[lin_split[0] + '\t' +
                                  lin_split[1] + '\t' +
                                  lin_split[2] + '\t' +
                                  strand] = lin_split[6]
            sj.close()
        except IOError:
            print 'Do you have SJ.out.tab files in your sample folder? DCC cannot find it.'
        return junctionReadCount

    def getskipjunctionCount(self, exonskipjunctions, junctionReadCount):
        # A list of CircRNA interval to exon skip (circ skip) junction read counts
        skipJctCount = {}
        for key in exonskipjunctions:
            junctions = exonskipjunctions[key]
            count = []
            for jct in junctions:
                try:
                    count.append(jct.split('\t')[0] + ':' +
                                 jct.split('\t')[1] + '-' +
                                 jct.split('\t')[2] +
                                 jct.split('\t')[3] + ':' +
                                 junctionReadCount[jct])
                except KeyError:
                    pass
            if len(count) > 0:
                counts = ';'.join((count))
                skipJctCount[key] = counts
        return skipJctCount

    def readcircCount(self, circRNACount):
        circCount = {}
        Countfile = open(circRNACount, 'r')
        for lin in Countfile:
            lin_split = lin.split('\t')
            itv = HTSeq.GenomicInterval(lin_split[0], int(lin_split[1]), int(lin_split[2]), lin_split[5])
            circCount[itv] = lin_split[4]
        Countfile.close()
        return circCount

    def printCirc_Skip_Count(self, circCount, skipJctCount, prefix):
        Circ_Skip_Count = []
        for key in skipJctCount:
            try:
                count = skipJctCount[key]
                Circ_Skip_Count.append([key.chrom, str(key.start), str(key.end), count, circCount[key], key.strand])
            except KeyError:
                pass

        # sort
        Circ_Skip_Count = sorted(Circ_Skip_Count, key=lambda x: (x[0], int(x[1]), int(x[2])))

        # print "Circ_Skip_Count:" + str(Circ_Skip_Count)
        Circ_Skip_Count_file = open(prefix + 'CircSkipJunction', 'w')

        for sublist in Circ_Skip_Count:
            # for list_item in sublist:
            Circ_Skip_Count_file.write("\t".join(sublist) + "\n")

        Circ_Skip_Count_file.close()

        return prefix + 'CircSkipJunction'
