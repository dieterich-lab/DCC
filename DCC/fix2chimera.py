#! /usr/bin/env python2

# If paire end data are used, this module fix the rolling circRNA could not detect by single run problem.

import os
import sys


class Fix2Chimera(object):
    def __init__(self, tmp_dir):
        self.tmp_dir = tmp_dir

    def fixreadname(self, chimeric_junction_file, output_file):
        # Because sometimes, for example -I flag of fastq-dump will add ".1" and ".2" read suffices.
        # "======== NOTE! ======= The default flag for matched pairedend reads is suffice '.1' and '.2'. "
        chimeric_junction = open(chimeric_junction_file, 'r')
        output = open(output_file, 'w')
        for line in chimeric_junction:
            # split field 10
            suffice = False
            s = [str(i).strip() for i in line.split('\t')]
            # print s
            if len(s[9].split('.')[-1]) == 1:
                suffice = True
            if suffice:
                s[9] = '.'.join(s[9].split('.')[:-1])
                output.write('\t'.join(s))
            else:
                output.write('\t'.join(s))

        chimeric_junction.close()
        output.close()

    def fixmate2(self, chimeric_junction_mate2, output_file):
        # Fix the orientation of mate2

        if not os.path.isfile(chimeric_junction_mate2):
            sys.exit("ERROR: File " + str(chimeric_junction_mate2) + " is missing!")

        chimeric_junction = open(chimeric_junction_mate2, 'r')

        def modify_junctiontype(junctiontype):
            if junctiontype == '0':
                return '0'
            else:
                return str(3 - int(junctiontype))

        linecnt = 1
        output = open(output_file, 'w')
        for line in chimeric_junction:
            line = line.rstrip()
            line_split = line.split('\t')

            if line_split[0] == "chr_donorA":
                continue
            # check if the row has all fields
            if len(line_split) < 14:
                print(("WARNING: File " + str(chimeric_junction_mate2) + ", line " + str(linecnt)
                       + " does not contain all features."))
                print(("WARNING: " + str(chimeric_junction_mate2) + " is probably corrupt."))
                print(("WARNING: Offending line: " + str(line)))

            linecnt += 1

            if line_split[2] == '+':
                line_split[2] = '-'
                line_split[5] = '-'
            else:
                line_split[2] = '+'
                line_split[5] = '+'
            line_split[6] = modify_junctiontype(line_split[6])
            output.write('\t'.join((line_split[0], line_split[4], line_split[2], line_split[3], line_split[1],
                                    line_split[5], line_split[6], line_split[7], line_split[8], line_split[9],
                                    line_split[10], line_split[11], line_split[12], line_split[13]))+'\n')

        chimeric_junction.close()
        output.close()

    def concatenatefiles(self, output_file, *fnames):
        import shutil
        destination = open(output_file, 'wb')
        for fname in fnames:
            if not os.path.isfile(fname):
                sys.exit("ERROR: File " + str(fname) + " is missing!")
            else:
                shutil.copyfileobj(open(fname, 'rb'), destination)
        destination.close()

    def fixchimerics(self, mate1, mate2, joined, output_file):
        # take chimeric.junction.out from two mate mapping and mate-joined mapping,
        # return fixed chimeric.out.junction

        # First, fix mate2
        self.fixmate2(mate2, mate2 + '.fixed')

        # Second, merge two mate files, select duplicates
        self.concatenatefiles(self.tmp_dir + 'tmp_merged', mate1, mate2 + '.fixed')  # does not care if files are empty
        self.printduplicates(self.tmp_dir + 'tmp_merged', self.tmp_dir + 'tmp_twochimera',
                             field=10)  # TODO: field is static?

        if not os.path.isfile(self.tmp_dir + 'tmp_twochimera'):
            sys.exit("PE-independent mode selected but all corresponding junctions files are empty!")

        self.concatenatefiles(output_file, self.tmp_dir + 'tmp_twochimera', joined)

    def printduplicates(self, merged, duplicates, field=10):
        inputfile = None
        dup = None

        if not os.path.isfile(merged):
            sys.exit("ERROR: File " + str(merged) + " is missing!")
        elif os.stat(merged).st_size == 0:
            print(("WARNING: File " + str(merged) + " is empty!"))
        else:
            try:
                inputfile = open(merged, 'r')
                dup = open(duplicates, 'w')
                seen = dict()
                # check read name suffice
                suffice = False
                if len(inputfile.readline().split('\t')[field - 1].split('.')[-1]) == 1:
                    suffice = True
                for line in inputfile:
                    line_split = line.split('\t')
                    if suffice:
                        key = line_split[field - 1][:-2]
                    else:
                        key = line_split[field - 1]
                    if key in seen:  # has been added to dict
                        if not seen[key][0]:  # but it has not yet been set to True (already written)
                            dup.write(seen[key][1])  # write the content of the key
                            seen[key] = (True, seen[key][1])  # set key to seen and written
                        dup.write(line)  # write the complete original line to duplicate file
                    else:
                        seen[key] = (False, line)  # executed when we see the read name first, sets false in dict
            finally:
                inputfile.close()
                dup.close()
