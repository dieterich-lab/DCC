#! /usr/bin/env python

# If paire end data are used, this module fix the rolling circRNA could not detect by single run problem.

class Fix2Chimera(object):
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
        chimeric_junction = open(chimeric_junction_mate2, 'r')

        def modify_junctiontype(junctiontype):
            if junctiontype == '0':
                return '0'
            else:
                return str(3 - int(junctiontype))

        output = open(output_file, 'w')
        for line in chimeric_junction:
            line_split = line.split('\t')
            if line_split[2] == '+':
                line_split[2] = '-'
                line_split[5] = '-'
            else:
                line_split[2] = '+'
                line_split[5] = '+'
            line_split[6] = modify_junctiontype(line_split[6])
            output.write('\t'.join((line_split[0], line_split[4], line_split[2], line_split[3], line_split[1],
                                    line_split[5], line_split[6], line_split[7], line_split[8], line_split[9],
                                    line_split[10], line_split[11], line_split[12], line_split[13])))

        chimeric_junction.close()
        output.close()

    #	def printduplicated(self,filein,basedon):

    def concatenatefiles(self, output_file, *fnames):
        import shutil
        destination = open(output_file, 'wb')
        for fname in fnames:
            shutil.copyfileobj(open(fname, 'rb'), destination)
        destination.close()

    def fixation(self, mate1, mate2, joined, output_file):
        # take chimeric.junction.out from two mate mapping and mate-joined mapping,
        # return fixed chimeric.out.junction

        # First, fix mate2
        self.fixmate2(mate2, mate2 + '.fixed')

        # Second, merge two mate files, select duplicates
        self.concatenatefiles('_tmp_DCC/tmp_merged', mate1, mate2 + '.fixed')
        self.printduplicates('_tmp_DCC/tmp_merged', '_tmp_DCC/tmp_twochimera', field=10)
        self.concatenatefiles(output_file, '_tmp_DCC/tmp_twochimera', joined)

    def printduplicates(self, infile, dupfile, field=10):
        inputfile = open(infile, 'r')
        dup = open(dupfile, 'w')
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
            if key in seen:
                if not seen[key][0]:
                    dup.write(seen[key][1])
                    seen[key] = (True, seen[key][1])
                dup.write(line)
            else:
                seen[key] = (False, line)
        inputfile.close()
        dup.close()
