#! /usr/bin/env python

# This script used to combine individual circRNA count files to a single count table

# invoke with column nr to extract as first parameter followed by
# filenames. The files should all have the same number of rows

import sys
#from optparse import OptionParser
import pybedtools
import re
import os

#usage = "usage: %prog -c column2select -t genecountOrcirccount -o outputfile"
#parser = OptionParser(usage=usage)
#parser.add_option("-k", "--col", dest="col",
#                  help="Column to select and paste together")
#parser.add_option("-c", "--circ", action='store_true', dest="type", default=True,
#                  help="The input are circRNA counts")
#parser.add_option("-h", "--host", action='store_false', dest="type", default=False,
#                  help="The input are host gene counts")
#parser.add_option("-o", "--output", dest="output",
#                  help="The output gene expression count table")                  
#(options, args) = parser.parse_args()


class Combine(object):
#    def __init__(self, Alist, genetype, output, debug=False):
#        ''' A model to combine circRNA coordinates from all samples and combine individual counts to one table'''
#        # Alist is a list of files to be combined, either mapped circRNA or host gene counts
#        self.Alist = Alist
#        self.genetype = genetype
#        self.output = output
#
#        if self.genetype == 'circ':
#            self.comb_coor(self.Alist)
#            self.map('tmp_coordinates', self.Alist)
#            res = self.combine(self.Alist, col=8, circ=True)
#        elif self.genetype == 'host':
#            res = self.combine(self.Alist, col=6, circ=False)
#        else:
#            print "genetype could only be 'circ' or 'host'!"  
#        self.writeouput(self.output, res)
#        
#        # Delete tmp file is not debug model
#        if debug:
#            pass
#        else:
#            # Delete intermidiate circRNA count files
#            pattern1 = '\.+circRNA\.*'
#            self.deletfile(os.getcwd(),pattern1)
#            pattern2 = '\.+hostgene$'
#            self.deletfile(os.getcwd(),pattern2)
        
    # Input a list of file names, which are the result of sorted circRNA file
    def comb_coor(self, circfiles, strand=True):
        '''
        Combine coordinates of all samples to one.
        '''
        coordinates = open('tmp_coordinates','w')
        coors = set()
        coorsDict = {} # Use all except the strand and junction type information as key, to uniq the list.
        for files in circfiles:
            onefile = open(files,'r')
            for lines in onefile:
                tmp = lines.split('\t')
                
                if strand:
                    coors.add(tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+'.'+'\t'+tmp[6]+'\t'+tmp[5]+'\n')
                else:
                    coorsDict[tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+'.'+'\t'] = tmp[6]+'\t'+tmp[5]+'\n'                    
            onefile.close()

        if not strand:
            coors = ['{}{}'.format(key,value) for key,value in coorsDict.iteritems()]

        coorsSorted = self.sortBed(coors,retList=True)
        for itm in coorsSorted:
            coordinates.write('\t'.join(itm))
    
    def map(self, coordinates, Alist, strand=True):
        '''
        Take the combined coordinates and list of circRNA files.
        '''
        bed1 = pybedtools.BedTool(coordinates)
        for files in Alist:
            bed2 = pybedtools.BedTool(files)
            if strand:
                mapped = bed1.map(bed2,s=True, f=1,r=True)
            else:
                mapped = bed1.map(bed2,s=False, f=1,r=True)
            mapped.moveto(files+'mapped')

    def deletfile(self, dirt, pattern):
        # First check wheter the input is a list of files or a regular expression string
        if isinstance(pattern,str):        
            # A list to store names of deleted files
            deleted = []
            for f in os.listdir(dirt):
                if re.search(pattern,f):
                    os.remove(os.path.join(dirt,f))
                    deleted.append(f)
        elif isinstance(pattern,list):
            for f in pattern:
                os.remove(os.path.join(dirt,f))
                deleted = pattern
        return deleted
    
    def sortBed(self, bedfile, retList=False):
        # The input could be a list with one element per line of bedfile, or a file
        # First check whether the input is a list
        if isinstance(bedfile, list) or isinstance(bedfile, tuple) or isinstance(bedfile, set):
            bedfile = list(bedfile)
            # check dimention
            if isinstance(bedfile[0],list):
                # Sort
                bedfileSorted = sorted(bedfile,key=lambda x: (x[0],int(x[1]),int(x[2]),x[5]))
            else:
                # Split
                for indx, elem in enumerate(bedfile):
                    bedfile[indx] = elem.split('\t')
                bedfileSorted = sorted(bedfile,key=lambda x: (x[0],int(x[1]),int(x[2]),x[5]))
            
        elif isinstance(bedfile, str):
            try:
                fileOpen = open(bedfile,'r').readlines()
                for indx, elem in enumerate(fileOpen):
                    fileOpen[indx] = elem.split('\t')
            except IOError:
                sys.exit('Input a quoted filename or a list.')
            else:
                bedfileSorted = sorted(fileOpen,key=lambda x: (x[0],int(x[1]),int(x[2]),x[5]))
                # Write to file
                if retList:
                    pass
                else:
                    output = open('bedfileSorted','w')
                    output.writelines(bedfileSorted)
        else:
            sys.exit('Can only sort list, set, tuple or file')
        return bedfileSorted
    
    
    def combine(self, Alist, col=7, circ=True ):
        # Alist is a list of files to be combined, either mapped circRNA or host gene counts
        col = int(col)
        res = {}
        
        for file_name in Alist:
            for line_nr, line in enumerate(open(file_name)):
                line_split = line.strip('\n').split('\t')
    
                if circ: # input are circRNA counts
                    # check whether the fifth field has any number or empty, if empty, substitute with 0
                    if line_split[col-1] == '.':
                        line_split[col-1] = '0'
                    #elif line_split[col-1] != '':
                        # res.setdefault(line_nr,[]) will set res[line_nr]=[] if res does not have key 'line_nr', 
                        #when the key 'line_nr' exist, it will return the value which is a list [something], and will append the read count later
                    #    pass
                    #elif len(line_split) != 5:
                    #    sys.exit("Number of fields in the individual count file is not correct.")
                    res.setdefault(line_nr, ['\t'.join(line_split[:3]),line_split[4]]).append(line_split[col-1])
                    
                else: # input are host gene counts
                    res.setdefault(line_nr, ['\t'.join(line_split[:3])]).append(line_split[col-1])
                
        return res

    def writeouput(self, output, res):        
        outfile = open(output,'w')
        for line_nr in sorted(res):
            outfile.write("%s\n" % '\t'.join(res[line_nr]))
        outfile.close()
        