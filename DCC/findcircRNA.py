#!/usr/bin/env python

import re
import collections

class Findcirc(object):
    # Initialize some parameters
    def __init__(self,endTol,minL,maxL):
   	#self.strand = strand
   	self.endTol = endTol
   	#self.output = output
   	self.maxL = int(maxL)
   	self.minL = int(minL)
   	#self.findcirc(Chim_junc)
    
    def cigarGenomicDist(self,cig):
   	C=re.findall('[a-zA-Z]',cig)
   	L=re.findall('\-?[0-9]+',cig)
   	n=len(L)
   	g=0
   	for i in range(0,n):
  		if C[i]!='S' and C[i]!="I":
 			g = g+int(L[i])
   	return g
	
    def printcircline(self,Chim_junc,output):
        junctionfile=open(Chim_junc,'r')
        outfile=open(output,'w')
        for line in junctionfile:
            L=line.split('\t')
            if  int(L[6])>=0 and L[0]==L[3] and L[2]==L[5] and ((L[2]=='-' and int(L[4])>int(L[1]) and self.minL<(int(L[4])-int(L[1]))<self.maxL) or (L[2]=='+' and int(L[1])>int(L[4]) and self.minL<(int(L[1])-int(L[4]))<self.maxL)):
                if (L[2]=='+' and (int(L[10])+self.endTol)>int(L[4]) and (int(L[12])+self.cigarGenomicDist(L[13])-self.endTol)<=int(L[1])) or (L[2]=='-' and (int(L[12])+self.endTol)>int(L[1]) and (int(L[10])+self.cigarGenomicDist(L[11])-self.endTol)<=int(L[4])):
                    outfile.write(line)
        outfile.close()
        junctionfile.close()
        
        
        
    # The bug here is: sepDuplicates, all the duplicates come from two mates are chimera, no problem, but the nonduplicates not all from joined mapping. For some reason, some reads chimerically mapped seperately, but not joinedly, although they are not two mates chimeria

    def sepDuplicates(self,Chim_junc,duplicates,nonduplicates):
        # seperate out duplicated reads, means both mates are chimera
        # Chim_junc is the output of printcircline function
        # Duplicated reads are appear in a scrambled order, the real strand read come later than the reversed strand
        
        junctionfile=open(Chim_junc,'r')
        dup = open(duplicates,'w')
        nondup = open(nonduplicates,'w')
        
        reads = []
        lines = [] # A list of all lines
        suffice = False
        for line in junctionfile:
            if len(line.split('\t')[9].split('.')[-1]) == 1:
                suffice = True
            if suffice:
                readname = '.'.join(line.split('\t')[9].split('.')[:-1])
            else:
                readname = line.split('\t')[9]
            reads.append(readname)
            lines.append(line)
        
        for indx,read in enumerate(reads):
            if reads.count(read) == 2:
                dup.write(lines[indx])
            elif reads.count(read) > 2:
                print 'Read %s has more than 2 count.' % read
                try:
                    logging.warning('Read %s has more than 2 count.' % read)
                except NameError:
                    pass
            else:
                nondup.write(lines[indx])

        junctionfile.close()
        dup.close()
        nondup.close()

    def smallcirc(self,duplicates,output,strand=True):
        '''
        Find small circRNAs in pairedend data where both mates are chimeric. Input is duplicates file by sepDuplicates.
        '''
        
        # The second occure read has the 'real' strand of that circle.
        
        dup = open(duplicates).readlines()
        outfile = open(output,'w')

        collect = []
        for line in dup:
            L = line.split('\t')
            if L[2] == '+':
                identifier = (L[0],L[1],L[3],L[4],L[9])
            elif L[2] == '-':
                identifier = (L[0],L[4],L[3],L[1],L[9])
            #print identifier
            if identifier in collect:
                if L[2] == '+':
                    if strand:
                        if L[6]=='0':
          		    res=[L[0],str(int(L[4])+1),str(int(L[1])-1),'-','0',L[7],L[8]] 
           	        else:
          		    res=[L[0],str(int(L[4])+1),str(int(L[1])-1),'-',str(3-int(L[6])),L[7],L[8]]
          		outfile.write(('\t').join(res) + '\n')
                    else:
                   	res=[L[0],str(int(L[4])+1),str(int(L[1])-1),L[2],L[6],L[7],L[8]]
                   	outfile.write(('\t').join(res)+ '\n')
                if L[2] == '-':
                    if strand:
                        if L[6]=='0':
                            res=[L[0],str(int(L[1])+1),str(int(L[4])-1),'+','0',L[7],L[8]]
                        else:
                            res=[L[0],str(int(L[1])+1),str(int(L[4])-1),'+',str(3-int(L[6])),L[7],L[8]]
			outfile.write(('\t').join(res)+ '\n')
                    else:
                        res=[L[0],str(int(L[1])+1),str(int(L[4])-1),L[2],L[6],L[7],L[8]]
			outfile.write(('\t').join(res)+ '\n')
            else:
                collect.append(identifier) # Identifiers for duplicates
	
            ## All switch to plus strand
            #if L[2] == '-':
            #    collect.append((L[0],L[4],'+',L[3],L[1],'+',L[9]))
            #else:
            #    collect.append(tuple(L[0:6]+[L[9]]))
        #
        #for line in dup:
        #    # Only print out with '+' strand
        #    L = line.split('\t')
        #    toTest = tuple(L[0:6]+[L[9]])
        #    if collect.count(toTest) == 2:
        #        #res = line
        #        res = [L[0],L[4],L[1],L[2],L[6],L[7],L[8]]
        #        outfile.write(('\t').join(res) + '\n')
        #        #outfile.write(res)

        outfile.close()


    def findcirc(self,Chim_junc,output,strand=True):
	junctionfile=open(Chim_junc,'r')
	outfile=open(output,'w')

	for line in junctionfile:
		L=line.split('\t')
		if  int(L[6])>=0 and L[0]==L[3] and L[2]==L[5] and ((L[2]=='-' and int(L[4])>int(L[1]) and self.minL < (int(L[4])-int(L[1])) < self.maxL) or (L[2]=='+' and int(L[1])>int(L[4]) and self.minL<(int(L[1])-int(L[4]))<self.maxL)):
			if (L[2]=='+' and (int(L[10])+self.endTol)>int(L[4]) and (int(L[12])+self.cigarGenomicDist(L[13])-self.endTol)<=int(L[1])):
				if strand:
					if L[6]=='0':
						res=[L[0],str(int(L[4])+1),str(int(L[1])-1),'-','0',L[7],L[8]]
					else:
						res=[L[0],str(int(L[4])+1),str(int(L[1])-1),'-',str(3-int(L[6])),L[7],L[8]]
					outfile.write(('\t').join(res) + '\n')
				else:
					res=[L[0],str(int(L[4])+1),str(int(L[1])-1),L[2],L[6],L[7],L[8]]
					outfile.write(('\t').join(res)+ '\n')
			if L[2]=='-' and (int(L[12])+self.endTol)>int(L[1]) and (int(L[10])+self.cigarGenomicDist(L[11])-self.endTol)<=int(L[4]):
				if strand:
					if L[6]=='0':
						res=[L[0],str(int(L[1])+1),str(int(L[4])-1),'+','0',L[7],L[8]]
					else:
						res=[L[0],str(int(L[1])+1),str(int(L[4])-1),'+',str(3-int(L[6])),L[7],L[8]]
					outfile.write(('\t').join(res)+ '\n')
				else:
					res=[L[0],str(int(L[1])+1),str(int(L[4])-1),L[2],L[6],L[7],L[8]]
					outfile.write(('\t').join(res)+ '\n')
	outfile.close()
        junctionfile.close()
#
