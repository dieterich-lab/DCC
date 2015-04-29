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
            if line.split('\t')[9].split('.')[-1] == '1' or line.split('\t')[9].split('.') == '2':
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
                print 'Read %s has more than 2 count, please check.' % read
                try:
                    logging.warning('Read %s has more than 2 count, please check.' % read)
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

# Create a instance example

# a=Findcirc(strand=True,endTol=0,Chim_junc='GSM1328501Chimeric.out.junction',output='tmp')

# Sort the result and add read counts
class Sort(object):
    
    def __init__(self):
        """ 
        The Findcirc_output option takes the output of Class Findcirc. 
        Output is the sorted and bed format table with read counts information.
        """
        #self.strand=strand
        #self.Findcirc_output=open(Findcirc_output,'r').readlines()
        #self.output=output
        #self.sort_count(self.Findcirc_output)
      
    def count(self,sortedlist,strand=True):
        '''
        This function takes a sorted list of circRNAs, will be called by function circ_sort, count each circRNAs and return a bed format count table 
        '''
        cnt=collections.Counter()
        tmp_count = [] # store the bed format count circRNA count, but duplicated lines are not yet removed
        for itm in sortedlist:
            if strand:
                circs = (itm[0],itm[1],itm[2],itm[3])
            elif not strand:
                circs = (itm[0],itm[1],itm[2])
            else:
                print "Please specify correct strand information."
            cnt[circs] += 1
            itm.append(str(cnt[circs]))
            #tmp_count.append( [itm[0],itm[1],itm[2],itm[3],itm[7],itm[4],itm[5],itm[6]] )
            tmp_count.append( [itm[0],itm[1],itm[2],'.',itm[7],itm[3],itm[4],itm[5],itm[6]] )
        # reverse the order of tmp_count, so when remove duplicates, the one counted all the counts
        tmp_count.reverse()
        # remove duplicated lines
        tmp_count_dedup = []
        lines_seen = set() # holds the lines already seen
        for itm in tmp_count:
            if strand:
                tmp_itm = tuple((itm[0],itm[1],itm[2],itm[5]))
            elif not strand:
                tmp_itm = tuple((itm[0],itm[1],itm[2]))
            if tmp_itm not in lines_seen:
                tmp_count_dedup.append(itm)
                lines_seen.add(tmp_itm)
        # return the deduplicated count table (list)
        tmp_count_dedup.reverse()
        return tmp_count_dedup
        
            
    def sort_count(self,findcircOut,output,strand=True):
        # Equal to sort command in linux
        # file_list is the object of file read by readlines()
        file_list=open(findcircOut,'r').readlines()
        #print 'debugging'
        output = open(output,'w')
        tmp_sort = []
        for itm in file_list:
            line_tmp = itm.strip().split('\t')
            tmp_sort.append(line_tmp)
        tmp_sorted = sorted(tmp_sort,key=lambda x: (x[0],int(x[1]),int(x[2]),x[5]))
        sorted_count = self.count(tmp_sorted,strand=strand)
        output.writelines( '\t'.join(j) + '\n' for j in sorted_count )
        output.close()


#b=Sort(True,'tmp','tmp_out')

###############
# tmp to remove: Findcirc_output, tmp_sorted, tmp_sort, tmp_count
#
#a=open('tmp','r').readlines()
#tmp_a=[]
#for itm in a:
#    line_tmp = itm.strip().split('\t')
#    tmp_a.append(line_tmp)
#
#tmp_asortd=sorted(tmp_a,key=lambda x: (x[0],x[1],x[2],x[3]))
#cnt=collections.Counter()
#tmp_count_a=[]
#for itm in tmp_asortd:
#    circs = (itm[0],itm[1],itm[2])
#    cnt[circs] += 1
#    itm.append(str(cnt[circs]))
#    tmp_count_a.append([itm[0],itm[1],itm[2],itm[3],itm[7],itm[4],itm[5],itm[6]])
#    


## usage
#a = Findcirc(strand=True,endTol=5,maxL=50000,minL=50)
#a.findcir('GSM','tmp')
#b = Sort(True)
#b.sort_count('tmp','tmp.sorted')

# If merge three Chimeric.junc files (two mates mapped seperately and togethor), then print out circlines and select only duplicated lines, 
# then it contains both normal circles and small circles, because normal circles have one Chimeric.jun from togethor mapping and one 
# (only one) from sigle mate mapping. Small circles have two chimerics from two mates mapped seperately.

# If merge two Chimeric.junc from two seperately mapped mates, print circles and select duplicated lines. This are circ candidates, and should 
# not present in the joined mapped Chimeric.junc file.

# Possible solution:
# Map two mates seperately, then cat mate1.Chimeric.out.junciton mate2.Chimeric.out.junciton | print duplicates, use this duplicated file merge
# with the joined Chimeric.out.junction file, input in to the program.
