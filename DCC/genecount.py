# This module count the host gene counts
# Required: pysam
# Input: circRNA coordinates
#        aligned bam files
#        reference sequence
#   
# coor: coordinates         
# rs: result

import pysam
import os
import warnings
import sys

class Genecount(object):
    #def __init__(self,circ_coor, bamfile, ref, output):
    #    """
    #    @circ_coor: quoted string, content with format "chr1\tstart\tend"
    #    @bamfile: quoted string
    #    @ref: quoted string
    #    @output: quoted string
    #    """
    #    
    #    self.circ_coor = circ_coor
    #    self.bamfile = bamfile
    #    self.ref = ref
    #    self.output = output
    #    self.comb_gen_count(self.circ_coor, self.bamfile, self.ref, self.output)
        

    def countmapped(self, string):
        # This function takes the 5th column (a string with mapping information),
        # and return the count of '.' and ','. Which are mapped read count
        
        rs = string.count(',') + string.count('.')
        
        return str(rs)
    
    def countspliced(self,string):
        # return the count of '>' and '<'
        rs = string.count('>') + string.count('<')
        return str(rs)
    
    def getreadscount(self, mpileup, countmapped=False):
        # Input a mpileup result, which is a list
        # The count information is in the 5th column
        # count the number of mapped reads represented by ',' and '.' or, spliced reads represented by the number of '>' and '<' in the string
        count = []
        
        print "length mpileup: %d" % (len(mpileup))
        
        for itm in mpileup:
            tmp = itm.split('\t')
            if countmapped:
                count.append(tmp[0]+'\t'+tmp[1]+'\t'+self.countmapped(tmp[4]))
            else:
                count.append(tmp[0]+'\t'+tmp[1]+'\t'+self.countspliced(tmp[4]))
        
        #print count
        return count
            
        
    def genecount(self, circ_coor, bamfile, ref):
        """
        @circ_coor: quoted string, content with format "chr1\tstart\tend"
        @bamfile: quoted string
        @ref: quoted string
        """
        
        # process the circ_coordinates to left circ position and right circ position
        coor = open(circ_coor,'r').readlines()
        start_coor = open('start_coor','w')
        end_coor = open('end_coor','w')
        
        for line in coor:
            tmp = line.split('\t')
            start_coor.write(tmp[0]+'\t'+tmp[1]+'\n')
            end_coor.write(tmp[0]+'\t'+tmp[2]+'\n')
            
        # close postion files:
        start_coor.close()
        end_coor.close()
        
        print 'Started linear gene expression counting %s' % bamfile 
        
        # mpileup get the read counts of the start and end positions
        mpileup_start = pysam.mpileup(bamfile,'-f',ref,'-l','start_coor') # A list object
        mpileup_end = pysam.mpileup(bamfile,'-f',ref,'-l','end_coor')
        
        # get count
        startcount = self.getreadscount(mpileup_start,countmapped=True)
        endcount = self.getreadscount(mpileup_end,countmapped=True)
        
        # remove tmp files
        os.remove('start_coor')
        os.remove('end_coor')
        
        return startcount, endcount
        
    
    def submpileup(self,mpileup1,mpileup2,left=True):
        # Genome region: mpileup1 < mpileup2
        # used to calculate the spliced read counts difference
        # Input mpileup is the result of getreadscount function. 
        # Input in a list of elements with format str\tstr\tstr
        new_mpileup = []
        # need to consider that mpileup1 and mpileup2 are not the same length
        if len(mpileup1)==len(mpileup2):
            for itm in range(len(mpileup1)):
                mpileup1_count = int(mpileup1[itm].split('\t')[2])
                mpileup2_count = int(mpileup2[itm].split('\t')[2])
                shift = mpileup1_count - mpileup2_count
                if shift<0:
                    shift=0
                new_mpileup.append(mpileup1[itm].split('\t')[0]+'\t'+mpileup1[itm].split('\t')[1]+'\t'+str(shift))
        else:
            mpileup1_dict = {}
            mpileup2_dict = {}
            for itm in mpileup1:
                mpileup1_dict[(itm.split('\t')[0]+'\t'+itm.split('\t')[1])] = int(itm.split('\t')[2])
            for itm in mpileup2:
                mpileup2_dict[(itm.split('\t')[0]+'\t'+str(int(itm.split('\t')[1])-1))] = int(itm.split('\t')[2])

            for itm in mpileup1_dict:
                if itm in mpileup2_dict:
                    if left:
                        shift = mpileup1_dict[itm] - mpileup2_dict[itm]
                        if shift<0:
                            shift=0
                        new_mpileup.append(itm.split('\t')[0]+'\t'+str(int(itm.split('\t')[1])+1)+'\t'+str(shift))
                    else:
                        shift = mpileup2_dict[itm] - mpileup1_dict[itm]
                        if shift<0:
                            shift=0
                        new_mpileup.append(itm.split('\t')[0]+'\t'+itm.split('\t')[1]+'\t'+str(shift))
                    
        return new_mpileup
    
    def linearsplicedreadscount(self,circ_coor,bamfile,ref,header=True):
        # Count linear spliced reads
        # process the circ_coordinates to left circ position and right circ position
        if header:
            coor = open(circ_coor,'r').readlines()[1:]
        else:
            coor = open(circ_coor,'r').readlines()
        start_coor = open('start_coor','w')
        start_coor_1 = open('start_coor_1','w')
        end_coor = open('end_coor','w')
        end_coor_1 = open('end_coor_1','w')
        
        for line in coor:
            tmp = line.split('\t')
            start_coor.write(tmp[0]+'\t'+tmp[1]+'\n')
            start_coor_1.write(tmp[0]+'\t'+str(int(tmp[1])-1)+'\n')
            end_coor.write(tmp[0]+'\t'+tmp[2]+'\n')
            end_coor_1.write(tmp[0]+'\t'+str(int(tmp[2])+1)+'\n')
            
        # close postion files:
        start_coor.close()
        start_coor_1.close()
        end_coor.close()
        end_coor_1.close()
        print 'Started linear spliced reads count for %s' % bamfile
        
        # mpileup get the number of spliced reads at circle start position and (start-1) postion.
        mpileup_start = pysam.mpileup(bamfile,'-f',ref,'-l','start_coor') # A list object
        mpileup_start_1 = pysam.mpileup(bamfile,'-f',ref,'-l','start_coor_1')
                
        # mpileup get the number of spliced reads at circle end position and (end+1) postion.
        mpileup_end = pysam.mpileup(bamfile,'-f',ref,'-l','end_coor')
        mpileup_end_1 = pysam.mpileup(bamfile,'-f',ref,'-l','end_coor_1')
        
        # get count
        startcount = self.submpileup(self.getreadscount(mpileup_start_1),self.getreadscount(mpileup_start))
        endcount = self.submpileup(self.getreadscount(mpileup_end),self.getreadscount(mpileup_end_1),left=False)

        # remove tmp files
        os.remove('start_coor')
        os.remove('start_coor_1')
        os.remove('end_coor')
        os.remove('end_coor_1')
        
        return startcount, endcount
    
    def comb_gen_count(self, circ_coor, bamfile, ref, output, countlinearsplicedreads=True):
        indx = open(circ_coor,'r').readlines()[1:] # make sure are tab delimited
        #tmp_start = open(options.startfile,'r') # make sure are tab delimited
        #tmp_end = open(options.endfile,'r') # make sure are taa delimited
        coordinates_indx_start = [] # This is the chromosome name and start position of the original bed file list, like Lvr_F_104_filtered_candid
        coordinates_indx_end = [] # This is the chromosome name and end position of the original bed file list, like Lvr_F_104_filtered_candid
        coordinates_start = {} # This is the chromosome name and start position of the counted read counts bed file list, like tmp_readcounts_start
        coordinates_end = {} # This is the chromosome name and start position of the counted read counts bed file list, like tmp_readcounts_end
        count_start=[] # Store the read counts of start positions
        count_end=[] # Store the read counts of the end positions
        coordinates=[] # Store chr, start, end information from circRNAs candidates file
        
        if countlinearsplicedreads:
            tmp_start, tmp_end = self.linearsplicedreadscount(circ_coor, bamfile, ref)
        else:
            # call genecount to get the start and end positon read counts
            tmp_start, tmp_end = self.genecount(circ_coor, bamfile, ref)
        
        print 'Ended linear gene expression counting %s' % bamfile
        
        for line in tmp_start:
            tmp_start_split=line.split('\t')
            coordinates_start[(tmp_start_split[0],tmp_start_split[1])] = tmp_start_split[2].strip()
    
        for line in tmp_end:
            tmp_end_split=line.split('\t')
            coordinates_end[(tmp_end_split[0],tmp_end_split[1])] = tmp_end_split[2].strip()
    
        for line in indx:
            indx_split = line.split('\t')
            coordinates_indx_start.append((indx_split[0],indx_split[1]))
            coordinates_indx_end.append((indx_split[0],indx_split[2].strip('\n')))
            coordinates.append((indx_split[0],indx_split[1],indx_split[2].strip('\n'))) 
    
        for itm in coordinates_indx_start:
            if itm in coordinates_start:
                count_start.append(coordinates_start[itm])
            else:
                count_start.append('0')
                warnings.warn('WARNING: circRNA start position '+str(itm)+' does not have mapped read counts, treated as 0', stacklevel=3)
                
        for itm in coordinates_indx_end:
            if itm in coordinates_end:
                count_end.append(coordinates_end[itm])
            else:
                count_end.append('0')
                warnings.warn('WARNING: circRNA end position '+str(itm)+' does not have mapped read counts, treated as 0', stacklevel=3)
        
        # write count table
        count_table = open(output,'w')
        if (len(coordinates)==len(count_start)==len(count_end)):
            for i in range(len(coordinates)):
                res = list(coordinates[i])+[count_start[i],count_end[i],str(int(round(float(count_start[i])+float(count_end[i]))/2))]
                count_table.write(('\t').join(res)+'\n')
        else:
            sys.exit('read count number does not match with number of circRNAs andidates')
    
        # close files
        #tmp_start.close()
        #tmp_end.close()
        count_table.close()

                
#a=Genecount('test2','GSM1328501dedupped.bam','rn6.fa','outfile')