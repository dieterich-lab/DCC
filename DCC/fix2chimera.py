#!/bin/python

from findcircRNA import Findcirc
# If paire end data are used, this module fix the rolling circRNA could not detect by single run problem.
fc = Findcirc(1,1,1)

class Fix2Chimera(object):

	def fixreadname(self,chimeric_junction_file,output_file):
		# Because sometimes, for example -I flag of fastq-dump will add ".1" and ".2" read suffices.
		# "======== NOTE! ======= The default flag for matched pairedend reads is suffice '.1' and '.2'. If not, please change the code!!! "
		chimeric_junction = open(chimeric_junction_file,'r')
		output = open(output_file,'w')
		for line in chimeric_junction:
			#split field 10
			suffice = False
			s = [str(i).strip() for i in line.split('\t')]
			#print s
			if len(s[9].split('.')[-1])==1:
				suffice = True
			if suffice:
				s[9] = '.'.join(s[9].split('.')[:-1])
				output.write( '\t'.join(s) )
			else:
				output.write( '\t'.join(s) )

		chimeric_junction.close()
		output.close()

	def fixmate2(self,chimeric_junction_mate2,output_file):
		# Fix the orientation of mate2
		chimeric_junction = open(chimeric_junction_mate2,'r')

		def modify_junctiontype(junctiontype):
			if junctiontype == '0':
				return '0'
			else:
				return str(3-int(junctiontype))

		output = open(output_file,'w')
		for line in chimeric_junction:
			line_split = line.split('\t')
			if line_split[2] == '+':
				line_split[2] = '-'
				line_split[5] = '-'
				line_split[6] = modify_junctiontype(line_split[6])
				output.write( '\t'.join((line_split)) )
			else:
				line_split[2] = '+'
				line_split[5] = '+'
				line_split[6] = modify_junctiontype(line_split[6])
				output.write( '\t'.join((line_split)) )

		chimeric_junction.close()
		output.close()

#	def printduplicated(self,filein,basedon):

	def concatenatefiles(self,output_file,*fnames):
		with open(output_file,'w') as output:
			for fname in fnames:
				with open(fname) as infile:
					for line in infile:
						output.write(line)

	def fixation(self,mate1,mate2,joined,output_file):
		# take chimeric.junction.out from two mate mapping and mate-joined mapping, 
		# return fixed chimeric.out.junction

		# First, fix mate2
		self.fixmate2(mate2,mate2+'.fixed')

		# Second, merge two mate files, select duplicates
		self.concatenatefiles('tmp_merged',mate1,mate2+'.fixed')
		fc.sepDuplicates('tmp_merged','tmp_twochimera','tmp_onechimera')
		self.concatenatefiles(output_file,'tmp_twochimera',joined)





