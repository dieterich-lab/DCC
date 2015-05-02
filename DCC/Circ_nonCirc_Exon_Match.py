#!/usr/bin/python

# given circRNA coordinates file and exon annotation gtf file, print out the adjacent exon coordinates and exon_id

import pybedtools
import HTSeq

class CircNonCircExon(object):

	def print_start_end_file(self, circcoordinates):
		circ = open(circcoordinates,'r')
		start_bed = open('start.bed','w')
		end_bed = open('end.bed','w')
		for lin in circ:
			lin_split = lin.split('\t')
			start_bed.write(lin_split[0]+'\t'+lin_split[1]+'\t'+lin_split[1]+'\n')
			end_bed.write(lin_split[0]+'\t'+lin_split[2]+'\t'+lin_split[2]+'\n')
		circ.close()
		start_bed.close()
		end_bed.close()


	def modifyExon_id(self, gtf_file):
		# write custom_exon_id as transcript_id:exon_number
		gtf = HTSeq.GFF_Reader(gtf_file, end_included=True)
		new_gtf = open(gtf_file.replace('.gtf','')+'.modified.gtf','w')
		for feature in gtf:
			custom_exon_id = '; custom_exon_id'+' '+'"'+feature.attr['transcript_id']+':'+feature.attr['exon_number']+'"'
			new_gtf.write( feature.get_gff_line().strip('\n')+custom_exon_id+'\n' )
		new_gtf.close()

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

	def intersectcirc(self, circ_file, modified_gtf_file):
		# imput the result file of print_start_end_file
		import pybedtools
		#intersectBed -a start.bed -b Drosophila_melanogaster.BDGP5.75.exon_id.dedup.gtf -wa -wb -loj > tmpintersect.2
		circ = pybedtools.BedTool(circ_file)
		gtf = pybedtools.BedTool(modified_gtf_file)
		intersectfile = circ.intersect(gtf,wa=True,wb=True,loj=True)
		# Store circExons as: circle start or end intervals as key, custom_exon_id as value
		circExons = {}
		for lin in intersectfile:
			lin_split = str(lin).split('\t')
			if lin_split[11].strip('\n') == '.':
				#lin_split[11] = ''
				pass
			else:
				circExons.setdefault( HTSeq.GenomicInterval(lin_split[0],int(lin_split[1]),int(lin_split[2]),lin_split[9]), set() ).add( HTSeq.parse_GFF_attribute_string(lin_split[11])['custom_exon_id'] )
			#circExons.setdefault( HTSeq.GenomicInterval(lin_split[0],int(lin_split[1]),int(lin_split[2]),lin_split[9]), [] ).append( { HTSeq.GenomicInterval(lin_split[3],int(lin_split[6]),int(lin_split[7]),lin_split[9]):HTSeq.parse_GFF_attribute_string(lin_split[11]) })
		return circExons
		#intersectfile.moveto('intersectfile')

	def printuniq(self, Infile):
		f = open(Infile,'r').readlines()
		keys = []
		for lin in f:
			lin_split = lin.split('\t')
			keys.append(lin_split[0]+'\t'+lin_split[1]+'\t'+lin_split[2])
		for lin in f:
			lin_split = lin.split('\t')
			if keys.count(lin_split[0]+'\t'+lin_split[1]+'\t'+lin_split[2]) == 1:
				print lin.strip('\n')

	def readgtf(self, gtf_file):
		# store nonCircExons based on transcript_id and exon_number with all its annotations from different transcripts
		gtf = HTSeq.GFF_Reader(gtf_file, end_included=True)
		nonCircExons = {} # A list of dictionaries
		for feature in gtf:
			nonCircExons.setdefault(feature.iv,[]).append(feature.attr)
		return nonCircExons

	def getAdjacent(self, custom_exon_id, start=True):
		if start:
			exon_number = int(custom_exon_id.split(':')[1])-1
		else:
			exon_number = int(custom_exon_id.split(':')[1])+1
		if exon_number==0:
			return 'None'
		else:
			transcript_id = custom_exon_id.split(':')[0]
			return transcript_id+':'+str(exon_number)

	def readHTSeqCount(self, HTSeqCount, exon_id2custom_exon_id):
		Count = open(HTSeqCount,'r').readlines()
		Count_custom_exon_id = {}
		for lin in Count:
			lin_split = lin.split('\t')
			if lin_split[0] == '__no_feature':
				break
			else:
				Count_custom_exon_id[exon_id2custom_exon_id[lin_split[0]]] = int(lin_split[1])
		return Count_custom_exon_id
	
	def findcircAdjacent(self, circExons,Custom_exon_id2Iv,Iv2Custom_exon_id):
		circAdjacentExons = {}
		for key in circExons.keys():
			for ids in circExons[key]:
				try:
					interval = Custom_exon_id2Iv[self.getAdjacent(ids)]
					interval2ids = Iv2Custom_exon_id[interval]
				except KeyError: # CircExon is the start or end of that transcript
					interval = 'None'
					interval2ids = 'None'				
				circAdjacentExons.setdefault(key,[]).extend(interval2ids)
				# From custom_exon_id find out the interval, from interval find out all the custom_exon_ids for that interval
		return circAdjacentExons

	def printCounts(self,Exons,Count_custom_exon_id,Custom_exon_id2Length):
		# Print the counts of circexons and adjacentexons
		# Exons: dictionaries with intervals as key, custom_exon_id as values
		ExonCounts = {}
		for key in Exons.keys():
			counts = []
			for ids in Exons[key]: # If for circAdjacentExons, ids here is a list
				try:
					counts.append(ids+' FPKM='+str(float(Count_custom_exon_id[ids])/Custom_exon_id2Length[ids]))
				except KeyError:
					pass
			ExonCounts[key] = counts
			#print str(key)+'\t'+'\t'.join((counts))
		return ExonCounts

		# For circExons, no problem because all the circExons are selected from deduplicates, but those one region has more 
		# than one count (one region has more than one distinct exon) shoul left out, because the reads counting process is ambigous.
    
	def readNonUniqgtf(self, NonUniqgtf):
		gtf = HTSeq.GFF_Reader(NonUniqgtf, end_included=True)
		Iv2Custom_exon_id = {}
		Custom_exon_id2Iv = {}
		Custom_exon_id2Length = {}
		for feature in gtf:
			Iv2Custom_exon_id.setdefault(feature.iv,[]).append(feature.attr['custom_exon_id'])
			Custom_exon_id2Iv.setdefault(feature.attr['custom_exon_id'],feature.iv)
			Custom_exon_id2Length.setdefault(feature.attr['custom_exon_id'],feature.iv.length)
		return (Iv2Custom_exon_id, Custom_exon_id2Iv, Custom_exon_id2Length)

	def printresults(self,circCount,circAdjacentCount,circExons,circAdjacentExons):
		result = open('exonFPKM','w')
		result_clean = open('exonFPKM_clean','w')
		for key in circCount:
			if len(circCount[key])>1 or len(circAdjacentCount[key])>1:
				pass
			else:
				result.write( str(key)+': '+'\t'+str(circCount[key])+'\t'+str(circAdjacentCount[key])+'\n' )
				try:
					circ = str(circCount[key][0]).split('=')[1]
				except IndexError:
					circ = 'NA'
				try:
					circAdjacent = str(circAdjacentCount[key][0]).split('=')[1]
				except IndexError:
					circAdjacent = 'NA'

				result_clean.write( circ+'\t'+ circAdjacent+'\n' )
		result.close()
		result_clean.close()




	#circExons or circAdjacentExons

    # def findAdjacent(self, circExons, nonCircExons):
    # 	# Scann through circExons and find the corresponding adjacent nonCircExons, store in dictionary of dictionaries with circ start or end coordinates 
    # 	for circ in circExons.keys():
    # 		for exons in circExons(circ):
    # 			# exons are also dictionaries with keys being exons intervals
    # 			exons.keys[0] # based on this look for adjacent exons





