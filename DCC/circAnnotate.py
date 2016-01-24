# This module annotate the circRNA with gene names, and also filter the circRNA candidates by requiring "CircRNA could not from more than two annotated linear gene."
import warnings
import logging
import os
import re
import HTSeq
from IntervalTree import IntervalTree

class CircAnnotate(object):
    def __init__(self,strand=True):
        self.strand = strand
        
    def selectGeneGtf(self,gtf_file):
        # select gene features for gtf or gff annotation file
        gtf = HTSeq.GFF_Reader(gtf_file, end_included=True)
        annotation_tree = IntervalTree()
        for feature in gtf:
            # Select only exon line
            iv = feature.iv
            try:
                row = feature.attr
                row['type'] = feature.type
            except:
                row = feature.get_gff_line()
            annotation_tree.insert(iv, annotation=row)
        return annotation_tree

    def annotate_one_interval(self,interval,annotation_tree,what='gene'):
        out = []
        annotation_tree.intersect(interval,lambda x: out.append(x.annotation))
        annotation = self.searchGeneName(out,what=what)
        return annotation

    def annotate(self,circfile,annotation_tree,output):
        # the circRNA file should be in a bed format, have chr\tstart\tend\t'.'\tjunctiontype\tstrand
        # The annotation tree should be a IntervalTree object
        
        # check the input
        with open(circfile,'r') as tmpcirc:
            tmpsplit = tmpcirc.readline().split('\t')
            if len(tmpsplit) != 6:
                warnings.warn('Input circRNA file is not the desired bed6 format!')
                logging.warning('Input circRNA file is not the desired bed6 format!')
            ncol = len(tmpsplit)
        
        # Annotate with Interval tree algorithm
        out = open(output,'w')
        circ_reagions = HTSeq.BED_Reader(circfile)
        for circ in circ_reagions:
            annotation = self.annotate_one_interval(circ.iv,annotation_tree,what='gene')
            out.write('\t'.join([circ.iv.chrom,str(circ.iv.start),str(circ.iv.end),annotation,str(int(circ.score)),circ.iv.strand])+'\n')
        out.close()

        #self.printbycolumns('_tmp_DCC/tmp_AnnotatedUnsorted',output,order=[1,2,3,6,4,5])
        #os.remove('_tmp_DCC/tmp_AnnotatedUnsorted')
          
    def uniqstring(self,strings,sep=','):
        string = set(strings.split(sep))
        string = sep.join(string)
        return string

    def annotateregions(self,circfile, annotation_tree,output):
        # Annotate with regions (Exon, intron, intergenic)
        # create left and right circle bundary bedfiles: chr\tstart\tstart  chr\tend\tend
        circ = open(circfile,'r').readlines()
        new_CircCoordinates = open(output,'w')
        for line in circ:
            line_split = line.split('\t')
            iv_left = HTSeq.GenomicInterval(line_split[0],int(line_split[1]),int(line_split[1]),line_split[5].strip('\n'))
            iv_right = HTSeq.GenomicInterval(line_split[0],int(line_split[2]),int(line_split[2]),line_split[5].strip('\n'))
            iv_left_annotation = self.annotate_one_interval(iv_left,annotation_tree,what='reagion')
            if not iv_left_annotation:
                iv_left_annotation = '.'
            iv_right_annotation = self.annotate_one_interval(iv_right,annotation_tree,what='reagion')
            if not iv_right_annotation:
                iv_right_annotation = '.'
            overall_annotation = self.uniqstring(iv_right_annotation+','+iv_right_annotation)
            iv_left_annotation = self.readRegionAnnotate(iv_left_annotation)
            iv_right_annotation = self.readRegionAnnotate(iv_right_annotation)
            new_CircCoordinates.write(line.strip('\n')+'\t'+iv_left_annotation+'-'+iv_right_annotation+'\t'+overall_annotation+'\n')
        new_CircCoordinates.close()

    
    def readRegionAnnotate(self,annotations):
        if 'exon' in annotations:
            return 'exon'
        elif len(annotations) > 1 and annotations != 'region':
            return 'intron'
        elif len(annotations)==1 or annotations == 'region':
            return 'intergenic'
                    
    def filtbygene(self,circ2filter,output):
        # This funtion filter the circs base on: circRNAs should not come from two genes.
        out = open(output,'w')
        with open(circ2filter,'r') as circ:
	   for line in circ:
		tmp = line.split('\t')
		n=tmp[3].split(',')
		try:
		  if len(n)==1:
		      out.write(line)
		except IndexError:
			pass
	out.close()


    def printbycolumns(self,fileIn,fileOut,order=[],sep='\t',fillempty=True):
        tmpIn = open(fileIn,'r')
        tmpOut = open(fileOut,'w')
        for lines in tmpIn:
            tmpsplit = [x.strip() for x in lines.split(sep)]
            if fillempty:
                tmpsplit = ['.' if x=='' else x for x in tmpsplit]
            # Get gene_id or gene_name annotation
            tmpsplit[5] = self.searchGeneName(tmpsplit[5])                 
            tmpOut.write('\t'.join([tmpsplit[int(i)-1] for i in order])+'\n')
        tmpIn.close()
        tmpOut.close()
        
        
    def searchGeneName1(self,annotations):
        # Search for gene_name in gtf annotation, if gene_name cannot be found, look for gene_id
        # input example: gene_id "ENSG00000187634"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
        ann = ','.join(list(set(re.findall(r'gene_name\=?\s?"([^;]*)"\;',annotations))))
        if len(ann)==0:
            # Look for "gene=", which is used in gff3 format
            ann = ','.join(list(set(re.findall(r'gene\=?\s?"([^;]*)"\;',annotations))))
            if len(ann)==0:
                # Look for gene_id
                ann = ','.join(list(set(re.findall(r'gene_id\=?\s?"([^;]*)"\;',annotations))))
                if len(ann)==0:
                    # Look for transcript_id
                    ann = ','.join(list(set(re.findall(r'transcript_id\=?\s?"([^;]*)"\;',annotations))))
        if len(ann)==0:
            ann = 'N/A'
        return ann
        
    def searchGeneName(self,annotations,what='gene'):
        if annotations == '.':
            genes = 'N/A'
        elif what=='gene':
            collect = set()
            for annotation in annotations:
                # Search for gene_name which is used by ensembl gtf annotation
                try:
                    gene = annotation['gene_name']
                except TypeError:
                    gene = self.searchGeneName1(annotation)
                except KeyError:
                    # Search for gene, which might used in GFF annotation
                    try:
                        gene = annotation['gene']
                    except:
                        # Search for gene_id
                        try:
                            gene = annotation['gene_id']
                        except:
                            try:
                                gene = annotation['transcript_id']
                            except:
                                gene = 'N/A'
                collect.add(gene)
            # Collapse all genes togethor
            if len(collect) > 1:
                try:
                    collect.remove('N/A')
                except KeyError:
                    pass
        else:
            # annotate reagion
            collect = set()
            for annotation in annotations:
                try:
                    gene = annotation['type']
                except TypeError:
                    gene = annotation.split('\t')[2]
                except:
                    gene = 'N/A'
                collect.add(gene)
        genes = ','.join(collect)
        if not genes:
            # empty string
            genes = '.'
        return genes
        
        
        
