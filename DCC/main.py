#!/usr/bin/env python
# -*- coding: utf-8 -*-

import findcircRNA as FC
import CombineCounts as CC
import argparse
#from optparse import OptionParser, OptionGroup
import os
import sys
#import re
import logging
import circAnnotate
import time
#import circFilter


parser = argparse.ArgumentParser(prog='main',formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', description='Contact jun.cheng@age.mpg.de')


parser.add_argument('--version', action='version', version='%(prog)s 2.0')
parser.add_argument("Input", metavar='Input', nargs="+",
                  help="Input of the chimeric.out.junction file from STAR. Alternatively a file with each sample name with path per line.")
#parser.add_argument("-O", "--output", dest="output", 
#                  help="Tab delimited outputfile, order the same with input: \
#                  chr\tstart\tend\tstand\tcount\tjunctiontype")
parser.add_argument("-deb", "--debug", dest="debug", action='store_true', default=False,
                  help="Once specified, temp files will not be deleted.")
 
                                                                                                     
group = parser.add_argument_group("Find circRNA Options","Options to find circRNAs from STAR output.")
group.add_argument("-D", "--detect", action='store_true', dest="detect", default=False,
                  help="If specified, the program will start detecting circRNAs from chimeric junction.")                    
group.add_argument("-S", action='store_true', dest="strand", default=True,
                  help="Specify when the library is stranded [Default].")
group.add_argument("-N", action='store_false', dest="strand",
                  help="Specify when the library is non-stranded")
group.add_argument("-E", "--endTol", dest="endTol", type=int, default= 5, choices=range(0,10),
                  help="Maximum tolerance of reads extrend over junction sites: Interger")
group.add_argument("-m", "--maximum", dest="max", type=int, default = 1000000,
                  help="The maximum range length of candidate circRNA allowed.")   
group.add_argument("-n", "--minimum", dest="min", type=int, default = 30,
                  help="The minimum range length of candidate circRNA allowed.")
group.add_argument("-an", "--annotate",dest="annotate",
                  help="Gene annotation file in gtf format, to annotate circRNAs，if provided, the circRNA intervals will be annotated (default with gene_id).")
#group.add_argument("-gf", "--getfasta", dest="getfasta",
#                  help="Get fasta file of circular RNAs. If a exon annotation file is provided, the circular RNA sequence will only contain annotated exons, otherwise whole sequence.")
group.add_argument("-P", "--pe", action='store_true', dest="pairedend", default=False,
                 help="Whether or not the data is paired end, if yes, paired end information is used for filtering:\
                 Boolean")
parser.add_argument_group(group)                  
  
      
group = parser.add_argument_group("Filtering Options", "Options to filter the circRNA candidates.")                              
group.add_argument("-F", "--filter", action='store_true', dest="filter", default=False,
                  help="If specified, the program will start filtering model.")
group.add_argument("-M", "--chrM", action='store_true', dest="chrM",default=False,
                  help="If specified, candidates from mitochondria chromosome will be removed.")
#group.add_argument("-J", "--junction", dest="junction",
#                  help="Provide a coustom junction file in gtf format, if only specify as True, only GT/AG or CT/AC junction will be considered.")
group.add_argument("-R", "--rep_file", dest="rep_file",
                  help="Coustom repetitive region file in gtf format.") 
group.add_argument("-L", "--Ln", dest="length", type=int, default=50,
                  help="Minimum length to check for repetitive regions, default 50.")                                                      
group.add_argument('-Nr', nargs=4, type=int, metavar=('level0', 'level1','threshold0','threshold1'), default=[4,2,5,5], help='Minimum read counts required for circRNAs with junction type 0 and junction type 1, \
                    Minimum number of samples above the expression corresponding level')
group.add_argument("-fg", "--filterbygene", action='store_true', dest="filterbygene", default=False,
                  help="Filter by gene annotation, candidates are not allowed to cover more than one gene.")                                        
parser.add_argument_group(group)


group = parser.add_argument_group("Host gene count Options", "Options to count host gene expression.")
group.add_argument("-G", "--gene", action='store_true', dest="gene", default=False,
                  help="If specified, the program will count host gene expression given circRNA coordinates.")
group.add_argument("-C", "--circ", dest="circ",
                  help="circRNA coordinates, or any tab delimited file with first three columns circRNA coordinates: chr\tstart\tend. If not specified, program take output of previous detection.")                  
group.add_argument("-B", "--bam", dest="bam", nargs = '+',
                  help="The mapped bam file where host gene count counted from, the same oder with the input file.")
group.add_argument("-A", "--refseq", dest="refseq",         
                  help="Reference sequnce fasta file.") 
#group.add_argument("-seq", "--seq", dest="seq",         
#                  help="Get circRNA sequence as fasta file.")
parser.add_argument_group(group)

                                    
options= parser.parse_args()

def main():
    ## Check whether the input is a list or a file
    #if isinstance(options.Input,list):
    #    pass
    #print type(options.Input)
                                            
    #if (options.detect + options.gene) > 1:
    #    sys.exit("Only one model can be specified, detect (-D), host gene count (-G).")
    
    timestr = time.strftime("%Y%m%d-%H%M%S")
    logging.basicConfig(filename='main.log'+timestr, filemode='w',level=logging.DEBUG,format='%(asctime)s %(message)s')
    
    #root = logging.getLogger()
    #root.setLevel(logging.DEBUG)
    #
    #ch = logging.StreamHandler(sys.stdout)
    #ch.setLevel(logging.DEBUG)
    #ch.filename='main.log'
    #ch.filemode='w'
    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    #ch.setFormatter(formatter)
    #root.addHandler(ch)

    logging.info(' '.join(sys.argv))
    logging.info('Program Started')

    # Get input file names
    filenames = [getfilename(name) for name in options.Input]
    samplelist = '\t'.join(filenames)
    
    # check whether the junction file names have duplicates
    same = False
    if len(set(filenames)) != len(options.Input):
        logging.info('Input filenames have duplicates, add number suffix in input order to output files for distinction.')
        print ('Input filenames have duplicates, add number suffix in input order to output files for distinction.')
        same = True
    
    # Make instance
    cm = CC.Combine()
    circAnn = circAnnotate.CircAnnotate(strand=options.strand)
    
    if options.detect:
        logging.info('Program start to detect circRNAs')
        
        # Start de novo circular RNA detection model
        # Create instances
        f = FC.Findcirc(endTol=options.endTol,maxL=options.max,minL=options.min)
        sort = FC.Sort()
        
        circfiles = [] # A list for .circRNA file names
        
        if options.pairedend:
            print '===== Please make sure that you mapped both the paired mates togethor and seperately, and run the fixation scripts!!! ====='
            logging.info("===== Please make sure that you mapped both the paired mates togethor and seperately, and run the fixation scripts!!! =====")
        
        
        def wrapfindcirc(files,output,strand=True,pairdend=True):
            if pairdend:
                f.printcircline(files,'tmp_printcirclines')
                f.sepDuplicates('tmp_printcirclines','tmp_duplicates','tmp_nonduplicates')
                # Find small circles
                f.smallcirc('tmp_duplicates','tmp_smallcircs')
                if strand:
                    # Find normal circles
                    f.findcirc('tmp_nonduplicates','tmp_normalcircs',strand=True)
                else:
                    f.findcirc('tmp_nonduplicates','tmp_normalcircs',strand=False)
                # Merge small and normal circles
                mergefiles(('tmp_smallcircs','tmp_normalcircs'),'tmp_findcirc')
            else:
                if strand:
                    f.findcirc(files,'tmp_findcirc',strand=True)
                else:
                    f.findcirc(files,'tmp_findcirc',strand=False)
            # Sort
            if strand:
                sort.sort_count('tmp_findcirc',output,strand=True)
            else:
                sort.sort_count('tmp_findcirc',output,strand=False)
            
                            
        for indx, files in enumerate(options.Input):
            
            if same:
                circfilename = getfilename(files)+str(indx)+'.circRNA'
            else:
                circfilename = getfilename(files)+'.circRNA'
            circfiles.append(circfilename)
    
                                                    
            if options.strand:
                logging.info( 'strand' )
                print 'strand'
                
                if options.pairedend:
                    wrapfindcirc(files,circfilename,strand=True,pairdend=True)
                else:
                    wrapfindcirc(files,circfilename,strand=True,pairdend=False)

            elif not options.strand:
                logging.info( 'nonstrand' )
                print 'nonstrand'
       	        
                if options.pairedend:
                    wrapfindcirc(files,circfilename,strand=False,pairdend=True)
                else:
                    wrapfindcirc(files,circfilename,strand=False,pairdend=False)
        #        
        #try:
        #    os.remove('tmp_findcirc')
        #    os.remove('tmp_printcirclines')
        #    os.remove('tmp_duplicates')
        #    os.remove('tmp_nonduplicates')
        #    os.remove('tmp_smallcircs')
        #except OSError:
        #    pass
            
            
        ### Combine the individual count files
        # Create a list of '.circRNA' file names
        logging.info('Start to combine individual circRNA read counts')

        if options.strand:
            cm.comb_coor(circfiles,strand=True)
            cm.map('tmp_coordinates', circfiles, strand=True)
        else:
            cm.comb_coor(circfiles,strand=False)
            cm.map('tmp_coordinates', circfiles, strand=False)
        res = cm.combine([files+'mapped' for files in circfiles],col=7,circ=True)
        
        if options.filter:
            cm.writeouput('tmp_circCount', res)
            if options.annotate:
                logging.info('Write in annotation.')
                circAnn.annotate('tmp_coordinates',options.annotate,'tmp_coordinatesannotated')
                os.remove('tmp_coordinates')
                os.rename('tmp_coordinatesannotated','tmp_coordinates')
        else:
            cm.writeouput('CircRNACount', res, samplelist, header=True)
            if options.annotate:
                logging.info('Write in annotation.')
                circAnn.annotate('tmp_coordinates',options.annotate,'CircCoordinates')
                circAnn.annotateregions('CircCoordinates',options.annotate)
            else:
                os.rename('tmp_coordinates','CircCoordinates')

        if not options.debug:
            deleted=cm.deletfile(os.getcwd(),circfiles+[files+'mapped' for files in circfiles])
            logdeleted(deleted)
        
    ### Filtering        
    if options.filter:
        logging.info('Program start to do filering')
        
        import circFilter as FT
        filt = FT.Circfilter(length=options.length, level0=options.Nr[0], level1=options.Nr[1], threshold0=options.Nr[2], threshold1=options.Nr[3])
        
        if not options.detect and len(options.Input) == 2: # Only use the program for filtering, need one coordinates file (bed6), one count file
            try:
                file2filter = options.Input[0]
                coorfile = options.Input[1]
                logging.info('Take file %s and %s for filtering' % (options.Input[0],options.Input[1]) )
                print 'Take file %s and %s for filtering' % (options.Input[0],options.Input[1])
            except IndexError:
                sys.exit('Please check the input. If only use the program for filtering, a coordinate file in bed6 format and a count file is needed.')
                logging.error('Program exit because input error. Please check the input. If only use the program for filtering, a coordinate file in bed6 format and a count file is needed.')
                
        elif options.detect:
            file2filter = 'tmp_circCount'
            coorfile = 'tmp_coordinates'
            logging.info('Take file tmp_circCount and tmp_coordinates for filtering') 
            print 'Take file tmp_circCount and tmp_coordinates for filtering'
            
        if options.rep_file:
            count,indx = filt.readcirc(file2filter,coorfile)
            logging.info('Filter by read counts.')
            count0,indx0 = filt.filtercount(count,indx) # result of first filtering by read counts
            filt.makeregion(indx0)
            logging.info('Filter by non repetitive region.')
            nonrep_left,nonrep_right = filt.nonrep_filter('tmp_left','tmp_right',options.rep_file)
            filt.intersectLeftandRightRegions(nonrep_left,nonrep_right,indx0,count0)
            if not options.chrM and not options.filterbygene:
                filt.sortOutput('tmp_unsortedWithChrM','CircRNACount',samplelist,'CircCoordinates',split=True)
        else:
            logging.error( 'A repetitive annotation file needed.' )
            sys.exit( 'A repetitive annotation file needed.' )
            
        # Filter chrM, if no further filtering, return 'CircRNACount' and 'CircCoordinates', else return 'tmp_unsortedNoChrM'
        if options.chrM:
            logging.info('Deleting circRNA candidates from Mitochondria chromosome.')
            filt.removeChrM('tmp_unsortedWithChrM')
            if not options.filterbygene:
                filt.sortOutput('tmp_unsortedNoChrM','CircRNACount',samplelist,'CircCoordinates',split=True)
        else:
            os.rename('tmp_unsortedWithChrM','tmp_unsortedNoChrM') # Note in this case 'tmp_unsortedNoChrM' actually has chrM
        
        # Filter by gene annotation, require one circRNA could not from more than one gene. return final 'CircRNACount' and 'CircCoordinates'
        if options.filterbygene:
            if options.annotate:
                logging.info('Filter by gene annotation. CircRNA candidates from more than one genes are deleted.')
                circAnn.filtbygene('tmp_unsortedNoChrM','tmp_unsortedfilterbygene')
                filt.sortOutput('tmp_unsortedfilterbygene','CircRNACount',samplelist,'CircCoordinates',split=True)
            else:
                logging.warning('To filter by gene annotation, a annotation file in bed 4 format needed, skiped filter by gene annotation.')
                filt.sortOutput('tmp_unsortedNoChrM','CircRNACount',samplelist,'CircCoordinates',split=True)
        
        # Add annotation of regions
        circAnn.annotateregions('CircCoordinates',options.annotate)
        
        logging.info('Filtering finished')
                      
        
    if options.gene:
        import genecount as GC
        # import the list of bamfile names as a file
        if not options.bam:
            print 'Please provide bam files, program will not count host gene expression.'
            logging.warning( 'Please provide bam files, program will not count host gene expression.' )
            
        if not options.refseq:
            print 'Please provide reference sequence, program will not count host gene expression.'
            logging.warning( 'Please provide reference sequence, program will not count host gene expression.' )
            
        if options.bam and options.refseq:
            # check whether the number of bamfiles is equale to the number of chimeric.junction.out files
            if len(options.bam) != len(options.Input):
                logging.error( "The number of bam files does not match with chimeric junction files." )
                sys.exit("The number of bam files does not match with chimeric junction files.")
            else:
                # For each sample (each bamfile), do one host gene count, and then combine to a single table
                gc = GC.Genecount()

                linearfiles = [] # A list for .linear file names
            
                for indx, files in enumerate(options.Input):
                    
                    if same:
                        linearfilename = getfilename(files)+str(indx)+'.linear'
                    else:
                        linearfilename = getfilename(files)+'.linear'
                    linearfiles.append(linearfilename)


                for indx, bamfile in enumerate(options.bam):
                    if options.circ:
                        logging.info('Counting linear gene expression based on provided circRNA coordinates')
                        print 'Counting linear gene expression based on provided circRNA coordinates'

                        gc.comb_gen_count(options.circ, bamfile, options.refseq, linearfiles[indx], countlinearsplicedreads=False)
                    else:
                        logging.info('Counting host gene expression based on detected and filtered circRNA coordinates')
                        print 'Counting host gene expression based on detected and filtered circRNA coordinates'
                        gc.comb_gen_count('CircRNACount', bamfile, options.refseq, linearfiles[indx], countlinearsplicedreads=False)
                
                logging.info("Finished linear gene expression counting, start to combine individual sample counts")
            
            # Combine all to a individual sample host gene count to a single table
            res = cm.combine(linearfiles,col=6,circ=False)
            cm.writeouput('LinearCount', res, samplelist, header=True)
            logging.info('Finished combine individual linear gene expression counts')
                
            if not options.debug:
                deleted=cm.deletfile(os.getcwd(),linearfilename)
                logdeleted(deleted)

    # Delte temp files
    if not options.debug:
        p3 = r'tmp\.*'
        deleted = cm.deletfile(os.getcwd(),p3)
        logdeleted(deleted)


def getfilename(namestring):
    tmp = namestring.split('/')
    filename = tmp[-1].strip('\n')
    return filename 

def logdeleted(deleted):
    for itm in deleted:
        logging.info('File'+' '+itm+' '+'Deleted!')

def mergefiles(files,output):
    # files is a list of file names
    with open(output,'w') as outfile:
        for fname in files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
        
        
if __name__ == "__main__":
    main()
