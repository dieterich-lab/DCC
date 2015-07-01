#! /usr/bin/env python
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
from fix2chimera import Fix2Chimera
#import pdb

def main():
    
    parser = argparse.ArgumentParser(prog='DCC',formatter_class=argparse.RawDescriptionHelpFormatter, fromfile_prefix_chars='@', description='Contact jun.cheng@age.mpg.de')
    
    
    parser.add_argument('--version', action='version', version='%(prog)s 0.3.1c')
    parser.add_argument("Input", metavar='Input', nargs="+",
                    help="Input of the chimeric.out.junction file from STAR. Alternatively, a sample sheet specifying where your chimeric.out.junction files are, each sample per line, provide with @ prefix (e.g. @samplesheet).")
    #parser.add_argument("-O", "--output", dest="output", 
    #                  help="Tab delimited outputfile, order the same with input: \
    #                  chr\tstart\tend\tstand\tcount\tjunctiontype")
    parser.add_argument("-temp", "--temp", dest="temp", action='store_true', default=False,
                    help="Once specified, temp files will not be deleted.")
                                                                                                        
    group = parser.add_argument_group("Find circRNA Options","Options to find circRNAs from STAR output.")
    group.add_argument("-D", "--detect", action='store_true', dest="detect", default=False,
                    help="Always specify if you want detect circRNAs from chimeric junctions.")                    
    #group.add_argument("-S", action='store_true', dest="strand", default=True,
    #                help="Specify when the library is stranded [Default].")
    group.add_argument("-N", "--nonstrand", action='store_false', dest="strand", default=True,
                    help="Specify when the library is non-stranded [default stranded].")
    group.add_argument("-E", "--endTol", dest="endTol", type=int, default= 5, choices=range(0,10),
                    help="Maximum base pair tolerance of reads extending over junction sites. [Interger, default 5]")
    group.add_argument("-m", "--maximum", dest="max", type=int, default = 1000000,
                    help="The maximum range of candidate circRNA allowed (including introns). [default 1000000]")   
    group.add_argument("-n", "--minimum", dest="min", type=int, default = 30,
                    help="The minimum range of candidate circRNA allowed (including introns). [default 30]")
    group.add_argument("-an", "--annotation",dest="annotate",
                    help="Gene annotation file in GTF/GFF3 format, to annotate circRNAs by their host gene name/identifier.")
    #group.add_argument("-gf", "--getfasta", dest="getfasta",
    #                  help="Get fasta file of circular RNAs. If a exon annotation file is provided, the circular RNA sequence will only contain annotated exons, otherwise whole sequence.")
    group.add_argument("-Pi", "--PE-independent", action='store_true', dest="pairedendindependent", default=False,
                    help="Specify when you have maped the PE data mates seperately. If specified, -mt1 and -mt2 should also be provied. [default False]")
    group.add_argument("-mt1", "--mate1", dest="mate1", nargs = '+',
                    help="For paired end data, Chimeric.out.juntion files from mate1 independent mapping result.")  
    group.add_argument("-mt2", "--mate2", dest="mate2", nargs = '+',
                    help="For paired end data, Chimeric.out.juntion files from mate2 independent mapping result.")
    parser.add_argument_group(group)
    
        
    group = parser.add_argument_group("Filtering Options", "Options to filter the circRNA candidates.")
    group.add_argument("-F", "--filter", action='store_true', dest="filter", default=False,
                    help="If specified, the program will do filtering on the detection results.")
    group.add_argument("-M", "--chrM", action='store_true', dest="chrM",default=False,
                    help="If specified, candidates from mitochondria chromosome will be removed.")
    #group.add_argument("-J", "--junction", dest="junction",
    #                  help="Provide a coustom junction file in gtf format, if only specify as True, only GT/AG or CT/AC junction will be considered.")
    group.add_argument("-R", "--rep_file", dest="rep_file",
                    help="Coustom repetitive region file in GTF format, to filter out circRNAs candidates in repetitive regions.") 
    group.add_argument("-L", "--Ln", dest="length", type=int, default=50,
                    help="Minimum length to check for repetitive regions. [default 50]")                                                      
    group.add_argument('-Nr', nargs=2, type=int, metavar=('level1','threshold1'), default=[2,5], help='Minimum read counts required for circRNAs; \
                        Minimum number of samples above the corresponding expression level')
    group.add_argument("-fg", "--filterbygene", action='store_true', dest="filterbygene", default=False,
                    help="If specified, filter by gene annotation. candidates are not allowed to span more than one gene.")                                        
    parser.add_argument_group(group)
    
    
    group = parser.add_argument_group("Host gene count Options", "Options to count host gene expression.")
    group.add_argument("-G", "--gene", action='store_true', dest="gene", default=False,
                    help="If specified, the program will count host gene expression given circRNA coordinates. By default, use the circRNA candidates detected from the same run.")
    group.add_argument("-C", "--circ", dest="circ",
                    help="User specified circRNA coordinates, any tab delimited file with first three columns as circRNA coordinates: chr\tstart\tend, which DCC will use to count host gene expression.")                  
    group.add_argument("-B", "--bam", dest="bam", nargs = '+',
                    help="A file specify where the mapped bam files are, from which host gene expression are counted from. Should be the same oder with the input chimeric junction files.")
    group.add_argument("-A", "--refseq", dest="refseq",
                    help="Reference sequnce fasta file.")
    #group.add_argument("-seq", "--seq", dest="seq",         
    #                  help="Get circRNA sequence as fasta file.")
    parser.add_argument_group(group)
                                        
    options= parser.parse_args()
    
    timestr = time.strftime("%Y%m%d-%H%M%S")
    logging.basicConfig(filename='DCC.log'+timestr, filemode='w',level=logging.DEBUG,format='%(asctime)s %(message)s')
    
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
    
    try:
        os.mkdir('_tmp_DCC')
    except OSError:
        from shutil import rmtree
        rmtree('_tmp_DCC/')
        os.mkdir('_tmp_DCC')
    
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
        logging.info('Program start to detect circRNAs.')
        if options.strand:
            logging.info('Strand data.')
        else:
            logging.info('nonstrand data, the strand of circRNAs guessed from the strandness of host genes.')
            print 'WARNING: nonstrand data, the strand of circRNAs guessed from the strandness of host genes.'
        
        # Start de novo circular RNA detection model
        # Create instances
        f = FC.Findcirc(endTol=options.endTol,maxL=options.max,minL=options.min)
        sort = FC.Sort()
        
        circfiles = [] # A list for .circRNA file names      
        
        def wrapfindcirc(files,output,strand=True,pairdendindependent=True):
            if pairdendindependent:
                f.printcircline(files,'_tmp_DCC/tmp_printcirclines')
                f.sepDuplicates('_tmp_DCC/tmp_printcirclines','_tmp_DCC/tmp_duplicates','_tmp_DCC/tmp_nonduplicates')
                # Find small circles
                f.smallcirc('_tmp_DCC/tmp_duplicates','_tmp_DCC/tmp_smallcircs')
                if strand:
                    # Find normal circles
                    f.findcirc('_tmp_DCC/tmp_nonduplicates','_tmp_DCC/tmp_normalcircs',strand=True)
                else:
                    f.findcirc('_tmp_DCC/tmp_nonduplicates','_tmp_DCC/tmp_normalcircs',strand=False)
                # Merge small and normal circles
                mergefiles('_tmp_DCC/tmp_findcirc','_tmp_DCC/tmp_smallcircs','_tmp_DCC/tmp_normalcircs')
            else:
                if strand:
                    f.findcirc(files,'_tmp_DCC/tmp_findcirc',strand=True)
                else:
                    f.findcirc(files,'_tmp_DCC/tmp_findcirc',strand=False)
            # Sort
            if strand:
                sort.sort_count('_tmp_DCC/tmp_findcirc',output,strand=True)
            else:
                sort.sort_count('_tmp_DCC/tmp_findcirc',output,strand=False)
        
        if options.pairedendindependent:
            print '===== Please make sure that you mapped both the paired mates togethor and seperately!!! ====='
            logging.info("===== Please make sure that you mapped both the paired mates togethor and seperately!!! =====")
  
            # Fix2chimera problem by STAR
            print ('Collect chimera from mates-seperate mapping.')
            logging.info('Collect chimera from mates-seperate mapping.')
            Input = fixall(options.Input,options.mate1,options.mate2)
        else:
            Input = options.Input

        for indx, files in enumerate(Input):
            logging.info( 'started detect circRNA from %s' % files )
            print 'started detect circRNA from %s' % files
            if same:
                circfilename = '_tmp_DCC/'+getfilename(files)+str(indx)+'.circRNA'
            else:
                circfilename = '_tmp_DCC/'+getfilename(files)+'.circRNA'
            circfiles.append(circfilename)
    
                                                    
            if options.strand:
                if options.pairedendindependent:
                    wrapfindcirc(files,circfilename,strand=True,pairdendindependent=True)
                else:
                    wrapfindcirc(files,circfilename,strand=True,pairdendindependent=False)

            else:      	        
                if options.pairedendindependent:
                    wrapfindcirc(files,circfilename,strand=False,pairdendindependent=True)
                else:
                    wrapfindcirc(files,circfilename,strand=False,pairdendindependent=False)
        #        
        #try:
        #    os.remove('_tmp_DCC/tmp_findcirc')
        #    os.remove('_tmp_DCC/tmp_printcirclines')
        #    os.remove('_tmp_DCC/tmp_duplicates')
        #    os.remove('_tmp_DCC/tmp_nonduplicates')
        #    os.remove('_tmp_DCC/tmp_smallcircs')
        #except OSError:
        #    pass
            
            
        ### Combine the individual count files
        # Create a list of '.circRNA' file names
        print('Start to combine individual circRNA read counts.')
        logging.info('Start to combine individual circRNA read counts.')

        cm.comb_coor(circfiles,strand=options.strand)
        cm.map('_tmp_DCC/tmp_coordinates', circfiles, strand=options.strand)

        res = cm.combine([files+'mapped' for files in circfiles],col=7,circ=True)

        if options.filter:
            cm.writeouput('_tmp_DCC/tmp_circCount', res)
            if options.annotate:
                logging.info('Write in annotation.')
                logging.info('Select gene features in Annotation file.')
                circAnn.selectGeneGtf(options.annotate)
                circAnn.annotate('_tmp_DCC/tmp_coordinates','_tmp_DCC/tmp_'+getfilename(options.annotate)+'.gene','_tmp_DCC/tmp_coordinatesannotated')
                os.remove('_tmp_DCC/tmp_coordinates')
                os.rename('_tmp_DCC/tmp_coordinatesannotated','_tmp_DCC/tmp_coordinates')
        else:
            cm.writeouput('CircRNACount', res, samplelist, header=True)
            if options.annotate:
                logging.info('Write in annotation.')
                logging.info('Select gene features in Annotation file.')
                circAnn.selectGeneGtf(options.annotate)
                circAnn.annotate('_tmp_DCC/tmp_coordinates','_tmp_DCC/tmp_'+getfilename(options.annotate)+'.gene','CircCoordinates')
                circAnn.annotateregions('CircCoordinates',options.annotate)
            else:
                os.rename('_tmp_DCC/tmp_coordinates','CircCoordinates')
        
    ### Filtering        
    if options.filter:
        logging.info('Program start to do filering')
        
        import circFilter as FT
        filt = FT.Circfilter(length=options.length,level1=options.Nr[0], threshold1=options.Nr[1])
        
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
            file2filter = '_tmp_DCC/tmp_circCount'
            coorfile = '_tmp_DCC/tmp_coordinates'
            logging.info('Take file _tmp_DCC/tmp_circCount and _tmp_DCC/tmp_coordinates for filtering') 
            print 'Take file _tmp_DCC/tmp_circCount and _tmp_DCC/tmp_coordinates for filtering'
            
        if options.rep_file:
            rep_file = options.rep_file
        else:
            from pkg_resources import resource_filename          
            rep_file = resource_filename('DCC', 'data/DCC.Repeats')
        count,indx = filt.readcirc(file2filter,coorfile)
        logging.info('Filter by read counts.')
        count0,indx0 = filt.filtercount(count,indx) # result of first filtering by read counts
        filt.makeregion(indx0)
        logging.info('Filter by non repetitive region.')
        nonrep_left,nonrep_right = filt.nonrep_filter('_tmp_DCC/tmp_left','_tmp_DCC/tmp_right',rep_file)
        filt.intersectLeftandRightRegions(nonrep_left,nonrep_right,indx0,count0)
        if not options.chrM and not options.filterbygene:
            filt.sortOutput('_tmp_DCC/tmp_unsortedWithChrM','CircRNACount',samplelist,'CircCoordinates',split=True)
            
        # Filter chrM, if no further filtering, return 'CircRNACount' and 'CircCoordinates', else return '_tmp_DCC/tmp_unsortedNoChrM'
        if options.chrM:
            logging.info('Deleting circRNA candidates from Mitochondria chromosome.')
            filt.removeChrM('_tmp_DCC/tmp_unsortedWithChrM')
            if not options.filterbygene:
                filt.sortOutput('_tmp_DCC/tmp_unsortedNoChrM','CircRNACount',samplelist,'CircCoordinates',split=True)
        else:
            os.rename('_tmp_DCC/tmp_unsortedWithChrM','_tmp_DCC/tmp_unsortedNoChrM') # Note in this case '_tmp_DCC/tmp_unsortedNoChrM' actually has chrM
        
        # Filter by gene annotation, require one circRNA could not from more than one gene. return final 'CircRNACount' and 'CircCoordinates'
        if options.filterbygene:
            if options.annotate:
                logging.info('Filter by gene annotation. CircRNA candidates from more than one genes are deleted.')
                circAnn.filtbygene('_tmp_DCC/tmp_unsortedNoChrM','_tmp_DCC/tmp_unsortedfilterbygene')
                filt.sortOutput('_tmp_DCC/tmp_unsortedfilterbygene','CircRNACount',samplelist,'CircCoordinates',split=True)
            else:
                logging.warning('To filter by gene annotation, a annotation file in GTF/GFF format needed, skiped filter by gene annotation.')
                filt.sortOutput('_tmp_DCC/tmp_unsortedNoChrM','CircRNACount',samplelist,'CircCoordinates',split=True)
        
        # Add annotation of regions
        if options.annotate:
            circAnn.annotateregions('CircCoordinates',options.annotate)
        
        logging.info('Filtering finished')
                      
        
    if options.gene:
        import genecount as GC
        # import the list of bamfile names as a file
        if not options.bam:
            #print 'Please provide bam files, program will not count host gene expression.'
            logging.info( 'Look for mapped bam files in the same directory as chimeric.out.junction files.' )
            bamfiles = convertjunctionfile2bamfile(options.Input)
        else:
            bamfiles = options.bam
            
        if not options.refseq:
            print 'Please provide reference sequence, program will not count host gene expression.'
            logging.warning('Please provide reference sequence, program will not count host gene expression.')
            
            
        if options.refseq:
            # check whether the number of bamfiles is equale to the number of chimeric.junction.out files
            if len(bamfiles) != len(options.Input):
                logging.error( "The number of bam files does not match with chimeric junction files." )
                sys.exit("The number of bam files does not match with chimeric junction files.")
            else:
                # For each sample (each bamfile), do one host gene count, and then combine to a single table
                gc = GC.Genecount()

                linearfiles = [] # A list for .linear file names
            
                for indx, files in enumerate(options.Input):
                    
                    if same:
                        linearfilename = '_tmp_DCC/'+getfilename(files)+str(indx)+'.linear'
                    else:
                        linearfilename = '_tmp_DCC/'+getfilename(files)+'.linear'
                    linearfiles.append(linearfilename)


                for indx, bamfile in enumerate(bamfiles):
                    if options.circ:
                        logging.info('Counting linear gene expression based on provided circRNA coordinates for %s' %bamfile)
                        #print 'Counting linear gene expression based on provided circRNA coordinates'

                        gc.comb_gen_count(options.circ, bamfile, options.refseq, linearfiles[indx], countlinearsplicedreads=False)
                    else:
                        logging.info('Counting host gene expression based on detected and filtered circRNA coordinates for %s' %bamfile)
                        print 'Counting host gene expression based on detected and filtered circRNA coordinates'
                        gc.comb_gen_count('CircRNACount', bamfile, options.refseq, linearfiles[indx], countlinearsplicedreads=False)
                
                logging.info("Finished linear gene expression counting, start to combine individual sample counts")
            
            # Combine all to a individual sample host gene count to a single table
            res = cm.combine(linearfiles,col=6,circ=False)
            cm.writeouput('LinearCount', res, samplelist, header=True)
            logging.info('Finished combine individual linear gene expression counts')
                
            if not options.temp:
                deleted=cm.deletfile(os.getcwd(),linearfiles)
                logdeleted(deleted)
    
    # CircSkip junction
    if options.annotate and not options.circ:
        logging.info('Count CircSkip junctions.')
        print('Count CircSkip junctions.')
        SJ_out_tab = getSJ_out_tab(options.Input)
        CircSkipfiles = findCircSkipJunction('CircCoordinates',options.annotate,circfiles,SJ_out_tab,strand=options.strand,same=same)
        fin = open('CircCoordinates','r').readlines()[1:]
        with open('_tmp_DCC/tmp_CircCoordinatesNoheader','w') as fout:
            fout.writelines(fin)
        cm.map('_tmp_DCC/tmp_CircCoordinatesNoheader',CircSkipfiles,strand=options.strand,col=4)
        CircSkipfilesmapped = [ fname+'mapped' for fname in CircSkipfiles ]
        res = cm.combine(CircSkipfilesmapped, col = 9)
        cm.writeouput('CircSkipJunctions',res,samplelist,header=True)
    else:
        logging.info('CircSkip junctions cannot be count.')
        
    # Delte temp files
    if not options.temp:
        p3 = r'^tmp_\.*'
        deleted = cm.deletfile(os.path.join(os.getcwd(),'_tmp_DCC/'),p3)
        logdeleted(deleted)
        deleted=cm.deletfile(os.getcwd(),circfiles+[files+'mapped' for files in circfiles])
        logdeleted(deleted)
        deleted=cm.deletfile('',CircSkipfiles)
        logdeleted(deleted)
        deleted=cm.deletfile('',CircSkipfilesmapped)
        logdeleted(deleted)    
    logging.info('Finished!')
        

def fixall(joinedfnames,mate1filenames,mate2filenames):
    # Fix all 2chimera in one read/read paire for all inputs
    # outputs as a list of fixed filenames, with .fixed end
    outputs = []
    fx = Fix2Chimera()
    # check mate1 and mate2 input
    if len(mate1filenames) == len(mate2filenames) == len(joinedfnames):
        for i in range(len(joinedfnames)):
            fx.fixation(mate1filenames[i],mate2filenames[i],joinedfnames[i],joinedfnames[i]+'.fixed')
            outputs.append(joinedfnames[i]+'.fixed')
    else:
        logging.error('The number of input mate1, mate2 and joined mapping files are different.')
        sys.exit('The number of input mate1, mate2 and joined mapping files are different.')

    return outputs
    

def getfilename(namestring):
    tmp = namestring.split('/')
    filename = tmp[-1].strip('\n')
    return filename 

def logdeleted(deleted):
    for itm in deleted:
        logging.info('File'+' '+itm+' '+'Deleted!')

def mergefiles(output,*fnames):
    import shutil
    destination = open(output,'wb')
    for fname in fnames:
        shutil.copyfileobj(open(fname,'rb'), destination)
    destination.close()        

def convertjunctionfile2bamfile(junctionfilelist):
    
    def getbamfname(junctionfname):
        import re
        import os
        # Get the stored directory
        dirt = '/'.join((junctionfname.split('/')[:-1]))+'/'
        p = r'.*Aligned\..*bam'
        for fname in os.listdir(dirt):
            if re.match(p,fname):
                bamfname = dirt+re.findall(p,fname)[0]
        return bamfname
        
    bamfnames = []       
    for fname in junctionfilelist:
        bamfnames.append(getbamfname(fname))
    return bamfnames 

# CircSkip junctions
def findCircSkipJunction(CircCoordinates,gtffile,circfiles,SJ_out_tab,strand=True,same=False):
    from Circ_nonCirc_Exon_Match import CircNonCircExon
    CircSkipfiles = []
    CCEM=CircNonCircExon()
    # Modify gtf file
    CCEM.select_exon(gtffile)
    if CCEM.modifyExon_id('_tmp_DCC/tmp_'+getfilename(gtffile)+'.exon.sorted'):
        # Start and end coordinates
        start2end=CCEM.print_start_end_file(CircCoordinates)
        Iv2Custom_exon_id, Custom_exon_id2Iv, Custom_exon_id2Length = CCEM.readNonUniqgtf('_tmp_DCC/tmp_'+getfilename(gtffile)+'.exon.sorted.modified')
        if strand:
            circStartExons = CCEM.intersectcirc('_tmp_DCC/tmp_start.bed','_tmp_DCC/tmp_'+getfilename(gtffile)+'.exon.sorted.modified') #Circle start or end to corresponding exons
        else:
            circStartExons = CCEM.intersectcirc('_tmp_DCC/tmp_start.bed','_tmp_DCC/tmp_'+getfilename(gtffile)+'.exon.sorted.modified',strand=False)
        circStartAdjacentExons, circStartAdjacentExonsIv = CCEM.findcircAdjacent(circStartExons,Custom_exon_id2Iv,Iv2Custom_exon_id,start=True)
        if strand:
            circEndExons = CCEM.intersectcirc('_tmp_DCC/tmp_end.bed','_tmp_DCC/tmp_'+getfilename(gtffile)+'.exon.sorted.modified') #Circle start or end to corresponding exons
        else:
            circEndExons = CCEM.intersectcirc('_tmp_DCC/tmp_end.bed','_tmp_DCC/tmp_'+getfilename(gtffile)+'.exon.sorted.modified',strand=False)
        circEndAdjacentExons, circEndAdjacentExonsIv = CCEM.findcircAdjacent(circEndExons,Custom_exon_id2Iv,Iv2Custom_exon_id,start=False)
        exonskipjunctions = CCEM.exonskipjunction(circStartAdjacentExonsIv,circEndAdjacentExonsIv,start2end)
        for indx, fname in enumerate(SJ_out_tab):
            if same:
                path = '_tmp_DCC/'+getfilename(fname).replace('SJ.out.tab',str(indx))
            else:
                path = '_tmp_DCC/'+getfilename(fname).replace('SJ.out.tab','')
            junctionReadCount = CCEM.readSJ_out_tab(fname)
            if len(junctionReadCount) == 0:
                logging.error('Do you have SJ.out.tab files in your sample folder? DCC cannot find it.')
                logging.info('Cannot fine SJ.out.tab files, please check the path. circSkip will not be output.')
                break
            else:
                skipJctCount = CCEM.getskipjunctionCount(exonskipjunctions,junctionReadCount)
                circCount=CCEM.readcircCount(circfiles[indx])
                CircSkipfile = CCEM.printCirc_Skip_Count(circCount,skipJctCount,path)
                CircSkipfiles.append(CircSkipfile)
        return CircSkipfiles

def getSJ_out_tab(chimeralist):
    SJ_out_tab = []
    for fname in chimeralist:
        SJ_out_tab.append(fname.replace('Chimeric.out.junction','SJ.out.tab'))
    return SJ_out_tab
    
    
if __name__ == "__main__":
    main()
