#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright© 2015 IBP of CHU of Grenoble
# DeCoRA.py is a free software : you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANYWARRANTY
# ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program. If not, see http://www.gnu.org/licenses/.
# For questions, please contact Alice CHOURY at alice.choury@hotmail.fr

# Note of Eméric BANKOLE (20170111)
# NB: For any debug don't forget to check always PATHWAY when switching between local computer or local server
# before launching DeCoRA

import Queue, threading, sys, globalvar, os, argparse, re
from sample import *
from html import *
import time


#----------------------------- GENERATING BED WITH INTRONICS REGIONS -------------------------------

def generating_bed_with_introns():
    IN = open(globalvar.BEDFILENAME,"r")
    lines = iter(IN.read().splitlines())
    output_filename = globalvar.OUTPUT_DIRECTORY+"/ROI.bed"
    OUT = open(output_filename,"w")
    for line in lines:
        infos = line.split("\t")
        start = int(infos[1])
        end = int(infos[2])
        label = infos[3]
        sens = label.split("(")[1][:1]
        start_output = start
        end_output = end
        if sens == "+":
            start_output = start - globalvar.LIMIT_INTRONIC_REGIONS_5P
            end_output = end + globalvar.LIMIT_INTRONIC_REGIONS_3P
        else:
            start_output = start - globalvar.LIMIT_INTRONIC_REGIONS_3P
            end_output = end + globalvar.LIMIT_INTRONIC_REGIONS_5P
        OUT.write(infos[0]+"\t"+str(start_output)+"\t"+str(end_output)+"\t"+label+"\n")
    OUT.close()
    IN.close()
    globalvar.BEDFILENAME = output_filename


#----------------------------- SEARCH FOR CNV -----------------------------

def mediane(L):
    L = [float(x) for x in L]
    L.sort()
    N = len(L)
    n = N/2.0
    p = int(n)
    if n == p:
        return (L[p-1]+L[p])/2.0
    else:
        return L[p]
        

#----------------------------- Add somme informations in each results file header (analysis_* and report_*)-------------

def header_info():
    
    pass
        
# ------------------------------------------ RESEARCH OF CNV  ------------------------------------------

# For all analysis of a run (analysis file in a directory)
# Calcul of median of depth of every analysis
# if a depth is twice the median, it's maybe a amplification
# if a depth is the half of the median, it's maybe a deletion
def CNV_hypothesis():
    if globalvar.VERBEUX:
        print 'Start of analysis for search CNV...'
    result_UNIX = os.popen('ls '+globalvar.OUTPUT_DIRECTORY+'/*_analysis.txt')
    filenames = result_UNIX.read().split('\n') 
    minimal_depths={}
    maximal_depths={}
    average_depths={}
    median_of_min_depths={} 
    median_of_max_depths={} 
    median_of_average_depths={} 
    values={} 
    for filename in filenames[:-1]:
        # sample represent barcode 
        file_sample = filename.split('/')[-1] 
        #sample = filename.split('/')[-1][9:-4]
        # modif line : 2 barcodes type (Example: IonDx-2 and IonXpress_002)
        if ("IonDx" in file_sample):
            sample = file_sample.split('_')[1]
        elif ("IonXpress" in file_sample):
            sample = file_sample.split('_')[1]+"_"+file_sample.split('_')[2]
        else:
            print "No barcode name, may be you are not on pgm server, please check it!"
            sample = "None"
        ind = Sample(sample)


        if globalvar.VERBEUX:
            print 'Reading of...'+sample+'_analysis.txt...'
        nb_reads = os.popen('samtools view '+globalvar.INPUT_DIRECTORY+'/'+sample+'_rawlib.bam | wc -l').read().split('\n')[0]  
        values[filename] = {}
        IN = open(filename, 'r')
        # We have 8 line in the header of results file
        line=IN.readline()
        line=IN.readline()
        line=IN.readline()
        line=IN.readline()
        line=IN.readline()
        line=IN.readline()
        line=IN.readline()
        line=IN.readline()
        line=IN.readline() # "\n exist here"
        line=IN.readline() # "\n exist here"

        while line != '':
            infos=line.split()
            min_depth_norm = float(infos[1]) / float(nb_reads)
            average_depth_norm = float(infos[2]) / float(nb_reads)
            max_depth_norm = float(infos[3]) / float(nb_reads)
            values[filename][infos[0]] = []
            values[filename][infos[0]].append(min_depth_norm)
            values[filename][infos[0]].append(average_depth_norm)
            values[filename][infos[0]].append(max_depth_norm)
            if not minimal_depths.has_key(infos[0]):
                minimal_depths[infos[0]] = []
                maximal_depths[infos[0]] = []
                average_depths[infos[0]] = []
            minimal_depths[infos[0]].append(min_depth_norm)
            maximal_depths[infos[0]].append(max_depth_norm)
            average_depths[infos[0]].append(average_depth_norm)
            line=IN.readline() 
        IN.close()

    for exon in minimal_depths:
        median_of_min_depths[exon]=mediane(minimal_depths[exon])
        median_of_max_depths[exon]=mediane(maximal_depths[exon])
        median_of_average_depths[exon]=mediane(average_depths[exon])
    RES = open(globalvar.OUTPUT_DIRECTORY+'CNV_predictions.txt', 'w')

    # Add header info in CNV file
    file_date = str(time.strftime('%Y-%m-%d %H:%M',time.localtime()))
    report_name = globalvar.OUTPUT_DIRECTORY.split("/")[-4] # Pathway on local computer
    #report_name = globalvar.OUTPUT_DIRECTORY.split("/")[-4] # Pathway for file on server pgm
    
    RunName = ""
    TeamName = ""

    with open(globalvar.INPUT_DIRECTORY+"/expMeta.dat") as data_file: # Pathway for expMeta.dat file (to check)
        lignes=data_file.readlines()
        RunName = lignes[0].split("=")[1]
        TeamProject = (lignes[3].split("=")[1])
        if re.findall("_\d+", TeamProject)[0]:
            m = re.findall("_\d+", TeamProject)[0]
            TeamName = (lignes[3].split("=")[1]).replace(m,"")
        sn2=re.findall("SN2-\d+", RunName)[0]
        if not re.search(sn2, report_name):
            print "No good report name or report name incomplete in the CNV analysis file header, may be you are not on pgm server, please check expMeta.dat file existing .."
            report_name = lignes[16].split("=")[1]
                
    # Write header info
    RES.write("Equipe : \t"+TeamName)
    RES.write("Run : \t"+RunName)
    RES.write("Report Name : \t"+report_name+"\n")
    RES.write("Decora version : \tv1.0.3"+"\n")
    RES.write("Decora CNV file date : \t"+file_date)
    RES.write("\n\n")

    if globalvar.VERBEUX:
        print 'Writting results...'
    for filename in values:
        for exon in values[filename]:
            if float(values[filename][exon][2]) > float(median_of_max_depths[exon]*2) and float(values[filename][exon][1]) > float(median_of_average_depths[exon]*3): 
                #sample = filename.split('_')[-1][:-4]
                # Modif line
                if ("IonDx-" in filename):
                    sample = filename.split('/')[-1].split('_')[0]+"\t"+filename.split('/')[-1].split('_')[1]
                elif ("IonXpress" in filename):
                    sample = filename.split('/')[-1].split('_')[0]+"\t"+filename.split('/')[-1].split('_')[1]+"_"+filename.split('/')[-1].split('_')[2]
                else:
                    print "No barcode name"
                    sample = "None"
                RES.write(sample+'\t'+exon+'\tamplifcation')
                RES.write('\n')

            if float(values[filename][exon][1]) < float(median_of_average_depths[exon]/2) and float(values[filename][exon][0])==0:
                #sample = filename.split('_')[-1][:-4]
                # Modif line
                if ("IonDx-" in filename):
                    sample = filename.split('/')[-1].split('_')[0]+"\t"+filename.split('/')[-1].split('_')[1]
                elif ("IonXpress" in filename):
                    sample = filename.split('/')[-1].split('_')[0]+"\t"+filename.split('/')[-1].split('_')[1]+"_"+filename.split('/')[-1].split('_')[2]
                else:
                    print "No barcode name"
                    sample = "None"
                RES.write(sample+'\t'+exon+'\tdeletion')
                RES.write('\n')
    RES.close()
    if globalvar.VERBEUX:
        print 'End of analysis for search CNV...'

# ------------------------------------------ Computing of no covered exons for all sample of a run  -------------------------------------

def Compute_no_covered_exon():
    if globalvar.VERBEUX:
        print 'Start of analysis for no covered exons...'
    #change path for debug
    #result_UNIX = os.popen('ls ../*_report.txt') # on my external drive usb
    result_UNIX = os.popen('ls '+globalvar.OUTPUT_DIRECTORY+'/*_report.txt') # on the pgm server
    filenames = result_UNIX.read().split('\n')
    num_patient = ""
    dict_values={}

    for filename in filenames[:-1]:
        num_patient = filename.split('/')[-1].split('_')[0]
        if globalvar.VERBEUX:
            print "starting for "+num_patient
        IN = open(filename, 'r')
        dict_values[num_patient]=[]
        list_tup = []
        lignes = IN.readlines()[9:]
        for line in lignes:
            colonnes = line.split("\t")
            list_exons = colonnes[5].split(",")
            gene = colonnes[0]
            if list_exons != ['']:
                # print line
                for exon in list_exons:
                    tup = (gene,exon)
                    if tup not in list_tup:
                        list_tup.append(tup)
                        dict_values[num_patient]=list_tup
        IN.close()
        print "end for "+num_patient

    #Just for debug
    #print dict_values

    all_gene_exon = []
    for list_gene_exon in dict_values.values():
        for gene_exon in list_gene_exon:
            #print gene_exon
            if gene_exon not in all_gene_exon:
                all_gene_exon.append(gene_exon)
    #Just for debug
    #print all_gene_exon
    
    res = ""
    for gene_exon in all_gene_exon:
        res+=gene_exon[0]+'\t'+gene_exon[1]+'\t'
        list_patient = []
        for patient in dict_values.keys():
            if (gene_exon in dict_values[patient]):
                if patient not in list_patient:
                    list_patient.append(patient)
        res+=" / ".join(list_patient)+'\t'+str(len(list_patient))+'\n'
    #Just for debug
    #print res

    # Add header info in CNV file
    file_date = str(time.strftime('%Y-%m-%d %H:%M',time.localtime()))
    report_name = globalvar.OUTPUT_DIRECTORY.split("/")[-4] # Pathway on local computer
    #report_name = globalvar.OUTPUT_DIRECTORY.split("/")[-4] # Pathway for file on server pgm
    
    RunName = ""
    TeamName = ""

    with open(globalvar.INPUT_DIRECTORY+"/expMeta.dat") as data_file: # Pathway for expMeta.dat file (to check)
        lignes=data_file.readlines()
        RunName = lignes[0].split("=")[1]
        TeamProject = (lignes[3].split("=")[1])
        if re.findall("_\d+", TeamProject)[0]:
            m = re.findall("_\d+", TeamProject)[0]
            TeamName = (lignes[3].split("=")[1]).replace(m,"")
        sn2=re.findall("SN2-\d+", RunName)[0]
        if not re.search(sn2, report_name):
            print "No good report name or report name incomplete in the CNV analysis file header, may be you are not on pgm server, please check expMeta.dat file existing .."
            report_name = lignes[16].split("=")[1]
    
    RES = open(globalvar.OUTPUT_DIRECTORY+'No_covered_exons_report.txt', 'w')       
    # Write header info
    RES.write("Equipe : \t"+TeamName)
    RES.write("Run : \t"+RunName)
    RES.write("Report Name : \t"+report_name+"\n")
    RES.write("Decora version : \tv1.0.3"+"\n")
    RES.write("Decora no covered exons file date : \t"+file_date)
    RES.write("\n\n")

    header = "Gene"+"\t"+"Exon"+"\t"+"Echantillon"+"\t"+"Nombre d'echantillons"+"\n"
    if globalvar.VERBEUX:
        print "writting results of no covered exons"
    
    RES.write(header)
    RES.write(res)
    RES.close()
    if globalvar.VERBEUX:
        print "End of Computing for no covoregae exons"



# --------------------------------------------- ANALYSIS ----------------------------------------------

# Function for one analysis
#
# Using SAMTools to have depth base por base
# options of SAMTools : 
# -q threshold of base quality
# -Q threshold of mapping quality
# here it's fixed at -q 6 -Q 4 because it's the threshold of Life for the VariantCaller 
#
# search the minimal depth and the maximal bias
# select the exon with a minimal depth under the depth threshold and/or with a maximal biais upper the bias threshold
#
# can analyse depths for a list of mutations (hotspots)
# select hotspot with a minimal depth under the depth threshold
#
# write a report of quality of sequencing and brut results
# write a wiggle file for IGV (if the option is used)
def analysis(path): 
    filename = path.split('/')[-1]
    sample = filename[:-11]    
    if sample != '':
        results = Sample(sample)  
        if globalvar.VERBEUX:
            print 'Start of analysis for '+sample+'...\n'
            print 'Start of running SAMTools for '+sample+'...\n'
        os.system('samtools depth -b '+globalvar.BEDFILENAME+' -q 6 -Q 4 '+globalvar.INPUT_DIRECTORY+'/'+ filename+' > '+globalvar.OUTPUT_DIRECTORY+'/tmp/coverage'+sample)
        if globalvar.VERBEUX:
            print 'End of running SAMTools for '+sample+'...\n'
        results.data_exons_calculation()
        results.select_exons_low_depth()
        results.data_genes_calculation()
        if globalvar.VERBEUX:
            print 'Start of bias analysis for '+sample+'...\n'
        os.system('samtools view -F 16 -b '+globalvar.INPUT_DIRECTORY+'/'+filename+' > '+globalvar.OUTPUT_DIRECTORY+'/tmp/'+sample+'_forward.bam')
        os.system('samtools depth -b '+globalvar.BEDFILENAME+' -q 6 -Q 4 '+globalvar.OUTPUT_DIRECTORY+'/tmp/'+sample+'_forward.bam > '+globalvar.OUTPUT_DIRECTORY+'/tmp/bias_forward'+sample)
        os.system('samtools view -f 16 -b '+globalvar.INPUT_DIRECTORY+'/'+filename+' > '+globalvar.OUTPUT_DIRECTORY+'/tmp/'+sample+'_reverse.bam')
        os.system('samtools depth -b '+globalvar.BEDFILENAME+' -q 6 -Q 4 '+globalvar.OUTPUT_DIRECTORY+'/tmp/'+sample+'_reverse.bam > '+globalvar.OUTPUT_DIRECTORY+'/tmp/bias_reverse'+sample)
        results.bias_calculation()
        if globalvar.BIAS_SELECTION:
            results.select_exons_high_bias()
        if globalvar.VERBEUX:
            print 'End of bias analysis for '+sample+'...\n'
        if not globalvar.HOTSPOT_FILE == '':
            if globalvar.VERBEUX:
                print 'Start of analysis of Hotspot for '+sample+'...\n'
            results.add_hotspots_of_mutationsfile()
            results.data_hotspots_calculation(sample)
            results.select_hotspots_low_coverage()
            results.save_hotspot_analysis()
        if globalvar.VERBEUX:
            print 'Saving for '+sample+'...\n'
        results.save_analysis()
        results.save_report()
        if globalvar.BED2WIG:
            results.bed2wig()
        if globalvar.VERBEUX:
            print 'End of analysis for '+sample+'...\n'
        
# ------------------- CLASSES FOR MULTI-THREADING ------------------

queueLock = threading.Lock()
workQueue = Queue.Queue(0)
threads = []

class Analysis (threading.Thread):
    def __init__(self, threadID, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.q = q

    def run(self):        
        process_analysis(self.q)

def process_analysis(q):
    while not globalvar.exitFlag:       
        queueLock.acquire()
        if not workQueue.empty():
            filename = q.get()
            queueLock.release()
            analysis(filename)
        else:
            queueLock.release()

# -------------------------------------------------- MAIN ---------------------------------------------------

def main():

# -------------------------------------------------- OPTIONS --------------------------------------------------
    parser = argparse.ArgumentParser(description='DeCoRA')    
    parser.add_argument('-i', '--input_directory', help='directory with BAM files to analysis',default='test')
    parser.add_argument('-r', '--roi_file', help='ROI file in BED format',default='test')
    parser.add_argument('-m', '--mutations_file', help='mutations file in BED format',default='')
    parser.add_argument('-o', '--output_directory', help='directory for results of analysis',default='results')  
    parser.add_argument('-w', '--wiggle_file', help='WIGGLE files for IGV',action='store_true',default=False) 
    parser.add_argument('-d', '--depth_threshold', help='depth threshold, by default: 30',type=int,default=30)
    parser.add_argument('-b', '--bias_threshold', help='bias threshold, by default : 0.2',type=float,default=0.2)    
    parser.add_argument('-B', '--bias_selection', help='show in the report exons with high bias',action='store_true',default=False)
    parser.add_argument('--CNV', help='deletion and amplification prediction',action='store_true',default=False)
    parser.add_argument('--CNVO', help='deletion and amplification prediction without depth analysis',action='store_true',default=False)
    parser.add_argument('-a', '--absolute', help='absolute minimum absolute average in report',action='store_true',default=False)
    parser.add_argument('--intronics_regions',nargs=3,help='in the BED file, presence of 10 bases of intronics regions in 5\' and 3\' of the ROI',default=["f",10,10])
    parser.add_argument('-v', '--verbose', help='debug version',action='store_true',default=False)
    args = parser.parse_args() 
    globalvar.INPUT_DIRECTORY =args.input_directory
    globalvar.OUTPUT_DIRECTORY =args.output_directory
    globalvar.DEPTH_THRESHOLD = args.depth_threshold
    globalvar.BIAS_THRESHOLD = args.bias_threshold
    globalvar.HOTSPOT_FILE = args.mutations_file    
    globalvar.BIAS_SELECTION = args.bias_selection
    globalvar.CNV_ANALYSIS = args.CNV
    globalvar.BEDFILENAME = args.roi_file
    globalvar.ASOLUTE_MIN = args.absolute
    globalvar.INTRONICS_REGIONS_ANALYSIS = args.intronics_regions[0][0] == "t"
    if (globalvar.INTRONICS_REGIONS_ANALYSIS):
        # modif line amandine
        globalvar.LIMIT_INTRONIC_REGIONS_5P = int(args.intronics_regions[1])
        globalvar.LIMIT_INTRONIC_REGIONS_3P = int(args.intronics_regions[2])
    globalvar.CNV_ANALYSIS_ONLY = args.CNVO
    globalvar.VERBEUX = args.verbose
    globalvar.BED2WIG = args.wiggle_file 

    if globalvar.CNV_ANALYSIS_ONLY:
        CNV_hypothesis()        
    else:    
        if globalvar.VERBEUX:
            print '\nStarting DeCoRA...\n'
            print '\nEnvironment variables used for the analysis :'
            print '\t- INPUT : '+globalvar.INPUT_DIRECTORY
            print '\t- OUTPUT : '+globalvar.OUTPUT_DIRECTORY
            print '\t- BED FILE NAME : '+globalvar.BEDFILENAME
            print '\t- DEPTH THRESHOLD : '+str(globalvar.DEPTH_THRESHOLD)
            if globalvar.BIAS_SELECTION:
                print '\t- BIAS THRESHOLD : '+str(globalvar.BIAS_THRESHOLD)
            if globalvar.HOTSPOT_FILE != "":
                print '\t- HOTSPOTS FILE NAME : '+globalvar.HOTSPOT_FILE
            rep = "no"
            if globalvar.CNV_ANALYSIS:  rep="yes"
            print '\t- CNV ANALYSIS : '+rep
            rep = "no"
            if globalvar.ASOLUTE_MIN:  rep="yes"
            print '\t- ABSOLUTE MIN IN REPORT : '+rep
            rep = "no"
            if globalvar.INTRONICS_REGIONS_ANALYSIS:  rep="yes"
            print '\t- PRESENCE OF INTRONIC REGIONS : '+rep
            rep = "no"
            if globalvar.BED2WIG:  rep="yes"
            print '\t- WIGGLE FILE : '+rep
            print "\n\n"
# -------------------------------------------------- INPUT --------------------------------------------------
        result_UNIX = os.popen(' ls '+globalvar.INPUT_DIRECTORY+'/*.bam')
        filenames = result_UNIX.read().split('\n')          
        os.system('mkdir '+globalvar.OUTPUT_DIRECTORY+'/tmp/')
        if (globalvar.INTRONICS_REGIONS_ANALYSIS):
            generating_bed_with_introns()
            
# -------------------------------------------------- FOR TMPFS --------------------------------------------------
        #os.system('sudo mount -t tmpfs -o size=2g tmpfs /results/plugins/DeCoRA/python/tmp')
        
# -------------------------------------------------- FOR MULTI-THREADING --------------------------------------------------
        

        # Create new threads
        for threadID in range (1,globalvar.NB_TREADS+1):
            thread = Analysis(threadID, workQueue)
            thread.start()
            threads.append(thread)
        # Fill the queue
        queueLock.acquire()
        for filename in filenames[:-1]:
            if filename != 'rawlib.bam':
                workQueue.put(filename)
        queueLock.release()
        # Wait for queue to empty
        while not workQueue.empty():
            pass
        # Notify threads it's time to exit
        globalvar.exitFlag = 1
        # Wait for all threads to complete
        for t in threads:
            t.join()
            
        if globalvar.CNV_ANALYSIS:
            CNV_hypothesis()

        Compute_no_covered_exon()    

        write_html()       
        os.system('rm -rf '+globalvar.OUTPUT_DIRECTORY+'/tmp/')
        #os.system('sudo umount /results/plugins/DeCoRA/python/tmp')
            
# Call main()
if __name__ == '__main__':
    main()
