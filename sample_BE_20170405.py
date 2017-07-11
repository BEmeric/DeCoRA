#!/usr/bin/python
# -*- coding: utf-8 -*-


# Copyright© 2015 Alice CHOURY
# sample.py is a part of a free software, DeCoRA : you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANYWARRANTY
# ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program. If not, see http://www.gnu.org/licenses/.


import math, globalvar, re, json, sys, os
from gene import *
import time

reload(sys)
sys.setdefaultencoding("utf-8")

class Sample:   

# ------------------------------------------- CONSTRUCTOR --------------------------------------------
    def __init__(self,sample):
        self.name = sample
        # Modif Eméric
        self.num_patient = "None"
        self.ReportName = "None"
        self.RunName = "None"
        self.TeamName = "None"


        self.no_covered_exons = []
        self.no_covered_hotspots = []
        self.genes = []
        IN = open(globalvar.BEDFILENAME, "r")
        line=IN.readline()
        
        while line != '':
            infos=line.split()
            num_chr = infos[0]
            start = int(infos[1])
            end = int (infos[2])
            if start > end:
                tmp = start 
                start = end
                end = tmp
            label = infos[3]
            names = label.split("_")
            # Modif EB
            # check this line below by bedfile format (Pgm life or Illumina or anoteher data sequencing)
            gene_name = names[0][:-3]

            if not self.present(gene_name):  
                gene = Gene(gene_name)
                self.add_gene(gene)
            # Modif EB
            # # check this line below by bedfile format (Pgm life or Illumina or anoteher data sequencing)
            info_exon = names[1][2:].split("(")

            exon_nb = info_exon[0]

            amplicon_name = ""
            if len(info_exon) > 1:
                amplicon_name = info_exon[1][:-1]
            exon = Exon(exon_nb,num_chr,start,end,amplicon_name) 
            self.add_exon(gene_name,exon)
            line=IN.readline()

    def add_gene(self, a):
        self.genes.append(a)

    def add_exon(self, b,a):
        for gene in self.genes:
            if gene.name == b:
                gene.add_exon(a)

    def present(self,b):
        for gene in self.genes:
            if gene.name == b:
                return True
        return False

    # Define getter and setter method for attribut self.num_patient
    # Save n°patient and n°barcode
    # self.name represent barcode name
    def get_NumPatient(self):
        if os.path.isfile(globalvar.OUTPUT_DIRECTORY+"/barcodes.json"):
            with open(globalvar.OUTPUT_DIRECTORY+"/barcodes.json") as data_file:
                data = json.load(data_file)
                if self.name in data.keys():
                    num_patient = data[self.name]["sample"]
                    self.num_patient = num_patient
                # In case of bam file from other data sequencing as Illumina or external lab's data
                else:
                    print "No Num patient corresponding to this barcode: "+self.name
                    self.num_patient = "None"
        else:
            print "file barcodes.json not exist"
            self.num_patient = "None"

        return self.num_patient

    def get_TeamName_RunName_ReportName(self):
        # Modif Emeric : saving for header info report file -------------------------------------------------------------------------------------------------------------
        report_name = globalvar.OUTPUT_DIRECTORY.split("/")[0] # Pathway on local computer or another server
        # report_name = globalvar.OUTPUT_DIRECTORY.split("/")[-4] # Pathway for file on server pgm
        TeamName = "None"
        TeamProject = "None"
        RunName = "None"
        
        # save Run name, report name and name of team
        if os.path.isfile(globalvar.INPUT_DIRECTORY+"/expMeta.dat"):
            with open(globalvar.INPUT_DIRECTORY+"/expMeta.dat") as data_file: # Pathway for expMeta.dat file on local computer
                lignes=data_file.readlines()
                RunName = lignes[0].split("=")[1]
                TeamProject = (lignes[3].split("=")[1])
                if re.findall("_\d+", TeamProject)[0]:
                    m = re.findall("_\d+", TeamProject)[0]
                    TeamName = (lignes[3].split("=")[1]).replace(m,"")
                else:
                    print "Incorrect Team name"
                    TeamName = "None"
                    
                if re.findall("SN2-\d+", RunName)[0]:
                    sn2=re.findall("SN2-\d+", RunName)[0]
                    if not re.search(sn2, report_name):
                        print "No report name in file result header , may be you are not on pgm server, please check expMeta.dat file for debug .."
                    else:
                        self.ReportName = report_name.strip("\n")
                        self.RunName = RunName.strip("\n")
                        self.TeamName = TeamName.strip("\n")
        else:
            print "file expMeta.dat not exist"
            self.ReportName = report_name.strip("\n")
            self.RunName = RunName.strip("\n")
            self.TeamName = TeamName.strip("\n")

        return self.ReportName, self.RunName, self.TeamName

        
        # End of saving for header info report file ---------------------------------------------------------------------------------------------------------------------




# ------------------------------------------- DATA LOADING --------------------------------------------

    def data_exons_calculation(self):
        if globalvar.VERBEUX:
            print "Start depth calculation for "+self.name+"...\n"
        IN = open(globalvar.OUTPUT_DIRECTORY+"/tmp/coverage"+self.name, "r")
        line=IN.readline() 
        while line != '':
            infos=line.split()
            num_chr = infos[0]
            position = int(infos[1])
            depth = int(infos[2])
            gene_name,num_ex = self.search_gene_exon_name(position,num_chr)
            self.add_in_depth_list(gene_name,num_ex,depth,position) 
            line = IN.readline()
        for gene in self.genes:
            for exon in gene.exons:
                exon.min_depth()
                exon.max_depth()
                exon.average_depth()                 
                exon.count_covered_bases()
        for gene in self.genes:
            gene.average_coverage_bases()
        if globalvar.VERBEUX:
            print "End depth calculation for "+self.name+"...\n"

    def bias_calculation(self):
        IN_forward = open(globalvar.OUTPUT_DIRECTORY+"/tmp/bias_forward"+self.name, "r")
        line_forward=IN_forward.readline() 
        while line_forward != '':
            infos_forward = line_forward.split()
            num_chr = infos_forward[0]
            position = int(infos_forward[1])
            gene_name,num_ex = self.search_gene_exon_name(position,num_chr)
            depth_forward = int(infos_forward[2])
            self.add_in_depth_list_forward(gene_name,num_ex,depth_forward,position)
            line_forward=IN_forward.readline() 
        IN_forward.close()
        IN_reverse = open(globalvar.OUTPUT_DIRECTORY+"/tmp/bias_reverse"+self.name, "r")
        line_reverse=IN_reverse.readline() 
        while line_reverse != '':
            infos_reverse = line_reverse.split()
            num_chr = infos_reverse[0]
            position = int(infos_reverse[1])
            gene_name,num_ex = self.search_gene_exon_name(position,num_chr)
            depth_reverse = int(infos_reverse[2])
            self.add_in_depth_list_reverse(gene_name,num_ex,depth_reverse,position)
            line_reverse=IN_reverse.readline() 
        IN_reverse.close()
        self.depths_verify()
        for gene in self.genes:
            for exon in gene.exons:
                for i in range(len(exon.list_depths_forward)):
                    bias = 0
                    if exon.list_depths_forward[i] != 0:
                        bias = exon.list_depths_forward[i] / float(exon.list_depths_forward[i] + exon.list_depths_reverse[i])
                    exon.max_bias(bias)

    def add_hotspots_of_mutationsfile(self):
        IN = open(globalvar.HOTSPOT_FILE, "r")
        # test
        #print "hotspot file: "+globalvar.HOTSPOT_FILE
        line=IN.readline() 
        line=IN.readline()
        while line != '':
            infos=line.split()
            num_chr = infos[0]
            start = int(infos[1])
            end = int (infos[2])
            label = infos[3]
            for gene in self.genes: 
                gene.add_hotspot(label,start,end,num_chr)
            line=IN.readline() 

    def data_hotspots_calculation(self,sample):
        IN = open(globalvar.OUTPUT_DIRECTORY+"/tmp/coverage"+sample, "r")
        line=IN.readline() 
        while line != '':
            infos=line.split()
            num_chr = infos[0]
            position = int(infos[1])
            depth = int(infos[2])
            gene_name,hotspot_names = self.search_gene_hotspot_name(position,num_chr)
            if gene_name != None: 
                self.add_in_depth_list_hs(depth,gene_name,hotspot_names,position)
            line=IN.readline()
        for gene in self.genes:
            for hotspot in gene.hotspots:
                hotspot.min_depth()

# ------------------------------------------- HOTSPOT FUNCTION -------------------------------------------- 
        
    def search_gene_hotspot_name(self,position,num_chr):
        for gene in self.genes:
            gene_name,hotspot_names = gene.search_gene_hotspot_name(position,num_chr)
            if gene_name != None:
                return gene.name,hotspot_names
        return None,None

    def add_in_depth_list_hs(self,depth,gene_name,hotspot_names,position):
        for gene in self.genes:
            if gene.name == gene_name:
                gene.add_in_depth_list_hs(hotspot_names,depth,position)

    def select_hotspots_low_coverage(self):
        for gene in self.genes:
            gene.select_hotspots_low_coverage()

# ------------------------------------------- FUNCTION FOR EXON --------------------------------------------

    def search_gene_exon_name(self,position,num_chr):
        for gene in self.genes:
            gene_name,exon_num = gene.search_gene_exon_name(position,num_chr)
            if gene_name != None and exon_num != None:
                return gene.name,exon_num
        return None,None

    def search_name_amplicon(self,num_ex,gene_name):
        for gene in self.genes:
            if gene.name == gene_name:
                for exon in gene.exons:
                    if exon.name == num_ex:
                        return exon.amplicon                          
        return ""
    
    def add_in_depth_list(self,gene_name,num_exon,depth,pos):
        for gene in self.genes:
            if gene.name == gene_name:
                gene.add_in_depth_list(num_exon,depth,pos)

    def add_in_depth_list_forward(self,gene_name,num_exon,depth,pos):
        for gene in self.genes:
            if gene.name == gene_name:
                gene.add_in_depth_list_forward(num_exon,depth,pos)

    def add_in_depth_list_reverse(self,gene_name,num_exon,depth,pos):
        for gene in self.genes:
            if gene.name == gene_name:
                gene.add_in_depth_list_reverse(num_exon,depth,pos)

    def select_exons_low_depth(self):
        for gene in self.genes:
            gene.select_exons_low_coverage()     

    def select_exons_high_bias(self):
        for gene in self.genes:
            gene.select_exons_high_bias()     


    def depths_verify(self):
        for gene in self.genes:
            for exon in gene.exons:
                exon.depths_verify()

# ------------------------------------------- FUNCTION FOR GENE --------------------------------------------   

    def data_genes_calculation(self):
        for gene in self.genes:
            gene.min_depth_gene()
            gene.max_depth_gene()
            gene.average_depth_gene()

# ------------------------------------------------- SAVING --------------------------------------------------

    def save_analysis(self):
        file_date = str(time.strftime('%Y-%m-%d %H:%M',time.localtime()))        
        report_name, RunName, TeamName = self.get_TeamName_RunName_ReportName()
        #OUT = open(globalvar.OUTPUT_DIRECTORY+"/analysis_"+self.name+".txt", "w")
        NumPatient = self.get_NumPatient()
        OUT = open(globalvar.OUTPUT_DIRECTORY+"/"+NumPatient+"_"+self.name+"_analysis.txt", "w")

        # report saved info in the report file header ----------------------------------------------------------------
        OUT.write("Equipe : \t"+TeamName+"\n")
        OUT.write("N°Patient : \t"+NumPatient+"\n")
        OUT.write("Barcode : \t"+self.name+"\n")
        OUT.write("Run : \t"+RunName+"\n")
        OUT.write("Report Name : \t"+report_name+"\n")
        OUT.write("Decora version : \tv1.0.3"+"\n")
        OUT.write("Decora analysis file date : \t"+file_date)
        OUT.write("\n\n")
        # End of reporting saved info in the report file header ------------------------------------------------------


        OUT.write("Label\tDepth min\tDepth moy\tDepth max\tBias")
        for gene in self.genes:
            for exon in gene.exons:
                label = gene.name+"_EX"+exon.name
                OUT.write("\n"+label+"\t"+str(exon.depth_min)+"\t"+str(round(exon.depth_average,2))+"\t"+str(exon.depth_max)+"\t"+str(round(exon.bias_max,2)))
        OUT.close()
        

    def save_hotspot_analysis(self):

        file_date = str(time.strftime('%Y-%m-%d %H:%M',time.localtime()))
        report_name, RunName, TeamName = self.get_TeamName_RunName_ReportName()
        #OUT = open(globalvar.OUTPUT_DIRECTORY+"/hotspot_analysis_"+self.name+".txt", "w")
        NumPatient = self.get_NumPatient()
        OUT = open(globalvar.OUTPUT_DIRECTORY+"/"+NumPatient+"_"+self.name+"_analysis_hotspot.txt", "w")

        # report saved info in the report file header ----------------------------------------------------------------
        OUT.write("Equipe : \t"+TeamName+"\n")
        OUT.write("N°Patient : \t"+NumPatient+"\n")
        OUT.write("Barcode : \t"+self.name+"\n")
        OUT.write("Run : \t"+RunName+"\n")
        OUT.write("Report Name : \t"+report_name+"\n")
        OUT.write("Decora version : \tv1.0.3"+"\n")
        OUT.write("Decora hotspot analysis file date : \t"+file_date)
        OUT.write("\n\n")
        # End of reporting saved info in the report file header ------------------------------------------------------

        OUT.write("Label\tDepth");
        for gene in self.genes:
            for hotspot in gene.hotspots:
                OUT.write("\n"+hotspot.name+"\t"+str(hotspot.depth_min)+"\t")
        OUT.close()
        

    def save_report(self):
        genes = []
        file_date = str(time.strftime('%Y-%m-%d %H:%M',time.localtime()))
        report_name, RunName, TeamName = self.get_TeamName_RunName_ReportName()
        #OUT = open(globalvar.OUTPUT_DIRECTORY+"/report_"+self.name+".txt", "w")
        NumPatient = self.get_NumPatient()
        OUT = open(globalvar.OUTPUT_DIRECTORY+"/"+NumPatient+"_"+self.name+"_report.txt", "w")

        # report saved info in the report file header ----------------------------------------------------------------
        OUT.write("Equipe : \t"+TeamName+"\n")
        OUT.write("N°Patient : \t"+NumPatient+"\n")
        OUT.write("Barcode : \t"+self.name+"\n")
        OUT.write("Run : \t"+RunName+"\n")
        OUT.write("Report Name : \t"+report_name+"\n")
        OUT.write("Decora version : \tv1.0.3"+"\n")
        OUT.write("Decora report file date : \t"+file_date)
        OUT.write("\n\n")
        # End of reporting saved info in the report file header ------------------------------------------------------

        OUT.write("Gene\tPercentage of covered bases")
        if (globalvar.ASOLUTE_MIN):
            OUT.write("\tAbsolute Minimal Depth\tAbsolute Average Depth\tMaximal Depth")
        else:
            OUT.write("\tMinimal Depth\tAverage Depth\tMaximal Depth")
        OUT.write("\tNot-fully-covered Exons")   
        if globalvar.INTRONICS_REGIONS_ANALYSIS:
            OUT.write("\tNot-fully-covered intronic regions")
        if globalvar.BIAS_SELECTION :
                OUT.write("\tExons with high bias") 
        if not globalvar.HOTSPOT_FILE == "":
                OUT.write("\tNot-covered HotSpots")   
        for gene in self.genes:
            OUT.write("\n"+gene.name+"\t"+str(round(gene.coverage,2))+"\t"+str(gene.depth_min)+"\t"+str(round(gene.depth_average,2))+"\t"+str(gene.depth_max)+"\t")
            no_covered_exons = sort_list(gene.no_covered_exons,gene.name,self)
            OUT.write(", ".join(no_covered_exons))
            if globalvar.INTRONICS_REGIONS_ANALYSIS:
                no_covered_introns = sort_list(gene.no_covered_introns,gene.name,self)
                OUT.write("\t")
                OUT.write(", ".join(no_covered_introns))
            if globalvar.BIAS_SELECTION :
                OUT.write("\t")
                wrong_bias_exons = sort_list(gene.wrong_bias_exons,gene.name,self)
                OUT.write(", ".join(wrong_bias_exons))
            if not globalvar.HOTSPOT_FILE == "":
                OUT.write("\t")
                for hs in gene.no_covered_hotspots:
                    min_depth_hs = gene.get_min_depth_hs(hs)
                    OUT.write(hs+" ("+str(min_depth_hs)+"), ")
        OUT.close()
        

    def bed2wig(self):
        NumPatient = self.get_NumPatient()
        OUT = open(globalvar.OUTPUT_DIRECTORY+"/"+NumPatient+"_"+self.name+".wig","w")
        #OUT = open(globalvar.OUTPUT_DIRECTORY+"/"+self.name+".wig","w")
        OUT.write("track type=wiggle_0 name=\""+self.name+"\" visibility=full\n")
        for gene in self.genes:
            for exon in gene.exons:
                OUT.write("\nfixedStep chrom="+exon.chr+" start="+str(exon.start+1)+" step=1 name="+gene.name+"_"+exon.name+"\n")
                for depth in exon.list_depths:
                    OUT.write(str(depth)+"\n")
        OUT.close()


# ------------------------------------------------- SORT --------------------------------------------------

def select_int_alpha_value(value):
    if not re.match('\d$', value[-1:]):
        val_int = int(value[:-1])
        val_alpha = str(value[-1:])
    else : 
        val_int = int(value)
        val_alpha = ""
    return val_int,val_alpha

def sort_list(T,gene,sample):
    tmp_tab = T
    tmp_tab_sorted = []
    while len(tmp_tab) > 0:       
        indice_min = 0
        int_min,alpha_min  = select_int_alpha_value(tmp_tab[indice_min])
        for i in range(1,len(tmp_tab)):
            val_int,val_alpha = select_int_alpha_value(tmp_tab[i])
            if val_int < int_min or (val_int == int_min and val_alpha < alpha_min):
                indice_min = i
                int_min = val_int
                alpha_min = val_alpha 
        amplicon = sample.search_name_amplicon(tmp_tab[indice_min],gene)
        if amplicon != "":
            tmp_tab_sorted.append(str(tmp_tab[indice_min])+"("+amplicon+")")
        else:
            tmp_tab_sorted.append(str(tmp_tab[indice_min]))
        tmp_tab.pop(indice_min)
    return tmp_tab_sorted    
