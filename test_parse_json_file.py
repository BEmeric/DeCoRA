#!/usr/bin/python

import Queue, threading, sys, globalvar, os, argparse, re
from datetime import *

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
	header = "Gene"+"\t"+"Exon"+"\t"+"Echantillon"+"\t"+"Nombre d'echantillons"+"\n"
	if globalvar.VERBEUX:
		print "writting results of no covered exons"
	RES = open('no_covered_exons_report.txt', 'w')
	RES.write(header)
	RES.write(res)
	RES.close()
	if globalvar.VERBEUX:
		print "End of Computing for no covoregae exons"

Compute_no_covered_exon()


















