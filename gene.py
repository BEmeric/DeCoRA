#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright© 2015 Alice CHOURY
# gene.py is a part of a free software, DeCoRA : you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANYWARRANTY
# ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program. If not, see http://www.gnu.org/licenses/.

import math, globalvar, sys
from exon import *
from hotspot import *

class Gene:
# ----------------------------------------------- CONSTRUCTEUR -----------------------------------------------
    def __init__(self,a):
        self.exons = []
        self.no_covered_exons = []
        self.no_covered_introns = []
        self.name = a
        self.depth_min_covered = float("inf")
        self.depth_min = float("inf")
        self.depth_max = 0
        self.depth_average_covered = 0
        self.depth_average = 0
        self.coverage = 0
        self.hotspots = []
        self.no_covered_hotspots = []
        self.wrong_bias_exons = []

# ----------------------------------------------- HOTSPOT FUNCTION -----------------------------------------------

    def add_hotspot(self,hotspot_label,start,end,num_chr):
        for exon in self.exons:
            if exon.chr == num_chr:
                if exon.start <= start and exon.end >= end:
                    hotspot = HotSpot(hotspot_label,num_chr,start,end)
                    self.hotspots.append(hotspot)

    def add_in_depth_list_hs(self,hotspot_names,depth,pos):
        for hotspot in self.hotspots:
            for hotspot_name in hotspot_names : 
                if hotspot.name == hotspot_name:
                    hotspot.add_in_depth_list(depth,pos)

    def search_gene_hotspot_name(self,position,num_chr):
        hotspot_names = []
        gene_name = None
        for hotspot in self.hotspots:
            hotspot_name = hotspot.search_hotspot_name(position,num_chr)
            if hotspot_name != None:
                hotspot_names.append(hotspot_name)
                gene_name = self.name
        return gene_name,hotspot_names

    def select_hotspots_low_coverage(self):
        for hotspot in self.hotspots:
            if hotspot.depth_min < globalvar.DEPTH_THRESHOLD:
                self.no_covered_hotspots.append(hotspot.name)

    def get_min_depth_hs(self,hs):
        for hotspot in self.hotspots:
            if hotspot.name == hs:
                return hotspot.depth_min
        return 0

# ----------------------------------------------- EXON FUNCTION -----------------------------------------------

    def add_exon(self, a):
        self.exons.append(a)

    def search_gene_exon_name(self,position,num_chr):
        for exon in self.exons:
            exon_num = exon.search_exon_name(position,num_chr)
            if exon_num != None:
                return self.name,exon_num
        return None,None

    def add_in_depth_list(self,num_exon,depth,pos):
        for exon in self.exons:
            if exon.name == num_exon:
                exon.add_in_depth_list(depth,pos)

    def add_in_depth_list_forward(self,num_exon,depth,pos):
        for exon in self.exons:
            if exon.name == num_exon:
                exon.add_in_depth_list_forward(depth,pos)

    def add_in_depth_list_reverse(self,num_exon,depth,pos):
        for exon in self.exons:
            if exon.name == num_exon:
                exon.add_in_depth_list_reverse(depth,pos)

    def select_exons_high_bias(self):
        for exon in self.exons:
            if exon.depth_min > globalvar.DEPTH_THRESHOLD and (exon.bias_max > 0.5 + globalvar.BIAS_THRESHOLD or exon.bias_max < 0.5 - globalvar.BIAS_THRESHOLD):
                self.wrong_bias_exons.append(exon.name)

    def select_exons_low_coverage(self):
        for exon in self.exons:
            if exon.depth_min < globalvar.DEPTH_THRESHOLD and exon.cds_no_covered:
                self.no_covered_exons.append(exon.name)
            if exon.depth_min < globalvar.DEPTH_THRESHOLD and not exon.cds_no_covered:
                self.no_covered_introns.append(exon.name)

# ----------------------------------------------- GENE FUNCTION -----------------------------------------------

    def min_depth_gene(self):
        for exon in self.exons:
            if globalvar.ASOLUTE_MIN:
                if self.depth_min > exon.depth_min:
                    self.depth_min = exon.depth_min
            else:
                if self.depth_min > exon.depth_min_covered:
                    self.depth_min = exon.depth_min_covered
        if self.depth_min == float("inf") or self.depth_min == "inf":
            self.depth_min=0

    def max_depth_gene(self):
        for exon in self.exons:
            if self.depth_max < exon.depth_max:
                self.depth_max = exon.depth_max

    def average_depth_gene(self):
        total_depth = 0      
        length_gene = 0  
        for exon in self.exons:
            if globalvar.ASOLUTE_MIN:
                total_depth = total_depth + exon.sum_depth
                length_gene = length_gene + exon.end - exon.start + 1
            else :
                total_depth = total_depth + exon.sum_depth_covered
                length_gene = length_gene + exon.nb_covered_bases
        
        # Modif EB: test debug
        if length_gene==0:
            # print self.name
            # print "A corriger dans le code : " +str(total_depth)+" / "+str(length_gene)
            length_gene=1
            # print "Correction apportée dans le code: "+str(total_depth)+" / "+str(length_gene)
            # sys.exit(1)
        self.depth_average = float(total_depth) / float(length_gene)

    def average_coverage_bases(self):
        sum_coverage = 0
        sum_len = 0
        for exon in self.exons:
            sum_coverage = sum_coverage + exon.nb_covered_bases
            sum_len = sum_len + (exon.end - exon.start)
        self.coverage = float(sum_coverage) / float(sum_len) * 100.0        
