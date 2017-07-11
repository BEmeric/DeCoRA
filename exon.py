#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright© 2015 Alice CHOURY
# exon.py is a part of a free software, DeCoRA : you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANYWARRANTY
# ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program. If not, see http://www.gnu.org/licenses/.

import math, globalvar 

class Exon:

# ------------------------------------------- CONSTRUCTOR --------------------------------------------   

    def __init__(self,a,b,c,d,e):
        self.name = a        
        self.amplicon = e
        self.chr = b
        self.start = c
        self.end = d
        self.depth_min_covered = float("inf")
        self.depth_min = float("inf")
        self.depth_max = 0
        self.depth_average = 0
        self.depth_average_covered = 0
        self.bias_max = 0.5
        self.sum_depth = 0
        self.sum_depth_covered = 0
        self.nb_covered_bases = 0
        self.cds_no_covered = False
        self.list_depths = []
        self.list_depths_forward = []
        self.list_depths_reverse = []

    def search_exon_name(self,position,num_chr):
        if self.chr == num_chr:
            if self.start <= position and self.end >= position:
                return self.name
        return None

# ------------------------------------------- FUNCTION FOR BIAS --------------------------------------------   

    def depths_verify(self):
        len_exon = self.end - self.start
        while len_exon > len(self.list_depths_forward):
            self.list_depths_forward.append(0)
        while len_exon > len(self.list_depths_reverse):
            self.list_depths_reverse.append(0)

# ------------------------------------------- FUNCTION FOR DEPTH --------------------------------------------   

    def min_depth(self):
        len_exon = self.end - self.start
        while len_exon > len(self.list_depths):
            self.list_depths.append(0)
        for position in range (len(self.list_depths)):
            depth = self.list_depths[position]
            if self.depth_min > depth :
                self.depth_min = depth
            if self.depth_min < globalvar.DEPTH_THRESHOLD: 
                if globalvar.INTRONICS_REGIONS_ANALYSIS:
                    if position >= globalvar.LIMIT_INTRONIC_REGIONS_5P and position <= (self.end - self.start - 1) - globalvar.LIMIT_INTRONIC_REGIONS_3P:
                        self.cds_no_covered = True
                else:
                    self.cds_no_covered = True
        if self.depth_min_covered > depth and depth >= globalvar.DEPTH_THRESHOLD:
            self.depth_min_covered = depth

    def max_depth(self):
        for position in range (len(self.list_depths)):
            depth = self.list_depths[position]
            if self.depth_max < depth :
                self.depth_max = depth

    def average_depth(self):
        for position in range (len(self.list_depths)):
            depth = self.list_depths[position]
            if depth >= globalvar.DEPTH_THRESHOLD:
                self.sum_depth_covered = self.sum_depth_covered + depth
            self.sum_depth = self.sum_depth + depth
        self.depth_average_covered = float(self.sum_depth_covered) / float(self.end - self.start + 1)
        self.depth_average = float(self.sum_depth) / float(self.end - self.start + 1)
    
# ------------------------------------------- FUNCTION FOR BIAS --------------------------------------------   

    def max_bias(self,bias):
        if math.fabs(0.5-self.bias_max) < math.fabs(0.5-bias):
            self.bias_max = bias

# ------------------------------------------- FUNCTION FOR DEPTHS LIST --------------------------------------------  
 
    def add_in_depth_list(self,depth,pos):
        while pos > self.start + 1 + len(self.list_depths):
            self.list_depths.append(0)
        self.list_depths.append(depth)
    
    def add_in_depth_list_forward(self,depth,pos):
        while pos > self.start + 1 + len(self.list_depths_forward):
            self.list_depths_forward.append(0)
        self.list_depths_forward.append(depth)
    
    def add_in_depth_list_reverse(self,depth,pos):
        while pos > self.start + 1 + len(self.list_depths_reverse):
            self.list_depths_reverse.append(0)
        self.list_depths_reverse.append(depth)

# ------------------------------------------- FUNCTION FOR COVERAGE -------------------------------------------- 

    def count_covered_bases(self):
        for position in range (len(self.list_depths)):
            depth = self.list_depths[position]
            if depth >= globalvar.DEPTH_THRESHOLD:
                self.nb_covered_bases = self.nb_covered_bases + 1
