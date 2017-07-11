#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright© 2015 Alice CHOURY
# hotspot.py is a part of a free software, DeCoRA : you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANYWARRANTY
# ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program. If not, see http://www.gnu.org/licenses/.

class HotSpot:
    def __init__(self,a,b,c,d):
        self.name = a
        self.chr = b
        self.start = c
        self.end = d
        self.depth_min = float("inf")
        self.list_depths = []

    def search_hotspot_name(self,position,num_chr):
        if self.chr == num_chr:
            if self.start <= position and self.end >= position:
                return self.name
        return None

    def min_depth(self):
        len_hs = self.end - self.start
        while len_hs > len(self.list_depths):
		    self.list_depths.append(0)
        for depth in self.list_depths :
            if self.depth_min > depth :
                self.depth_min = depth

    # add depth in list if the position is between the start and the end of the hotspot
    # add "0" if the edge of hotspot is not covered, then add the depth 
    def add_in_depth_list(self,depth,pos):
        while pos > self.start + 1 + len(self.list_depths):
            self.list_depths.append(0)
        self.list_depths.append(depth)
