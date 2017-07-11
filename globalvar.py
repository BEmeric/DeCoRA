#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright© 2015 Alice CHOURY
# globalvar.py is a part of a free software, DeCoRA : you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANYWARRANTY
# ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program. If not, see http://www.gnu.org/licenses/.

import os

INPUT_DIRECTORY ="."
OUTPUT_DIRECTORY ="results"
DEPTH_THRESHOLD = 30
BIAS_ANALYSIS = False
BIAS_THRESHOLD = 0.2
BIAS_SELECTION = False
HOTSPOT_FILE = ""
CNV_ANALYSIS = False
BEDFILENAME = "../data/ROI.bed"
INTRONIC_REGIONS_ANALYSIS = False
LIMIT_INTRONIC_REGIONS_5P = 0
LIMIT_INTRONIC_REGIONS_3P = 0
CNV_ANALYSIS_ONLY = False
VERBEUX = False
BED2WIG = False
ASOLUTE_MIN = False
exitFlag = 0
NB_TREADS = 6
