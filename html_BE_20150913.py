#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright© 2015 Alice CHOURY
# html.py is a part of a free software, DeCoRA : you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANYWARRANTY
# ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program. If not, see http://www.gnu.org/licenses/.

import globalvar, json, os

def write_html():
    dict_sample = read_json()

    OUT = open(globalvar.OUTPUT_DIRECTORY+"/annotatic_block.html", "w")
    
    OUT.write("<?xml version=\"1.0\" encoding=\"iso-8859-1\"?> \n")
    OUT.write("<!DOCTYPE html> \n")
    OUT.write("<html> \n")
    OUT.write("<head> \n")    
    OUT.write("<link rel=\"stylesheet\" media=\"all\" href=\"/site_media/resources/bootstrap/css/bootstrap.min.css\" />\n")
    OUT.write("<link href=\"/site_media/resources/kendo/styles/kendo.common.min.css\" rel=\"stylesheet\" />\n")
    OUT.write("<link href=\"/site_media/resources/less/kendo.tb.min.css\" rel=\"stylesheet\" />\n")
    #OUT.write("<link type=\"text/css\" rel=\"stylesheet\" href=\"/site_media/resources/styles/tb-layout.css\" />\n")
    OUT.write("<link type=\"text/css\" rel=\"stylesheet\" href=\"/site_media/resources/styles/tb-styles.min.css\" />\n")
    OUT.write("<link rel=\"stylesheet\" type=\"text/css\" href=\"/site_media/stylesheet.css\"/>\n")
    OUT.write("<link rel=\"stylesheet\" type=\"text/css\" href=\"/site_media/resources/styles/print.css\" media=\"print\" />\n")
    OUT.write("<link rel=\"stylesheet\" type=\"text/css\" href=\"/site_media/resources/styles/report.css\" media=\"screen\" />\n")
    OUT.write("</head> \n")
    OUT.write("<title>DeCoRA Summary</title> \n")
    OUT.write("<body> \n")
    OUT.write("<div class=\"topbar\"><div class=\"logoholder\"> \n")
    OUT.write("\t<a href=\"http\://www.iontorrent.com/\"><img src=\"/site_media/images/raw_name_small.png\" \n")
    OUT.write("\t\t alt=\"IonTorrent Systems, Inc.\" style=\"border\:none\;\"/></a> \n")
    OUT.write("</div>\n")
    OUT.write("<div style=\"width\:1040px\;margin-left\:auto\;margin-right\:auto\;height\:100%\"> \n")
    OUT.write("\t<h1><center>DeCoRA Summary</center></h1> \n")
    OUT.write("\t<div> \n")
    OUT.write("\t <br/> \n")
    OUT.write("\t <style type=\"text/css\"> \n")
    OUT.write("\t\tth {text-align\:center\;width=100%} \n")
    OUT.write("\t\ttd {text-align\:right\;width=100%} \n")
    OUT.write("\t\t.help {cursor\:help\; border-bottom\: 1px dotted \#A9A9A9} \n")
    OUT.write("\t </style> \n")
    OUT.write("<div class=\"k-widget k-grid\"> \n")
    OUT.write("\t<table class=\"table-striped\" style=\"width:100%\"> \n")
    OUT.write("\t\t<thead class=\"k-grid-header\"> \n")
    OUT.write("\t\t<tr style=\"width:100%\"> \n")
    OUT.write("\t\t <th class=\"k-header\"><span class=\"help\" title=\"The barcode ID for each set of reads.\">Barcode Name</span></th> \n")
    OUT.write("\t<th class=\"k-header\"><span class=\"help\" title=\"Sample name\">Sample Name</span></th> \n")    
    OUT.write("\t<th class=\"k-header\"><span class=\"help\" title=\"The minimum depth and the maximum bias of every ROI for the sample\">Detailed analysis</span></th>\n")
    OUT.write("\t<th class=\"k-header\"><span class=\"help\" title=\"Analysis of file of mutations\">Analysis of mutations</span></th>\n") 
    OUT.write("\t<th class=\"k-header\"><span class=\"help\" title=\"File for IGV\">Wiggle file</span></th> \n")
    OUT.write("\t<th class=\"k-header\"><span class=\"help\" title=\"Report of sequencing quality\">Report file</span></th> \n")
    OUT.write("\t\t</tr> \n")
    OUT.write("\t\t</thead> \n")
    OUT.write("\n")
    OUT.write("\n")

    write_html_line(OUT,dict_sample)
    
    OUT.write("\t </table> \n")
    OUT.write(" </div> \n")
    OUT.write("</div> \n")
    OUT.write("\n")
    OUT.write("<a class=\"btn\" href=\"/plugins/DeCoRA/data/aide_IGV.pdf\"> Help for WIGGLE files in IGV </a> \n")
    if globalvar.CNV_ANALYSIS:
        OUT.write("<a class=\"btn\" href=\"CNV_predictions.txt\"> CNV_predictions.txt </a> \n")
    OUT.write("<div class=\"footer\"> \n")
    OUT.write("\t<a href=/plugins/DeCoRA/data/terms-of-use.txt>Terms of Use</a> \n")
    OUT.write("\t<br/>Copyright &copy; 2015<a> Alice CHOURY</a> \n")
    OUT.write("</div> \n")
    OUT.write("</font></body></html>\n")
    OUT.close()

def write_html_line(OUT,dict_sample):
    for barcode in dict_sample:
        sample = dict_sample[barcode]
        OUT.write("\t\t<tr style=\"width:100%\"> \n")
        OUT.write("\t\t<td>"+barcode+"</td>\n")
        OUT.write("\t\t<td style=\"text-align\:left\">"+sample+"</td> \n")
        OUT.write("\t\t<td><a class=\"btn\" href=\"analysis_"+barcode+".txt\">analysis_"+barcode+".txt</a></td> \n")
        if (not globalvar.HOTSPOT_FILE == ""):
            OUT.write("\t\t<td><a class=\"btn\" href=\"hotspot_analysis_"+barcode+".txt\">hotspot_analysis_"+barcode+".txt</a></td> \n")
        else:
            OUT.write("\t\t<td></td> \n")
        if (globalvar.BED2WIG):
            OUT.write("\t\t<td><a class=\"btn\" href=\""+barcode+".wig\">"+barcode+".wig</a></td> \n")
        else:
            OUT.write("\t\t<td></td> \n")
        OUT.write("\t\t<td><a class=\"btn\" href=\"report_"+barcode+".txt\">report_"+barcode+".txt</a></td> \n")
        OUT.write("\t\t</tr> \n")
        OUT.write("\n")

def read_json():
    dict_sample = {}
    with open(globalvar.OUTPUT_DIRECTORY+"/barcodes.json") as data_file:    
        data = json.load(data_file)
        for barcode in data:
            if barcode != "nomatch":
                dict_sample[barcode]=data[barcode]["sample"]
    return dict_sample