#!/usr/bin/env python3
# Author:  Michael Heflin
# Date:  September 20, 2016
# Organization:  JPL, Caltech

from __future__ import print_function

prolog="""
**PROGRAM**
    getCoseismic.py
      
**PURPOSE**
    Make kml file from coseismic break estimates

**USAGE**
"""
epilog="""
**EXAMPLE**
    getCoseismic.py --lat 33 --lon -115 --width 2 --height 2 -o coseismic.kml -t 2010-04-08
               
**COPYRIGHT**
    | Copyright 2016, by the California Institute of Technology
    | United States Government Sponsorship acknowledged
    | All rights reserved

**AUTHORS**
    | Michael Heflin
    | Jet Propulsion Laboratory
    | California Institute of Technology
    | Pasadena, CA, USA
"""

# Import modules
import os
import sys
import math
import time
import datetime
import calendar
import argparse
import subprocess
import urllib.request

def runCmd(cmd):
    '''run a command'''

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,executable='/bin/bash')
    (out, err) = p.communicate()
    if p.returncode != 0:
        raise UserWarning('failed to run {}\n{}\n'.format(cmd.split()[0],
            err))
    return out

def _getParser():
    parser = argparse.ArgumentParser(description=prolog,epilog=epilog,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', action='store', dest='output',required=True,help='output kml file')
    parser.add_argument('-t', action='store', dest='epoch',required=True,help='epoch')
    parser.add_argument('--lat', action='store', dest='lat',required=True,help='center latitude in degress')
    parser.add_argument('--lon', action='store', dest='lon',required=True,help='center longitude in degrees')
    parser.add_argument('--width', action='store', dest='width',required=True,help='width in degrees')
    parser.add_argument('--height', action='store', dest='height',required=True,help='height in degrees')
    parser.add_argument('--dt', action='store', dest='dt',required=False,help='optional time window in years')
    parser.add_argument('--scale', action='store', dest='scale',required=False,help='optional scale for offsets')
    return parser

def main():

    # Read command line arguments
    parser = _getParser()
    results = parser.parse_args()

    # Set bounds
    latmin = float(results.lat) - float(results.height)/2
    latmax = float(results.lat) + float(results.height)/2
    lonmin = float(results.lon) - float(results.width)/2
    lonmax = float(results.lon) + float(results.width)/2

    # Set epoch
    if (len(results.epoch) == 10):
        results.epoch = datetime.datetime.strptime(results.epoch,"%Y-%m-%d").strftime("%y%b%d").upper()
        ntime = time.strptime(results.epoch,"%y%b%d")
        jtime = time.strptime("2000JAN01","%Y%b%d")
        ytime = float(calendar.timegm(ntime)-calendar.timegm(jtime))
        ytime = ytime/(86400.*365.25)
        ytime = ytime + 2000.

    # Set time window
    dt = 0.02
    if (results.dt != None):
        dt = float(results.dt)

    # Set scale
    scale = 320 
    if (results.scale != None):
        scale = float(results.scale)

    # Read table of positions
    response1 = urllib.request.urlopen('http://sideshow.jpl.nasa.gov/post/tables/table2.html')
    lines = response1.readall().decode('utf-8').splitlines()

    # Read table of breaks
    response2 = urllib.request.urlopen('http://sideshow.jpl.nasa.gov/post/tables/table3.html')
    breaks = response2.readall().decode('utf-8').splitlines()

    # Start kml file
    outFile = open(results.output,'w')
    print("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",file=outFile)
    print("<kml xmlns=\"http://www.opengis.net/kml/2.2\">",file=outFile)
    print(" <Folder>",file=outFile)

    # Add markers and vectors
    for i in range(0,len(lines)):
        test = lines[i].split()
        if (len(test) == 8):
            if (test[1] == 'POS'):
                lon = float(test[3])
                lat = float(test[2])
                if ((lon > lonmin) & (lon < lonmax) & (lat > latmin) & (lat < latmax)):
                    lonv = 0
                    latv = 0
                    for j in range(0,len(breaks)):
                        test2 = breaks[j].split()
                        if (len(test2) == 8):
                            if ((test2[0] == test[0]) & (float(test2[1]) > ytime-dt ) & (float(test2[1]) < ytime+dt)):
                                lonv = lonv + float(test2[3])
                                latv = latv + float(test2[2])

                    # Draw marker 
                    print("  <Placemark>",file=outFile)
                    print("   <description><![CDATA[",file=outFile)
                    print("    <a href=\"http://sideshow.jpl.nasa.gov/post/links/{:s}.html\">".format(test[0]),file=outFile)
                    print("     <img src=\"http://sideshow.jpl.nasa.gov/post/plots/{:s}.jpg\" width=\"300\" height=\"300\">".format(test[0]),file=outFile)
                    print("    </a>",file=outFile)
                    print("   ]]></description>",file=outFile)
                    print("   <Style><IconStyle>",file=outFile)
    #               print("    <color>FF00FF14</color>",file=outFile)
                    print("    <color>FF78FF78</color>",file=outFile)
                    print("    <scale>0.25</scale>",file=outFile)
    #               print("    <Icon><href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href></Icon>",file=outFile)
                    print("    <Icon><href>http://maps.google.com/mapfiles/kml/shapes/road_shield3.png</href></Icon>",file=outFile)
                    print("   </IconStyle></Style>",file=outFile)
                    print("   <Point>",file=outFile)
                    print("    <coordinates>",file=outFile)
                    print("     {:f},{:f},0".format(lon,lat),file=outFile)
                    print("    </coordinates>",file=outFile)
                    print("   </Point>",file=outFile)
                    print("  </Placemark>",file=outFile)

                    # Draw vector    
                    print("  <Placemark>",file=outFile)
                    print("   <Style><LineStyle>",file=outFile)
                    print("    <color>FF1400FF</color>",file=outFile)
                    print("    <width>2</width>",file=outFile)
                    print("   </LineStyle></Style>",file=outFile)
                    print("   <LineString>",file=outFile)
                    print("   <coordinates>",file=outFile)
                    print("   {:f},{:f},0".format(lon,lat),file=outFile)
                    print("   {:f},{:f},0".format(lon+(lonv)/scale,lat+(latv)/scale),file=outFile)
                    print("    </coordinates>",file=outFile)
                    print("   </LineString>",file=outFile)
                    print("  </Placemark>",file=outFile)

    # Finish kml file
    print(" </Folder>",file=outFile)
    print("</kml>",file=outFile)
    outFile.close()

if __name__ == '__main__':
    main()