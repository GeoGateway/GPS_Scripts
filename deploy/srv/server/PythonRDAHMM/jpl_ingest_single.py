#!/usr/local/bin/python
#==========================================================================
# Ingest the JPL dataset into the corresponding databases. 
# Destination data directory and temporary working directory are defined in
# properties.
#
# output: corresponding overall sqlite db file with all data ingested;
#   as well as duplicate-filled sqlite db file for individual stations 
#
# usage:
#   jpl_ingest_single.py 
#
# output:
#   /path/to/rdahmm/jpl_data.sqlite
#   /path/to/rdahmm/jplID.sqlite
#===========================================================================
import os, sys, string 
import sqlite3 as db
import urllib2
from datetime import date
from datetime import timedelta
from datetime import datetime

from properties import properties
from unrsites import unr_sites

from time import strptime
from time import mktime

url_prefix="ftp://sideshow.jpl.nasa.gov/pub/usrs/mbh/point/"
url_suffix=".series"

station_list = wnamsites()

data_path = properties('data_path')
temp_path = properties('temp_path')
model_path = properties('model_path')

datadir = data_path + "JPL_" + dataset.upper() + "/"
dbfile = datadir + "JPL_" + dataset.upper() + ".sqlite"
workdir = temp_path + "JPL_" + dataset.upper() + "/"
#print datadir, dbfile

#Make a data directory for the data set if needed.
if not os.path.exists(datadir):
    cmd = "mkdir -p " + datadir
    os.system(cmd) 
if not os.path.exists(workdir):
    cmd = "mkdir -p " + workdir
    os.system(cmd)

#if the same db file exists, drop it
if os.path.isfile(dbfile):
    #print "deleting old database " + dbfile
    os.remove(dbfile)

# creating/connecting the database 
conn = db.connect(dbfile)
# creating a Cursor
cur = conn.cursor()
# creating tables
# These all go to the main DB file for the network.
sql ="""CREATE TABLE GPSTimeSeries (
      StationID CHAR(4), 
      North Num,
      East Num,
      Up  Num,
      Nsig Num,
      Esig Num, 
      Usig Num,
      Timestamp TEXT,
      UNIQUE (StationID, Timestamp))"""
cur.execute(sql)
sql ="""CREATE TABLE ReferencePositions (
      StationID CHAR(4), 
      Latitude Num,
      Longitude Num, 
      Height Num, 
      UNIQUE (StationID))"""
cur.execute(sql)
conn.commit()

# clear working directory
cmd = "rm -f " + workdir + "*"
os.system(cmd)

#Get the list of stations
stationListUrl="https://sideshow.jpl.nasa.gov/post/tables/table2.html"
wgetcmd="wget -nv -P" + workdir + " "+url
os.system(wgetcmd)

# Read all the lines
with open(workdir + stationID + url_suffix, 'r') as f: 
    station_list = f.readlines()
    
    # Loop over all the stations in the station list.
    for entry in station_list[7:]:
        if 'POS' in entry:
            (stationID, lat, long, height) = entry
            print entry
        else if 'VEL' in entry:
            print "Velocity found:", entry
        else:
            print "weirdness:", entry
        

cur.close()
conn.close()
