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

# These are global variables
url_prefix="ftp://sideshow.jpl.nasa.gov/pub/usrs/mbh/point/"
url_suffix=".series"

data_path = properties('data_path')
temp_path = properties('temp_path')
model_path = properties('model_path')

datadir = data_path + "JPL" + "/"
dbfile = datadir + "JPL" + ".sqlite"
workdir = temp_path + "JPL" + "/"
#print datadir, dbfile

# The list of stations
stationList=[]
stationListUrl="https://sideshow.jpl.nasa.gov/post/tables/table2.html"

def prepareDirectories():
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

    # clear working directory
    cmd = "rm -f " + workdir + "*"
    os.system(cmd)

def createNetworkDBTables():
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

def readStationList():
    #Get the list of stations
    wgetcmd="wget -nv -P" + workdir + " " + stationListUrl
    os.system(wgetcmd)

    # Read all the lines
    with open(workdir+"/"+"table2.html", 'r') as f: 
        station_list = f.readlines()
        
        # Loop over all the stations in the station list.
        # The following is some ridiculous code to convert the string in the file to the data we need
        for entry in station_list[9:]:
            if "POS" in entry:
                cleanedEntry=""
                for entryElement in entry.split():
                    if entryElement != " ": 
                        cleanedEntry+=entryElement+ ","
                        # print cleanedEntry.split(",")        
                (stationID, pos, lat, long, height, sn, se, sv, end) = cleanedEntry.split(",")
                stationList.append(stationID)
                # Insert these now into the DB
                sql = "INSERT INTO ReferencePositions (StationID, Latitude, Longitude, Height) "
                sql += " VALUES ('%s', '%s', '%s', '%s')" % (stationID, lat, long, height)
                cur.execute(sql)
                conn.commit()
    f.close()

# Create a DB for each station    
def createDbForStations(stationID):
    station_dbfile = datadir + stationID + ".sqlite"

    #If there is already a db file, delete it.
    if os.path.isfile(station_dbfile):
        #print "deleting old station database " + station_dbfile
        os.remove(station_dbfile)

    # Create the table for the station's GPS time series.
    station_conn = db.connect(station_dbfile)
    station_cur = station_conn.cursor()
    station_sql ="""CREATE TABLE StationGPSTimeSeries (
           North Num,
           East Num,
           Up  Num,
           Nsig Num,
           Esig Num, 
           Usig Num,
           Timestamp TEXT,
           Interploated INT Default 0,
           UNIQUE(Timestamp))"""
    station_cur.execute(station_sql)
    station_conn.commit()

    # Download the station file
    remoteFile=url_prefix+stationID+url_suffix
    wgetcmd="wget -nv -P" + workdir + " " + remoteFile
    os.system(wgetcmd)
        
    #Insert each line of the file into the station's DB.
    with open(workdir + stationID + url_suffix, 'r') as f: 
        data = f.readlines()
	last_line = ""
        # Create two DB tables per station.
        # GPSTimeSeries: this is just the data as is.
        # StationGPSTimeSeries: this is the data with missing values inserted.
	for line in data:
	    record = string.split(line)
            # print record
            (decimalDate,east,north,vert,esig,nsig,vsig,encorr,evcorr,nvcorr,time2,year,month,day,hour, minute,sec)=record
            timeStamp=datetime(int(year),int(month),int(day),int(hour),int(minute),int(sec));
            sql = "INSERT INTO GPSTimeSeries (StationID, North, East, Up, Nsig, Esig, Usig, Timestamp) "
            sql += " VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')" % (stationID, north, east, vert, nsig, esig, vsig, timeStamp)        
            
	    try: 
                cur.execute(sql)
	    except db.IntegrityError as error: 
	    	print "Source data error: ", sql
		print "sqlite3.IntegrityError: ", error, "\n"
		continue

            # Create a separate table that fills in missing data.
            if last_line == "":
                last_line = line
            else:
	        last_record = string.split(last_line)
                (ldecimalDate,least,lnorth,lvert,lesig,lnsig,lvsig,lencorr,levcorr,lnvcorr,ltime2,lyear,lmonth,lday,lhour,lminute,lsec)=last_record
                lastTimeStamp=datetime(int(lyear),int(lmonth),int(lday),int(lhour),int(lminute),int(lsec));
                # if missing days from last to current, fill with last
                for i in range(1, (timeStamp - lastTimeStamp).days):
                    ts = lastTimeStamp + timedelta(days=i)
                    interploated = 1
                    station_sql = "INSERT INTO StationGPSTimeSeries (North, East, Up, Nsig, Esig, Usig, Timestamp, Interploated) "
                    station_sql += " VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')" % (lnorth, least, lvert, lnsig, lesig, lvsig, ts, interploated)
                    station_cur.execute(station_sql)
		last_line = line

    station_conn.commit()
    station_cur.close()
    station_conn.close()

#
# This is the code that is executed
#

prepareDirectories()

# creating/connecting the database 
conn = db.connect(dbfile)
# creating a Cursor
cur = conn.cursor()
                
createNetworkDBTables()
readStationList()

for station in stationList:
    createDbForStations(station)

# create index 
sql = "CREATE INDEX idx_StationID ON GPSTimeSeries(StationID)"
cur.execute(sql)
sql = "CREATE INDEX idx_Timestamp ON GPSTimeSeries(Timestamp)"
cur.execute(sql)
sql = "CREATE INDEX idx_RefStationID ON ReferencePositions(StationID)"
cur.execute(sql)
    
cur.close()
conn.close()
