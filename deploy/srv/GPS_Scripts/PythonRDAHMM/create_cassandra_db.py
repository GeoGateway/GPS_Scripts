# #!/usr/local/bin/python
#==========================================================================
# Script to load the final json file into Cassandra DB.

# usage:
	 # python create_cassandra_db scripps_dataset_name
	 # eg: python create_cassandra_db.py UNR_SPLICE

# output:
#   Caasandra DB keyspace named gps_timeseries_database which contains
#   three tables-
#   collections_meta_network: Data for that station
#   collections_meta_stations: Meta data for each station
#   collections_time_series_stations: Timeseries for each station
#   The database is indexed based on latitude and longitude
#===========================================================================
import os, sys, string, re, json
from datetime import date, datetime, timedelta, time


from properties import properties
from cassandra.cluster import Cluster
from cassandra.policies import TokenAwarePolicy
from cassandra.policies import DCAwareRoundRobinPolicy
from cassandra.policies import RetryPolicy

# Some useful global constants
today = datetime.today()
serverName = "gf9.ucs.indiana.edu"
updateTime = str(today.strftime("%Y-%m-%dT%H:%M:%S"))
beginDate = "1994-01-01"
endDate = str(today.strftime("%Y-%m-%d"))
# endDate = '2016-10-19'
centerLng = "-119.7713889"
centerLat = "36.7477778"
stateChangeNumTxtFile = "stateChangeNums.txt"
stateChangeNumJsInput = "stateChangeNums.txt.jsi"
allStationInputName = "all_stations.all.input"
filters = "Fill_Missing"

# Used to separate parts of the station name
SEPARATOR_CHARACTER="_"
NO_DATA_TIME_STAMP="22:22:22"
FINAL_PATH=properties('eval_path')


def setStationId(stationList, stationData):
    #Get the station name.
    stationName=stationList.split(SEPARATOR_CHARACTER)[2];

    stationData['id'] = stationName
    stationData['pro_dir'] = "daily_project_" + stationName + "_" + endDate
    stationData['AFile'] = "daily_project_" + stationName + ".A"
    stationData['BFile'] = "daily_project_" + stationName + ".B"
    stationData['InputFile'] = "daily_project_" + stationName + "_" + endDate + ".all.input"
    stationData['RawInputFile'] = "daily_project_" + stationName + "_" + endDate + ".all.raw"
    stationData['SwfInputFile'] = "daily_project_" + stationName + "_" + endDate + ".plotswf.input"
    stationData['DygraphsInputFile'] = "daily_project_" + stationName + "_" + endDate + ".dygraphs.js"
    stationData['LFile'] = "daily_project_" + stationName + ".L"
    stationData['XPngFile'] = "daily_project_" + stationName + "_" + endDate + ".all.input.X.png"
    stationData['YPngFile'] = "daily_project_" + stationName + "_" + endDate + ".all.input.Y.png"
    stationData['ZPngFile'] = "daily_project_" + stationName + "_" + endDate + ".all.input.Z.png"
    stationData['XTinyPngFile'] = "daily_project_" + stationName + "_" + endDate + ".all.input.X_tiny.png"
    stationData['YTinyPngFile'] = "daily_project_" + stationName + "_" + endDate + ".all.input.Y_tiny.png"
    stationData['ZTinyPngFile'] = "daily_project_" + stationName + "_" + endDate + ".all.input.Z_tiny.png"
    stationData['PiFile'] = "daily_project_" + stationName + ".pi"
    stationData['QFile'] = "daily_project_" + stationName + "_" + endDate + ".all.Q"
    stationData['MaxValFile'] = "daily_project_" + stationName + ".maxval"
    stationData['MinValFile'] = "daily_project_" + stationName + ".minval"
    stationData['RangeFile'] = "daily_project_" + stationName + ".range"
    stationData['ModelFiles'] = "daily_project_" + stationName + ".zip"
    stationData['RefFile'] = "daily_project_" + stationName + ".input.ref"
    return


def setStationStartDate(stationDir, stationData):
    startFileName = stationDir + "daily_project_" + stationData['id'] + ".input.starttime"
    if (os.path.isfile(startFileName)):
        with open(startFileName,"r") as startFile:
            startDate = startFile.readline().rstrip()
        startFile.close()
    else:
		startDate = "1994-01-01"
    stationData['start_date'] = startDate
    return

def setStationRefLatLonHgt(stationDir, stationData):
    refFileName = stationDir + stationData['RefFile']
    refLat=""
    refLon=""
    refHgt=""
    if (os.path.isfile(refFileName)):
        with open(refFileName,"r") as refFile:
            refParts=refFile.readline().split(" ")
            refLat=refParts[0]
            refLon=refParts[1]
            refHgt=refParts[2].rstrip() # Have to chomp off the final \n
        refFile.close()
    else:
        refLat="1.0"
        refLon="2.0"
        refHgt="-1.0"

    stationData['lat'] = refLat      
    stationData['long'] = refLon      
    stationData['height'] = refHgt
    return


def setStatusChanges(stationDir, stationData):
    # Open the .all.Q and the .all.raw files.  We get the state from the first and
    # the data from the second. 
    # TODO: for now, we assume these files always exist
    qFileName = stationDir + stationData['QFile']
    rawFileName = stationDir + stationData['RawInputFile']
   



    # Bail out if the required files don't exist   
    if((not os.path.isfile(qFileName)) or (not os.path.isfile(rawFileName))): 
        return 
    qFile = open(qFileName,"r")
    rawFile = open(rawFileName,"r")

    stateChanges = []
    changeCount = 0
    # Now step through the Q file looking for state changes
    # If we find a state change, get the date from the raw file
    # We will save these to the string stateChangeArray since we
    # need to record in latest-first order
    qline1 = qFile.readline()
    rline1 = rawFile.readline()        
    while True:
        eventData = {}
        qline2 = qFile.readline()
        rline2 = rawFile.readline()
        if not qline2: break
        
        # See if qline1 and qline2 are the same.  If so, extract the dates from rline1 and rline2
        # The line splits below are specific to the raw file line format.
        if (qline1.rstrip() != qline2.rstrip()):
            eventdate = rline2.split(" ")[1] 
            eventdate = eventdate.split("T")[0]
            oldstate = qline1.rstrip()
            newstate = qline2.rstrip()
            eventData['date'] = eventdate
            eventData['from'] = oldstate
            eventData['to'] = newstate 
            stateChanges.append(eventData)
            changeCount += 1

        # Make the previous "next" lines the "first" lines for the next comparison
        qline1=qline2
        rline1=rline2

    stationData['status_changes'] = stateChanges
    stationData['change_count'] = changeCount

    # Clean up
    qFile.close
    rawFile.close
    return


def setTimesNoData(stationDir, stationData):
    rawFileName = stationDir + stationData['RawInputFile']

    # Required file doesn't exist so bail out
    if(not os.path.isfile(rawFileName)): return
    rawFile = open(rawFileName, "r")
    
    noDataRanges = []
    noDataCount = 0
    noDataEvent = {}

    # We need to set a no-data range from beginDate (for the epoch, 1994-01-01) to the day before
    # our first data point for this station.  If the station has data before 1994-01-01, then 
    # ignore.
    firstDataDateParts=rawFile.readline().split(" ")[1].split("T")[0].split("-");

    beginEpoch=date(1994,1,1)

    #Convert this into a data object
    dayMinusOne=date(int(firstDataDateParts[0]),int(firstDataDateParts[1]),int(firstDataDateParts[2]))
    dayMinusOne-=timedelta(days=1)
    if(dayMinusOne > beginEpoch): 
        dayMinusOneString=dayMinusOne.isoformat()
        noDataEvent['to'] = dayMinusOneString
        noDataEvent['from'] = beginDate
        noDataCount += 1

    #Reset the "raw" file to the beginning
    rawFile.seek(0)

    # Step through the file to find the starting and ending dates with no data.
    # By convention, this occurs when the line has a timestamp T22:22:22.  Also, by
    # convention, we will record the latest to earliest dates with no data.

    while True:
        noDataEvent = {}
        nodata=False
        rline1=rawFile.readline()
        if not rline1: break

        # Get the date and timestamp, following format conventions
        fulleventdate1=rline1.split(" ")[1]
        eventdate1=fulleventdate1.split("T")[0]
        timestamp1=fulleventdate1.split("T")[1]

        # See if we have detected a no-data line
        if(timestamp1==NO_DATA_TIME_STAMP):
            nodata=True
            #Keep eventdate1 in case this is an isolated no-data line.
            eventdate_keep=eventdate1

            # We have a no-data line, so step ahead until the 
            # no-data line ends.
            while(nodata):
                rline2=rawFile.readline()
                if not rline2: break
                fulleventdate2=rline2.split(" ")[1]
                eventdate2=fulleventdate2.split("T")[0]
                timestamp2=fulleventdate2.split("T")[1]
                if(timestamp2!=NO_DATA_TIME_STAMP):
                    # Data exists for the second time stamp, so break out
                    # The last no-data line was the previous line
                    nodata=False
                    break
                else:
                    # No data for this line either, so keep this timestamp
                    # and start the while(nodata) loop again
                    eventdate_keep=eventdate2

            # We now know the range of no-data values, so insert this range, latest first
	    noDataEvent['to'] = eventdate_keep
	    noDataEvent['from'] = eventdate1
	    noDataRanges.append(noDataEvent)
	    noDataCount += 1
            
    # Finally, prepend the data-not-yet-available date range, from the last day of data
    # until today's date.
    today=date.today()
    formattedToday=today.isoformat() 
    
    #Reread the last event
    rawFile.seek(0)
    lastRawLine=rawFile.readlines()[-1]
    lastRawDate=lastRawLine.split(" ")[1].split("T")[0]
    lastDataDateParts=lastRawDate.split("-")  # This is the last date
    #Create a new date object out of the string we get from the file.
    lastDataDatePlus1=date(int(lastDataDateParts[0]),int(lastDataDateParts[1]),int(lastDataDateParts[2]))
    #Now increment this date one day.
    lastDataDatePlus1+=timedelta(days=1)    
    #Now convert to a string
    lastDataDataP1String=lastDataDatePlus1.isoformat()

    noDataEvent = {}
    noDataEvent['to'] = formattedToday
    noDataEvent['from'] = lastDataDataP1String
    noDataRanges.append(noDataEvent)
    noDataCount += 1
    
    stationData['time_nodata'] = noDataRanges
    stationData['nodata_count'] = noDataCount
    rawFile.close
    return

numargv = len(sys.argv)
if numargv == 1:
    sys.exit("usage: create_mongodb.py scripps_dataset_name")
elif numargv == 2:
    dataSet = sys.argv[1]
else:
    sys.exit("Invalid number of parameters!")

projectDir = FINAL_PATH +dataSet

if(os.path.isdir(projectDir)):
# Open the JSON file that will contain the results
    outputPath = FINAL_PATH+dataSet + "_FILL.json"
    summaryData = {}
        
    summaryData['update_time'] = updateTime
    summaryData['data_source'] = dataSet
    summaryData['begin_date'] = beginDate
    summaryData['end_date'] = endDate
    summaryData['center_longitude'] = centerLng
    summaryData['center_latitude'] = centerLat
    summaryData['server_url'] = "http://" + serverName + "/daily_rdahmmexec/daily/" + dataSet
    summaryData['stateChangeNumTxtFile'] = stateChangeNumTxtFile
    summaryData['stateChangeNumJsInput'] = stateChangeNumJsInput
    summaryData['allStationInputName'] = allStationInputName
    summaryData['Filters'] = filters
    summaryData['video_url'] = ""

    stations = []
    stationCount = 0


for stationList in os.listdir(projectDir):
    stationPath = projectDir + "/" + stationList + "/"
    if (os.path.isdir(stationPath)):
        stationData = {}

        setStationId(stationList, stationData)
        setStationStartDate(stationPath, stationData)
        setStationRefLatLonHgt(stationPath, stationData)
        setStatusChanges(stationPath, stationData)
        setTimesNoData(stationPath, stationData)

        stations.append(stationData)
        stationCount += 1 

summaryData['stations'] = stations
summaryData['station_count'] = stationCount



##################################################################
# 
#	New Code to add data to Cassandra DB instead of saving to json
#
##################################################################
def convert_string_to_date(string):
#     return datetime.strptime(string, '%Y-%m-%d').now().date()
    return datetime.strptime(string, '%Y-%m-%d')


def checkDateForData(selectedDate,noDataDates):
    dataOnDate=True;
    #   Selected date is after the last no-data date.
    if selectedDate > convert_string_to_date(noDataDates[len(noDataDates) -1]['to']):
        dataOnDate=True;

    # Otherwise, check each no-data interval to see if the date falls within.
    else:
        for elem in noDataDates:
            startDate = convert_string_to_date(elem['from'])
            endDate = convert_string_to_date(elem['to'])

            if (startDate <= selectedDate) & (endDate >= selectedDate):
                dataOnDate=False
                break

    return dataOnDate



# Would be nice to throw an exception here 
def getPrecedingStateChange(selectedDate,statusChanges):
    
    stateLastDate = selectedDate
    latestPossibleDate=convert_string_to_date(statusChanges[len(statusChanges)-1]['date'])
    earliestPossibleDate=convert_string_to_date(statusChanges[0]['date'])

#     This should actually throw an erorr since there is no earlier state change.
    if(selectedDate <= earliestPossibleDate):
#         stateLastDate=earliestPossibleDate;
        stateLastDate=selectedDate
    
    elif (selectedDate >= latestPossibleDate):
        stateLastDate=latestPossibleDate
        
    else:
        for i, e in reversed(list(enumerate(statusChanges))):
            stateChangeDate1 = convert_string_to_date(statusChanges[i-1]['date'])
            stateChangeDate2 = convert_string_to_date(statusChanges[i]['date'])    
#             The last state change date to find is the one
#             on or before the curren date.
            if (selectedDate >= stateChangeDate1) & (selectedDate < stateChangeDate2):
                stateLastDate=stateChangeDate1
#                 Dates are in order, so we can stop
                break;
        
    return stateLastDate


def getStationState(date,gpsStation):
    #     theState=gpsStationState[0]  #This is the default.
    theState=0
    today=date
    lastMonth=today - timedelta(30)
    dayBefore=today - timedelta(1)

    statusChanges=gpsStation['status_changes']
    noDataDates=gpsStation['time_nodata']

    earliestDataDate = convert_string_to_date(gpsStation['start_date'])
    dataOnDate = checkDateForData(date,noDataDates)

    #     Hopefully this big if-else construction correctly captures all the case.
    #     It can be simplified later.

    #     Provided date is before any data available for that station, so state is light blue
    if (today < earliestDataDate):
    #         theState=gpsStationState[3]  #light blue
        theState=3
    #     print ("c1")


    #     Date falls within data range, there are no status changes, and data is available on date
    elif (today >= earliestDataDate) & (len(statusChanges)==0) & (dataOnDate==True):
    #         theState=gpsStationState[0]  # green
        theState=0
    #     print ("c2")

    #     Date falls within data range, there are no state changes, and no data on selected date.
    elif (today > earliestDataDate) & (len(statusChanges)==0) &(dataOnDate==False):
    #         theState=gpsStationState[3]  # light blue
        theState=3
    #     print ("c3")


    #     We have data on the date, but it proceeds the earliest state change
    elif (today < convert_string_to_date(statusChanges[0]['date'])):
    #         theState=gpsStationState[0]  # green
        theState=0
    #     print ("c4")

    #     See if the date falls within 1 day or 1 month of a state change.
    elif (len(statusChanges) > 0):
        #     print ("c5")
        #    Get nearest preceding state change date
        stateLastDate=getPrecedingStateChange(date,statusChanges)

    #         See if state change was yesterday.
        if(stateLastDate > dayBefore):
            #             theState=gpsStationState[1]  # red
            theState=1

    #         See if the station has changed state in between the
    #         selected date and 30 days prior to the selected date.
        elif (stateLastDate >= lastMonth) & (stateLastDate <= today):
            #       See if we have no data within 24 hours of the selected date.            
            if (dataOnDate == False):
                #       Data is missing within last 24 hours and state has changed within last 30 days.
                #                 theState=gpsStationState[4]  # blue
                theState=4
            else:
                #       We have data on the date and state has changed within a 30 day window.
                #                 theState=gpsStationState[2] # yellow
                theState=2
        #         No data is available on selected date for this station, and station is
        #         not within either the 1 day or 30 day window.
        elif(dataOnDate==False):
            #             theState=gpsStationState[3]  # light blue
            theState=3

    return theState



#########SAVE TO CASSANDRA DB###########
cluster = Cluster(contact_points=['127.0.0.1'], \
                  port=9042, \
                  load_balancing_policy= TokenAwarePolicy(DCAwareRoundRobinPolicy(local_dc='datacenter1')), \
                  default_retry_policy = RetryPolicy()
                 )
session = cluster.connect('system')
rows = cluster.metadata.keyspaces.keys()

KEYSPACE='GPS_'+dataSet
KEYSPACE=KEYSPACE.lower()

# Create database
keyspace_drop_stmt = session.prepare(
    "DROP KEYSPACE " + KEYSPACE
)

if KEYSPACE in rows:
    print ("dropping existing keyspace...")
    session.execute(keyspace_drop_stmt)

print("creating keyspace...")
session.execute("""
    CREATE KEYSPACE %s
    WITH replication = { 'class': 'SimpleStrategy', 'replication_factor': '1' }
    """ % KEYSPACE)

print("setting keyspace...")
session = cluster.connect(KEYSPACE)

print("creating table...")
session.execute("""
    CREATE TABLE time_series_stations (
        date timestamp,
        station_id text,
        status int,
        PRIMARY KEY (date, station_id)
    )
""")

session.execute("""
    CREATE TABLE meta_stations (
        station_id text PRIMARY KEY,
        pro_dir text,
        "AFile" text,
        "BFile" text,
        "InputFile" text,
        "RawInputFile" text,
        "SwfInputFile" text,
        "DygraphsInputFile" text,
        "LFile" text,
        "XPngFile" text,
        "YPngFile" text,
        "ZPngFile" text,
        "XTinyPngFile" text,
        "YTinyPngFile" text,
        "ZTinyPngFile" text,
        "PiFile" text,
        "QFile" text,
        "MaxValFile" text,
        "MinValFile" text,
        "RangeFile" text,
        "ModelFiles" text,
        "RefFile" text,
        start_date timestamp,
        lat double,
        lon double,
        height double
    );
""")

print "creating index..."
session.execute("""
    CREATE INDEX on %s .meta_stations (lat);
""" % KEYSPACE)

session.execute("""
    CREATE INDEX on %s .meta_stations (lon);
""" % KEYSPACE)


session.execute("""
    CREATE TABLE meta_network (
        network_name text PRIMARY KEY,
        update_time timestamp,
        data_source text,
        begin_date timestamp,
        end_date timestamp,
        center_longitude double,
        center_latitude double,
        server_url text,
        "stateChangeNumTxtFile" text,
        "stateChangeNumJsInput" text,
        "allStationInputName" text,
        "Filters" text,
        video_url text,
        station_count int
    );
""")


def get_station_element_as_list(station):
    dictionary = station.copy()
    dictionary.pop('status_changes', None)
    dictionary.pop('time_nodata', None)
    dictionary.pop('change_count', None)
    dictionary.pop('nodata_count', None)
    dictlist = []
    for key, value in sorted(dictionary.iteritems()):
#         print key, " : " ,value
        if (key == 'height') | (key == 'lat') | (key == 'long'):
            dictlist.append(float(value))
        elif key == 'start_date':
            dictlist.append(datetime.strptime(value, "%Y-%m-%d"))
        else:
            dictlist.append(value)

    return dictlist


start_date=datetime.strptime(beginDate, '%Y-%m-%d')
end_date=datetime.strptime(endDate, '%Y-%m-%d')



# http://stackoverflow.com/questions/1060279/iterating-through-a-range-of-dates-in-python
from datetime import timedelta, date
# http://www.datastax.com/dev/blog/4-simple-rules-when-using-the-datastax-drivers-for-cassandra
from cassandra.query import BatchStatement

timeseries_insert_stmt = session.prepare("""
        INSERT INTO Time_series_stations
        (date, station_id, status) 
        VALUES (?, ?, ?)
""")


meta_stations_insert_stmt = session.prepare("""
    INSERT INTO meta_stations
    ("AFile","BFile","DygraphsInputFile",
    "InputFile","LFile","MaxValFile", 
    "MinValFile", "ModelFiles", 
    "PiFile","QFile","RangeFile",
    "RawInputFile","RefFile",
    "SwfInputFile","XPngFile",
    "XTinyPngFile","YPngFile",
    "YTinyPngFile","ZPngFile",
    "ZTinyPngFile",height,
    station_id,lat,lon,
    pro_dir,start_date)
    VALUES (?, ?, ?, ?, ?, 
            ?, ?, ?, ?, ?, 
            ?, ?, ?, ?, ?, 
            ?, ?, ?, ?, ?, 
            ?, ?, ?, ?, ?, ?)
""")


meta_network_insert_stmt = session.prepare("""
    INSERT INTO meta_network
    (network_name, "Filters", "allStationInputName",
    begin_date, center_latitude, center_longitude, 
    data_source, end_date, server_url, 
    "stateChangeNumJsInput","stateChangeNumTxtFile", 
    station_count, update_time, video_url)
    VALUES (?, ?, ?, ?, ?, 
            ?, ?, ?, ?, ?, 
            ?, ?, ?, ?)   
""")

meta_network_list = [dataSet, \
                     summaryData['Filters'], \
                     summaryData['allStationInputName'], \
                     datetime.strptime(summaryData['begin_date'], "%Y-%m-%d"), \
                     float(summaryData['center_latitude']), \
                     float(summaryData['center_longitude']), \
                     summaryData['data_source'], \
                     datetime.strptime(summaryData['end_date'], "%Y-%m-%d"), \
                     summaryData['server_url'], \
                     summaryData['stateChangeNumJsInput'], \
                     summaryData['stateChangeNumTxtFile'], \
                     int(summaryData['station_count']), \
                     datetime.strptime(summaryData['update_time'], "%Y-%m-%dT%H:%M:%S"), \
                     summaryData['video_url']
                    ]

# Add data to meta_network
session.execute(meta_network_insert_stmt, meta_network_list)


# Add data to timeseries_stations
def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)
        
for single_date in daterange(start_date, end_date):
    
    for station in stations:
        session.execute(timeseries_insert_stmt, \
          [single_date, station['id'],getStationState(single_date, station)])


# Add data to meta_stations
batch = BatchStatement()

for station in stations:
#     print station['id']
    
    batch.add(meta_stations_insert_stmt, \
          get_station_element_as_list(station))
    
session.execute(batch)


# END SESSION
session.shutdown()