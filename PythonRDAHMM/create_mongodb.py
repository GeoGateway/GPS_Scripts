# #!/usr/local/bin/python
#==========================================================================
# Script to load the final json file into MongoDB.

# usage:
	 # python create_mongodb scripps_dataset_name
	 # eg: python create_mongodb.py UNR_SPLICE

# output:
# 	MongoDB database named gps_timeseries_database which contains
# 	two collections-
# 	collections_meta: Data commom to all stations
# 	collections_stations: Timeseries for each station
#	The database is indexed based on latitude and longitude
#===========================================================================
import sys
from properties import properties
import json
import pymongo
from pymongo import MongoClient

numargv = len(sys.argv)
if numargv == 1:
    sys.exit("usage: create_mongodb.py scripps_dataset_name")
elif numargv == 2:
    dataset = sys.argv[1]
else:
    sys.exit("Invalid number of parameters!")

# Path to json
json_file_name = dataset + '_FILL.json'
json_file = properties('eval_path') + json_file_name 

# Read the json
with open(json_file) as data_file:    
    data = json.load(data_file)

# Convert lat & long to float
for dictionary in data['stations']:
    dictionary['lat'] = float(dictionary['lat'])
    dictionary['long'] = float(dictionary['long'])

# Connect to MongoDB
client = MongoClient('localhost', 27017)

# Delete existing database
database_name = dataset+'_timeseries_database'
client.drop_database(database_name)

# Create new database
db =client[database_name]

# Create 2 collections
# for meta
collections_meta = db.collections_meta
# for stations
collections_stations = db.collections_stations 

# Add meta data
meta= {}
meta['update_time'] = data['update_time']
meta['data_source'] = data['data_source']
meta['begin_date'] = data['begin_date']
meta['end_date'] = data['end_date']
meta['center_longitude'] = data['center_longitude']
meta['center_latitude'] = data['center_latitude']
meta['server_url'] = data['server_url']
meta['stateChangeNumTxtFile'] = data['stateChangeNumTxtFile']
meta['stateChangeNumJsInput'] = data['stateChangeNumJsInput']
meta['allStationInputName'] = data['allStationInputName']
meta['Filters'] = data['Filters']
meta['video_url'] = data['video_url']
meta['station_count'] = data['station_count']

collections_meta.insert_one(meta)

# Add station data
collections_stations.insert_many(data['stations'])

# Create index based on latitude, longitude
db.collections_stations.create_index( [("lat", pymongo.ASCENDING), \
                                     ("lon", pymongo.ASCENDING)] )

# Close connection to database
client.close()