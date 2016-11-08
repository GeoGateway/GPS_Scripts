# mongo.py
#==========================================================================
# Starts the service on port 5000 that responds to HTTP get request.
# The service connects to mongodb which is running on port27017. 
# The service connects to UNR_SPLICE_timeseries_database.
# It returns the status asssociated to all stations within the bounding box in a json file.

# usage: to start the service
#   python mongo.py

# eg: request to service: http://localhost:5000/gps?lat_min=18.005611&lat_max=48.987386&lon_min=-124.626080&lon_max=-62.361014
#     output: Is a json file containinf information
#             of all stations in the bounding box

#==========================================================================


from flask import Flask
from flask import jsonify
from flask import request
from flask_pymongo import PyMongo
from flask.ext.jsonpify import jsonify

import json
from datetime import datetime

app = Flask(__name__)

database_name = 'GPS_UNR_SPLICE'
app.config['MONGO_DBNAME'] = database_name
app.config['MONGO_URI'] = 'mongodb://localhost:27017/'+database_name

mongo = PyMongo(app)

@app.route('/'+ database_name+'/time_series', methods=['GET'])
def get_data_in_bounding_box():

	lat_min=request.args.get('lat_min')
	lat_max=request.args.get('lat_max')
	lon_min=request.args.get('lon_min')
	lon_max=request.args.get('lon_max')
	year=request.args.get('year')
	month=request.args.get('month')
	day=request.args.get('day')

	if not lat_min:
		lat_min=-89.998914447
	else:
		lat_min=float(lat_min)   

	if not lat_max: 
		lat_max=83.64323665
	else:
		lat_max=float(lat_max)   

	if not lon_min:
		lon_min=-359.998835383    
	else:
		lon_min=float(lon_min)   

	if not lon_max:
		lon_max=-0.037155605
	else:
		lon_max=float(lon_max)

	today = datetime.today()

	if not year:
		year=str(today.year)

	if not month:
		month=str(today.month)

	if not day:
		day=str(today.day)


	collections_time_series_stations = mongo.db.collections_time_series_stations
	# collections_meta = mongo.db.collections_meta

	output = {}
	# cursor_meta = collections_meta.find()

	# for i in cursor_meta:
	#   i.pop('_id', None)
	#   i['network_station_count'] = i['station_count']
	#   i.pop('station_count', None)
	#   output = i

	status_to_find = 'status.' + year + '.' + month + '.' + day
	cursor_stations = collections_time_series_stations.find({ "loc" : \
	                                                  { "$geoWithin" : 
	                                                    { "$box" : [ [lon_min, lat_min],\
	                                                               [lon_max, lat_max]] \
	                                                    } }, \
	                                                 status_to_find : {'$exists': 1} \
	                                                }, {'station_id': 1, 'loc': 1, status_to_find:1})

	list_station = []
	if cursor_stations:
	    for station in cursor_stations:
	        try:
	            res = {'station_id' : station['station_id'],
	                   'status' : station['status'][year][month][day],
	                   'lat' : station['loc'][1],
	                   'lon' : station['loc'][0]
	                  }
	            list_station.append(res)
	        except:
	            continue

	output['station_count'] = len(list_station)
	output["stations"] = list_station
	# http://stackoverflow.com/questions/19877903/using-mongo-with-flask-and-python
	# return json.dumps(output, sort_keys=True, indent=2)
	return jsonify(output)

########### GET NETWORK META
@app.route('/'+ database_name+ '/network_meta', methods=['GET'])
def get_network_meta_data():

	collections_meta_network = mongo.db.collections_meta_network

	output = {}
	cursor_meta_network = collections_meta_network.find()

	for i in cursor_meta_network:
	    i.pop('_id', None)
	    i['network_station_count'] = i['station_count']
	    i.pop('station_count', None)
    	output = i

	return jsonify(output)

#########GET STATION META
@app.route('/'+ database_name+'/station_meta', methods=['GET'])
def get_station_meta_data():

	station_id_to_find=request.args.get('station_id_to_find')
	collections_meta_stations = mongo.db.collections_meta_stations

	output = {}
	cursor_meta_stations = collections_meta_stations.find( \
	                                        { '_id': station_id_to_find })

	for i in cursor_meta_stations:
	    output=i

	return jsonify(output)


if __name__ == '__main__':
    # app.run(debug=True)
    app.run(host= '0.0.0.0', threaded=True, debug=False)