# cassandra.py
#==========================================================================
# Starts the service on port 5000 that responds to HTTP get request.
# The service connects to mongodb which is running on port27017. 
# The service connects to UNR_SPLICE_timeseries_database.
# It returns the status asssociated to all stations within the bounding box in a json file.

# usage: to start the service
#   python cassandra.py

# eg: request to service: http://localhost:5000/gps?lat_min=18.005611&lat_max=48.987386&lon_min=-124.626080&lon_max=-62.361014
#     output: Is a json file containing information
#             of all stations in the bounding box

#==========================================================================


from flask import Flask
from flask import jsonify
from flask import request
from flask.ext.jsonpify import jsonify

from cassandra.cluster import Cluster
from cassandra.query import dict_factory
from cassandra.query import named_tuple_factory

import json
from datetime import datetime

app = Flask(__name__)

KEYSPACE = 'gps_unr_splice'
database_name = 'GPS_UNR_SPLICE'

# Establish connection to cassandra
cluster = Cluster()

# Initialize all queries to none
session_bounding_box = None
session_network = None
session_station = None
meta_stations_select_stmt = None
timeseries_select_stmt = None
meta_network_select_stmt = None
meta_stations_select_stmt = None

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
		year=today.year
	else:
		year=int(year)

	if not month:
		month=today.month
	else:
		month=int(month)

	if not day:
		day=today.day
	else:
		day=int(day)


	session_bounding_box = cluster.connect(KEYSPACE)
	session_bounding_box.row_factory = named_tuple_factory

	meta_stations_select_stmt = session_bounding_box.prepare("""
	    select station_id, lat, lon from meta_stations 
	    where lon >= ? and lon <= ? and lat >= ? and lat <= ?
	    allow FILTERING
	""")

	selected_stations = session_bounding_box.execute(meta_stations_select_stmt, \
	                        [lon_min, lon_max, lat_min, lat_max])


	futures = []
	timeseries_select_stmt = session_bounding_box.prepare("""
	    SELECT status  
	    FROM time_series_stations
	    WHERE date = ? AND station_id = ?
	""")

	list_station = list()
	for e in list(selected_stations):
	    res ={ 'station_id' : str(e.station_id), \
	         'lat' : str(e.lat), \
	         'lon' : str(e.lon) \
	         }
	    list_station.append(res)
	    futures.append(session_bounding_box.execute_async(timeseries_select_stmt, \
	                              [datetime(year,month,day), str(e.station_id)]))


	# wait for them to complete and use the results

	for idx, future in enumerate(futures):
	    record = future.result()
	    
	#     print list(record)
	    t = list_station[idx]
	    t['status']= str(list(record)[0][0])
	    
	output = {}

	output['station_count'] = len(list_station)
	output["stations"] = list_station
	# http://stackoverflow.com/questions/19877903/using-mongo-with-flask-and-python
	# return json.dumps(output, sort_keys=True, indent=2)
	return jsonify(output)

# # ########### GET NETWORK META
@app.route('/'+ database_name+ '/network_meta', methods=['GET'])
def get_network_meta_data():

	session_network = cluster.connect(KEYSPACE)
	session_network.row_factory = dict_factory

	meta_network_select_stmt = session_network.prepare("""
	    select * from meta_network 
	""")

	nw = session_network.execute(meta_network_select_stmt)

	return jsonify(nw[0])

# #########GET STATION META
@app.route('/'+ database_name+'/station_meta', methods=['GET'])
def get_station_meta_data():

	session_station = cluster.connect(KEYSPACE)
	session_station.row_factory = dict_factory

	station_id_to_find=request.args.get('station_id_to_find')
	
	meta_stations_select_stmt = session_station.prepare("""
	    select * from meta_stations 
	    where station_id = ?
	""")

	selected_station = session_station.execute(meta_stations_select_stmt, \
	                        [station_id_to_find])	

	return jsonify(selected_station[0])


if __name__ == '__main__':
    # app.run(debug=True)
    app.run(host= '0.0.0.0', threaded=True, debug=False)