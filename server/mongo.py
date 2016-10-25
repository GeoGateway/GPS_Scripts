# mongo.py
#==========================================================================
# Starts the service on port 5000 that responds to HTTP get request.
# The service connects to mongodb which is running on port27017. 
# The service connects to UNR_SPLICE_timeseries_database.
# It returns the status asssociated to all stations within the bounding box in a json file.

# usage: to start the service
#   python mongo.py

# eg: request to service: http://localhost:5000/gps?lat_min=18.005611&lat_max=48.987386&long_min=-124.626080&long_max=-62.361014
#     output: Is a json file containinf information
#             of all stations in the bounding box

#==========================================================================


from flask import Flask
from flask import jsonify
from flask import request
from flask_pymongo import PyMongo
import json
from datetime import datetime

app = Flask(__name__)

database_name = 'GPS_UNR_SPLICE'
app.config['MONGO_DBNAME'] = database_name
app.config['MONGO_URI'] = 'mongodb://localhost:27017/'+database_name

mongo = PyMongo(app)

@app.route('/gps', methods=['GET'])
def get_data_in_bounding_box():

	lat_min=request.args.get('lat_min')
	lat_max=request.args.get('lat_max')
	long_min=request.args.get('long_min')
	long_max=request.args.get('long_max')
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

	if not long_min:
		long_min=-359.998835383    
	else:
		long_min=float(long_min)   

	if not long_max:
		long_max=-0.037155605
	else:
		long_max=float(long_max)

	today = datetime.today()

	if not year:
		year=str(today.year)

	if not month:
		month=str(today.month)

	if not day:
		day=str(today.day)


	collections_stations = mongo.db.collections_stations
	# collections_meta = mongo.db.collections_meta

	output = {}
	# cursor_meta = collections_meta.find()

	# for i in cursor_meta:
	#   i.pop('_id', None)
	#   i['network_station_count'] = i['station_count']
	#   i.pop('station_count', None)
	#   output = i

	cursor_stations = collections_stations.find( { "loc" : 
                                                  { "$geoWithin" :
                                                    { "$box" : [ [long_min, lat_min],
                                                               [long_max, lat_max]]
                                             } } } )
	list_station = []
	if cursor_stations:
	    for i in cursor_stations:
	        i.pop('_id', None)
	        list_station.append(i)

	output['station_count'] = cursor_stations.count()

	for i in cursor_stations:
		i.pop('_id', None)
		list_station.append(i)

	list_status = []

	for station in list_station:
    
	    try:
	        res = {'id' : station['id'],
	               'status' : station['status'][year][month][day]
	              }
	        list_status.append(res)
	    except:
	        continue

	output["status"] = list_status
	# http://stackoverflow.com/questions/19877903/using-mongo-with-flask-and-python
	return json.dumps(output, sort_keys=True, indent=2)

if __name__ == '__main__':
    # app.run(debug=True)
    app.run(host= '0.0.0.0', threaded=True, debug=False)