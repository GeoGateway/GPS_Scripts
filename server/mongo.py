# mongo.py
#==========================================================================
# Starts the service that on port 5000 that responds to HTTP get request.
# The service connects to mongodb which is running on port27017. 
# The service connects to UNR_SPLICE_timeseries_database.
# It returns the timeseries asssociated to all stations within the bounding box in a json file.

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

app = Flask(__name__)

database_name = 'UNR_SPLICE_timeseries_database'
app.config['MONGO_DBNAME'] = database_name
app.config['MONGO_URI'] = 'mongodb://localhost:27017/'+database_name

mongo = PyMongo(app)

@app.route('/gps', methods=['GET'])
def get_data_in_bounding_box():

  lat_min=request.args.get('lat_min')
  lat_max=request.args.get('lat_max')
  long_min=request.args.get('long_min')
  long_max=request.args.get('long_max')

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

  collections_stations = mongo.db.collections_stations
  collections_meta = mongo.db.collections_meta

  output = {}
  cursor_meta = collections_meta.find()

  for i in cursor_meta:
      i.pop('_id', None)
      i['network_station_count'] = i['station_count']
      i.pop('station_count', None)
      output = i;

  cursor_stations = collections_stations.find( {"lat" : {"$gte" : lat_min, "$lte" : lat_max}, \
                            "long" : {"$gte" : long_min, "$lte" : long_max}})

  list_station = []
  output['station_count'] = cursor_stations.count()

  for i in cursor_stations:
      i.pop('_id', None)
      list_station.append(i)

  output["stations"] = list_station
  # http://stackoverflow.com/questions/19877903/using-mongo-with-flask-and-python
  return json.dumps(output, sort_keys=True, indent=2)

if __name__ == '__main__':
    # app.run(debug=True)
    app.run(host= '0.0.0.0', threaded=True, debug=False)