#!/bin/bash

echo "***********************************"
echo "Installing pre requisites"
echo "***********************************"
pip install -r requirements.txt 

echo 
echo
echo "***********************************"
echo "Not Starting MongoDB"
echo "***********************************"
#sudo service mongod start

echo 
echo
echo "***********************************"
echo "Changing permissions to Executables"
echo "***********************************"
cd /root/GPS_Scripts/RDAHMM
chmod -R 777 rdahmm3

echo 
echo
echo "***********************************"
echo "Adding data to MongoDB"
echo "***********************************"
python PythonRDAHMM/create_mongodb.py UNR_SPLICE


echo 
echo
echo "***********************************"
echo "Creating virtual environment"
echo "***********************************"
virtualenv flask-mongo --no-site-packages

echo 
echo
echo "***********************************"
echo "Enabling virtual environment"
echo "***********************************"
source flask-mongo/bin/activate
# alias activate=". flask/bin/activate"

echo 
echo
echo "***********************************"
echo "Installing server prerequisites"
echo "***********************************"
pip install -r requirements-server-mongo.txt

echo 
echo
echo "***********************************"
echo "Starting server"
echo "***********************************"
echo "Check if the service has started by typing or pasting following command in a new terminal"
echo curl "http://localhost:5000/gps?lat_min=18.005611&lat_max=48.987386&long_min=-124.626080&long_max=-62.361014"
echo "Should display a json file giving data for USA"
python start_server_mongo.py
