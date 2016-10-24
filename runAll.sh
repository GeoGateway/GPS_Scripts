#!/bin/bash

echo "***********************************"
echo "Installing pre requisites"
echo "***********************************"
pip install -r requirements.txt 

echo 
echo
echo "***********************************"
echo "Installing MongoDB"
echo "***********************************"
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 7F0CEB10  
echo "deb http://repo.mongodb.org/apt/ubuntu "$(lsb_release -sc)"/mongodb-org/3.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-3.0.list  
sudo apt-get update  
sudo apt-get install -y mongodb-org

echo 
echo
echo "***********************************"
echo "Starting MongoDB"
echo "***********************************"
sudo service mongod start

echo 
echo
echo "***********************************"
echo "Adding data to MongoDB"
echo "***********************************"
python PythonRDAHMM/create_mongodb.py UNR_SPLICE

echo 
echo
echo "***********************************"
echo "Next steps"
echo "***********************************"
echo "Change directory to server and execute start_server.sh"
echo "cd server"
echo "echo sudo source start_server.sh"
