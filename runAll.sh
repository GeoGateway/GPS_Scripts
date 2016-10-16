#!/bin/bash

echo "***********************************"
echo "Installing pre requisites"
echo "***********************************"
pip install -r requirements.txt 

echo 
echo
echo "***********************************"
echo "Creating MongoDB"
echo "***********************************"
python PythonRDAHMM/create_mongodb.py UNR_SPLICE

echo 
echo
echo "***********************************"
echo "Next steps"
echo "***********************************"
echo "Change directory to server and execute start_server.sh"
echo "echo server"
echo "echo source start_server.sh"