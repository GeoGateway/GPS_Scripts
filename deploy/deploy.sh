#!/bin/bash

echo 
echo
echo "***********************************"
echo "Installing prerequisites"
echo "***********************************"
sudo apt-get update
sudo apt-get install python
sudo apt-get python-pip
pip install --upgrade pip
pip install virtualenv


echo 
echo
echo "***********************************"
echo "Creating virtual environment"
echo "***********************************"
virtualenv deploy --no-site-packages

echo 
echo
echo "***********************************"
echo "Enabling virtual environment"
echo "***********************************"
source deploy/bin/activate

echo "***********************************"
echo "Installing pre requisites"
echo "***********************************"
pip install -r requirements-deploy.txt 

echo 
echo
echo "***********************************"
echo "Installing simplejson on server"
echo "***********************************"
ansible backend -i inventory -m raw -a "sudo yum install -y python-simplejson"


echo "***********************************"
echo "Deploying code on backend"
echo "***********************************"
ansible-playbook -i inventory main.yml
