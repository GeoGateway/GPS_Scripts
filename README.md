# GPS_Scripts

Summary
------------------
Server to respond to the [Geo-Gateway](http://geo-gateway.org/main.html) request to get data. Originally, entire data was being loaded to client side slowing down the system. This repository tries to speed up the process by storing the data in MongoDB and making it easy to query and fetch only relevant data.  
This role is not yet operational.  

Prerequisite
------------------

Please make sure the following steps are carried out before installation   

1. If using a cloud instance, you should be able to ssh into your cloud virtual machine and have <b>root access</b>.  

  Eg: If you are using amazon ec2 you can connect using following command
  ```
    ssh -i "*****.pem" ubuntu@ec2-publicIP.compute-1.amazonaws.com
  ```
2. [Python2.7](https://www.python.org/download/releases/2.7/) installed in your machine.  
3. [MongoDB](https://www.mongodb.com/) is installed. MongoDB can be installed by following the steps from this [tutorial](https://www.howtoforge.com/tutorial/install-mongodb-on-ubuntu-14.04/).
  ```
    sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 7F0CEB10  
    echo "deb http://repo.mongodb.org/apt/ubuntu "$(lsb_release -sc)"/mongodb-org/3.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-3.0.list  
    sudo apt-get update  
    sudo apt-get install -y mongodb-org
  ```

4. Start the MongoDB service by-  
  ```
   sudo service mongod start 
   
   OR  
   
   mongod  
  ```

Installation
------------------

1. Clone this repository  
   ```
   git clone https://github.com/GeoGateway/GPS_Scripts.git
   ```  
2. Change directory to this repository  
  ```
  cd GPS_Scripts
  ```  
3. Execute the runAll.sh file. It will load the json data to mongoDB and print out the instructions that need to be followed to start the server.  
  ```
  sudo runAll.sh
  ```

Author Information
------------------

Hrushikesh Dhumal (hrushikesh.dhumal@gmail.com)
