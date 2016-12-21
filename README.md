# GPS_Scripts

Summary
------------------
Server to respond to the frontend- [Cassandra DB](http://cassandra.apache.org/). The server is designed using micro service architecture to handle requests from the frontend. The current implementation can process 80 million data points in under 9 seconds using Cassandra DB. This repository is designed to store data in cassandraDB and start a flask service listening at port 5000.  

Prerequisite
------------------
Please make sure the following steps are carried out before installation   

1. If using a cloud instance, you should be able to ssh into your cloud virtual machine and have <b>root access</b>.  

  Eg: If using amazon ec2, connect using following command
  ```
    ssh -i "*****.pem" ubuntu@ec2-publicIP.compute-1.amazonaws.com
  ```
2. [Python2.7](https://www.python.org/download/releases/2.7/) installed in your machine.  
3. Git  
4. Pip    
5. [Cassandra DB](http://cassandra.apache.org/) installed and executing.  


Backend is an instance of 64 bit, 4 GB RAM, 100GB memory [CentOS 5.4](https://www.centos.org/) or less with Python 2.4 or less installed. An image can be found in community AMI of N.Virginia region of Amazon AWS EC2 (ami-7ea24a17).  

Installation
------------------

1. Clone this repository  

  ```shell
  git clone https://github.com/GeoGateway/GPS_Scripts.git
  ```  
2. Change directory to the GPS_Scripts folder   

  ```shell
  cd GPS_Scripts
  ```  
3.  Source the start_server_cassandra.sh file. It will create a virtual environment, install the prerequites and start the server on port 5000.  

  ```shell
  source start_server_cassandra.sh
  ```  
4. To create database execute the [create_cassandra_db.py](PythonRDAHMM/create_cassandra_db.py).  
```shell
 python create_cassandra_db.py UNR_SPLICE
```  

Author Information
------------------

Hrushikesh Dhumal (hrushikesh.dhumal@gmail.com)


