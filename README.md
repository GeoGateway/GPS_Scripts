# GPS_Scripts

Summary
------------------
Server to respond to the frontend- [Geo-Gateway](http://geo-gateway.org/main.html). The server is designed using micro service architecture to handle requests from the frontend. The current implementation can process 8 million data points in under 1.5 seconds using MongoDB. This repository is designed to configure and deploy code on the backend.

Prerequisite
------------------
The deployment is divided into 2 parts-  
1. A machine where you will clone this repository- it can be your machine or a machine on the cloud. A machine on cloud is preferred because deployement process will be faster.  
2. Backend machine.  

For first part you will need-  
1. [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)    
2. [Python2.7](https://www.python.org/download/releases/2.7/)  
3. [Pip] (https://pypi.python.org/pypi/pip)

For the second part (backend) please make sure the following steps are carried out-      

1. If using a cloud instance, you should be able to ssh into your cloud virtual machine and have <b>root access</b>.  

  Eg: If you are using amazon ec2 you can connect using following command
  ```shell
    ssh -i "~/.ssh/geo-gateway.pem" root@ec2-54-86-93-172.compute-1.amazonaws.com
  ```
2. Backend is an instance of 64 bit, 4 GB RAM, 100GB memory [CentOS 5.4](https://www.centos.org/) or less with Python 2.4 or less installed. An image can be found in community AMI of N.Virginia region of Amazon AWS EC2 (ami-7ea24a17).  

Installation
------------------

1. Clone this repository  

  ```shell
  git clone https://github.com/GeoGateway/GPS_Scripts.git
  ```  
2. Change directory to the deploy folder   

  ```shell
  cd GPS_Scripts/deploy
  ```  
3. Change the backend I.P. address in the inventory.  
  [backend]  
  ~~54.86.93.172~~  
  54.164.84.82  
4. Execute the deploy.sh as follows- 

  ```shell   
  source deploy.sh  
  ```  
  This will install [Ansible](https://www.ansible.com/) to automate the deployement process. On backend, it will install [MongoDB](https://www.mongodb.com/), [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git), [Python2.7](https://www.python.org/download/releases/2.7/), [laplack](http://www.netlib.org/lapack/), [tmux](https://tmux.github.io/). It takes some time for the code to execute so please wait.  
5. Ssh into the backend by  

  ```shell
    ssh -i "~/.ssh/geo-gateway.pem" root@ec2-54-86-93-172.compute-1.amazonaws.com
  ```
6. Clone this repository  

  ```shell  
  git clone https://github.com/GeoGateway/GPS_Scripts.git
  ```
7. Change directory to GPS_Scripts  

  ```shell
  cd GPS_Scripts/deploy/srv/GPS_Scripts/
  ```  
8. Start a tmux session by  

  ```shell  
  tmux
  ```  
9. Execute the start_server_mongo.sh file. It will start the server, however there is no database generated. 

 ```shell 
  source start_server_mongo.sh
  ```  
10. Detach from the tmux session and let the server run continously by pressing together Control, B and D keys.  
11. To generate the database first add the model files to the RDAHMM/Model folder  
  For eg:  
  
  ```shell
  rsync -au --progress -e "ssh -i ~/.ssh/geo-gateway.pem" \
       /media/hru/Data/_Active_Projects/IS/RDAHMM/Model/UNR_SPLICE.tar.gz \
       root@ec2-54-85-163-109.compute-1.amazonaws.com:
  ```  
12. Start a new tmux session by  

  ```shell
  tmux
  ```  
13.  The data can be added to the database by  

  ```shell
  sh PythonRDAHMM/cron_rdahmm_unr_all.sh
  ```  
14. Detach from the tmux session by pressing together Control, B and D key.  

 
Author Information
------------------

Hrushikesh Dhumal (hrushikesh.dhumal@gmail.com)


