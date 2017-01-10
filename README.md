# GPS_Scripts

Summary
------------------
Server to respond to the frontend- [Geo-Gateway](http://geo-gateway.org/main.html). The server is designed using micro service architecture to handle requests from the frontend. The current implementation can process 60 million data points in under 1.5 seconds using MongoDB. This repository is designed to configure and deploy code on the backend.

Prerequisite
------------------
The deployment is divided into 2 parts-  
1. Deploy machine- Where you will clone this repository. It can be your PC or a VM on the cloud. A VM on cloud is preferred because deployement process will be faster.  
2. Host machine - It will host the backend server. Make sure that it is an instance of 64 bit, 4 GB RAM or more, 100GB memory or more, [CentOS 5.4](https://www.centos.org/) or less with Python 2.4 or less installed. An image can be found in community AMI of N.Virginia region of Amazon AWS EC2 (ami-7ea24a17).   

The deploy machine first should have following installed-  
1. [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)    
2. [Python2.7](https://www.python.org/download/releases/2.7/)  
3. [Pip] (https://pypi.python.org/pypi/pip)

Please make sure the following steps are carried out-      

1. You should have <b>root access</b> of host. 
1. You should be able to <b>ssh</b> from deploy machine to host. This can be done by copying the public SSH key of deploy machine to authorized_keys of host.  

  Eg: If you are using amazon ec2 instance for host you can connect using following command-  
  ```shell
    ssh -i "~/.ssh/geo-gateway.pem" root@ec2-54-86-93-172.compute-1.amazonaws.com
  ```  
1.  For Ansible to run please setup a **ssh agent** on deploy machine and add the public key.    

  ```shell  
    ssh-agent bash  
    ssh-add ~/.ssh/id_rsa_gpsSetup
  ```  


Installation
------------------

1. Clone this repository on deploy machine.   

  ```shell
  git clone https://github.com/GeoGateway/GPS_Scripts.git
  ```  
2. Change directory to the deploy folder.     

  ```shell
  cd GPS_Scripts/deploy
  ```  
3. Change the backend I.P. address in the inventory to specify on which machine to install.  
  [backend]  
  ~~54.86.93.172~~  
  54.157.16.80  
4. Set the destination in [backend-code.yml](deploy/backend-code.yml) where you want to install the server. It should be in a partition where you have 100 GB space.  
~~server_dest: /mnt/data~~  
server_dest: /root  
5. Execute the deploy.sh as follows- 

  ```shell   
  source deploy.sh  
  ```  
6. Ssh into the backend by  

  ```shell
    ssh -i "~/.ssh/geo-gateway.pem" root@ec2-54-86-93-172.compute-1.amazonaws.com
  ```  
7. Change directory to server  

  ```shell
  cd server
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


