#==========================================================================
# Ingest JPL scripps datasets downloaded into proper directories and databases
#
# usage: python scripps_ingest_jpl.py
#
#===========================================================================
from ftplib import FTP
from properties import properties
import os

data_path = properties('data_path')
datadir = data_path + "JPL"
print "Wrinting to: ", datadir

def make_sure_path_exists(path): 
    '''

    Function to  make a new directory structure if it doesnt exists.
    Created by following: 
    # http://stackoverflow.com/questions/273192/how-to-check-if-a-directory-exists-and-create-it-if-necessary
    Args: 
        path (str): Directory path

    Returns: Creates directories if they dont exist

    '''
    if not os.path.exists(path):
        os.makedirs(path)

def connect(server='sideshow.jpl.nasa.gov'):
    '''
    Function to establish ftp connection to server
    
    Args: 
        server (str): url to server. 
            By default it is the JPL server url

    Returns: 
        ftp connection
    '''

    ftp = FTP(server)     # connect to host, default port
    ftp.login()                     # user anonymous, passwd anonymous@
    return ftp

def complete_data_downloaded(web_list, datadir):
    '''
    Function to check if all the files are downloaded. 
    It was implemented as the connection to the Nasa server
    gets timeout and drops many times.
    
    Args: 
        web_list (list): list containing the files that need to downloaded
        datadir (str): the directory where files are being downloaded
        
    Returns:
        Boolean indicating if download is complete

    '''
    
#     http://stackoverflow.com/questions/11697709/comparing-two-lists-in-python
    file_list = os.listdir(datadir)
    common_files = list(set(file_list) & set(filenames))    
    if  len(common_files)== len(web_list):
        return True
    else:
        print("Missing " + str(len(web_list) - len(common_files)) + " files")
        return False

make_sure_path_exists(datadir)

ftp = connect()
ftp.cwd('pub/usrs/mbh/point') # Change working directory to time series folder


filenames = ftp.nlst() # get filenames within the directory
    
complete_data_downloaded(filenames, datadir)

connect_count = 0
while not complete_data_downloaded(filenames, datadir):
    try:
        for filename in filenames:
            local_filename = os.path.join(datadir, filename)
            if not os.path.exists(local_filename):
                print("Writing to: "+ local_filename)
                file = open(local_filename, 'wb')
                ftp.retrbinary('RETR '+ filename, file.write)

                file.close() 
            else:
                print(local_filename+" exists")
    except:
        ftp = connect()
        ftp.cwd('pub/usrs/mbh/point') 
        connect_count +=1

ftp.quit()