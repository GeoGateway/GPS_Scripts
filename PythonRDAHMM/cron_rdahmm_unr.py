#!/usr/local/bin/python
#==========================================================================
# Ingest, and execute rdahmm evaluation for UNR datasets 
# Set up a cron job to run nightly
#
# usage: cron_rdahmm_unr.py
#
#===========================================================================
import os, subprocess, sys
from threading import Thread
from properties import properties

unr_cmd = properties('script_path') + "/unr_ingest_single.py"
eval_cmd = properties('script_path') + "/rdahmm_eval_single.py"
xml_cmd = properties('script_path') + "/create_summary_xmls.py"
json_cmd = properties('script_path') + "/create_summary_jsons.py"

class ThreadJob(Thread):

    def __init__(self, dataset):
        Thread.__init__(self)
        self.source = dataset
        self.dataset = "UNR_" + dataset.upper()

    def run(self):
	# ingest a given dataset: igs08 | fid
        print "+++Starting process UNR ", self.source, " ..."
        cmd = unr_cmd
	print "Command: " + cmd
        p = subprocess.Popen(["sudo", cmd, self.source], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        if p.returncode != 0:
            print p.stderr  
        print "+++Finished process UNR ", self.source
     
        # run rdahmm evaluation on the corresponding dataset 
        print "+++Starting process ", self.dataset, " ..."
        cmd = eval_cmd
	print "command: " + cmd
        #cmd = "echo"
        p = subprocess.Popen(["sudo", cmd, self.dataset], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        if p.returncode != 0:
            print p.stderr        
        print "+++Finished process ", self.dataset

        # # create summary xml on the corresponding dataset 
        # print "+++creating summary xml for  ", self.dataset, " ..."
        # cmd = xml_cmd
        # #cmd = "echo"
        # p = subprocess.Popen([cmd, self.dataset], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # (stdout, stderr) = p.communicate()
        # if p.returncode != 0:
        #     print p.stderr        
        # print "+++Finished creating summary xml for  ", self.dataset

        # # create summary json on the corresponding dataset 
        # print "+++creating summary json for  ", self.dataset, " ..."
        # cmd = json_cmd
        # #cmd = "echo"
        # p = subprocess.Popen([cmd, self.dataset], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # (stdout, stderr) = p.communicate()
        # if p.returncode != 0:
        #     print p.stderr        
        # print "+++Finished creating summary json for  ", self.dataset

for dataset in ['fid', 'igs08']:
    t = ThreadJob(dataset)
    t.start()
