#!/bin/bash

# ~/GPS_Scripts/PythonRDAHMM/cron_rdahmm_unr.py
python PythonRDAHMM/unr_ingest_single.py igs08
python PythonRDAHMM/unr_ingest_single.py fid

python PythonRDAHMM/unr_splice.py

python PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE

#python PythonRDAHMM/create_cassandra_db.py UNR_SPLICE
python PythonRDAHMM/create_mongodb.py UNR_SPLICE


# /home/yuma/PythonRDAHMM/cron_rdahmm_unr.py
# /home/yuma/PythonRDAHMM/unr_splice.py
# /home/yuma/PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_xmls.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_jsons.py UNR_SPLICE