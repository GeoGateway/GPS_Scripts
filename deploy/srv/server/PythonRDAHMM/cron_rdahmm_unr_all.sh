#!/bin/bash

# ~/GPS_Scripts/PythonRDAHMM/cron_rdahmm_unr.py
pypy PythonRDAHMM/unr_ingest_single.py igs08

# Need to use python instead of pypy for some obscure bug workaround
python PythonRDAHMM/unr_ingest_single.py fid

pypy PythonRDAHMM/unr_splice.py

pypy PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE

#python PythonRDAHMM/create_cassandra_db.py UNR_SPLICE
pypy PythonRDAHMM/create_mongodb.py UNR_SPLICE


# /home/yuma/PythonRDAHMM/cron_rdahmm_unr.py
# /home/yuma/PythonRDAHMM/unr_splice.py
# /home/yuma/PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_xmls.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_jsons.py UNR_SPLICE
