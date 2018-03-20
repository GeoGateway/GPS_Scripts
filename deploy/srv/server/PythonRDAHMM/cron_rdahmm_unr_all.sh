#!/bin/bash

# ~/GPS_Scripts/PythonRDAHMM/cron_rdahmm_unr.py
time python PythonRDAHMM/unr_ingest_single.py igs08
echo "IGS08 ingesting complete"

# Need to use python instead of pypy for some obscure bug workaround
time python PythonRDAHMM/unr_ingest_single.py fid
echo "FID ingestion complete"

time python PythonRDAHMM/unr_splice.py
echo "Splicing complete"

time python PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE
echo "RDAHMM evaluation complete"

#python PythonRDAHMM/create_cassandra_db.py UNR_SPLICE
#time python PythonRDAHMM/create_mongodb.py UNR_SPLICE
time python PythonRDAHMM/mongodb_load_delta.py UNR_SPLICE
echo "MongoDB update complete"

# /home/yuma/PythonRDAHMM/cron_rdahmm_unr.py
# /home/yuma/PythonRDAHMM/unr_splice.py
# /home/yuma/PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_xmls.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_jsons.py UNR_SPLICE
