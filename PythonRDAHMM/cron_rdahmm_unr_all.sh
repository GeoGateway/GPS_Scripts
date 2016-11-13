#!/bin/bash

# ~/GPS_Scripts/PythonRDAHMM/cron_rdahmm_unr.py
python PythonRDAHMM/unr_ingest_single.py igs08
python PythonRDAHMM/unr_ingest_single.py fid

python PythonRDAHMM/rdahmm_eval_single.py igs08
python PythonRDAHMM/rdahmm_eval_single.py fid

~/GPS_Scripts/PythonRDAHMM/unr_splice.py
python ~/GPS_Scripts/PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE
~/GPS_Scripts/PythonRDAHMM/create_mongodb.py UNR_SPLICE


# /home/yuma/PythonRDAHMM/cron_rdahmm_unr.py
# /home/yuma/PythonRDAHMM/unr_splice.py
# /home/yuma/PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_xmls.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_jsons.py UNR_SPLICE
