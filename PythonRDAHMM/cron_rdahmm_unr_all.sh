#!/bin/bash

~/GPS_Scripts/PythonRDAHMM/cron_rdahmm_unr.py
~/GPS_Scripts/PythonRDAHMM/unr_splice.py
~/GPS_Scripts/PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE
~/GPS_Scripts/PythonRDAHMM/create_mongodb.py UNR_SPLICE


# /home/yuma/PythonRDAHMM/cron_rdahmm_unr.py
# /home/yuma/PythonRDAHMM/unr_splice.py
# /home/yuma/PythonRDAHMM/rdahmm_eval_single.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_xmls.py UNR_SPLICE
# /home/yuma/PythonRDAHMM/create_summary_jsons.py UNR_SPLICE
