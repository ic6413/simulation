#!/usr/bin/env python
import sys
import read_setting as rr

id_i = int(sys.argv[1])
step1 = sys.argv[2]
step2 = sys.argv[3]

# check if atom i in trace id
if int(rr.logfile['freq_dump_single_trace']) == 1 and id_i in [int(rr.logfile['trace_id1']), int(rr.logfile['trace_id2']), int(rr.logfile['trace_id3'])]:
    use_trace_or_all = "trace"
# check if freq_dump_single_all is 1
elif int(rr.logfile['freq_dump_single_all']) == 1:
    use_trace_or_all = "all"

print(use_trace_or_all)