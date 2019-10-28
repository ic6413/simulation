#!/usr/bin/env python
import sys
import read_setting.read_setting as rr
variable_name_string = sys.argv[1]
result = rr.logfile[variable_name_string]
print(result)