#!/usr/bin/env python
import sys
import read_setting as rr
variable_name_string = sys.argv[1]
result = rr.log_variable[variable_name_string]
print(result)