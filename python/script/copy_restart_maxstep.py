#!/usr/bin/env python
import re
import subprocess
import os
import os.path

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

# create folders
for root, dirnames, filenames in os.walk("../output/rst"):
    filenames_resort = natural_sort(filenames)
    filepath = os.path.join(root, filenames_resort[-1])
        
# copy diagram wall    
subprocess.run(["cp", filepath, filenames_resort[-1]])