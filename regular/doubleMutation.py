import subprocess
import sys
import os
param1= sys.argv[1] 
param2= sys.argv[2] 
param3= sys.argv[3] 
subprocess.call(['./multiMutant.sh', param1, param2, param3])