import subprocess
import sys
import os
import time
param1= sys.argv[1] 
param2= sys.argv[2] 
param3= sys.argv[3] 
subprocess.call(['./multiMutant.sh', param1, param2, param3])

#Dealing with race condition
time.sleep(.5)

cwd = os.getcwd()

#Go into the folder we make
os.chdir(param1+param2+param3+"_out")

#Get any pdb file and transport it back to /regular
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if os.path.splitext(root+name)[-1].lower()==".pdb":
            #os.replace(os.path.dirname(os.path.realpath(__file__)), cwd)
            print(os.path.join(root, name))