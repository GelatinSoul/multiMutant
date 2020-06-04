import sys
import os
from os import path

##
#Helper functions
##
def getRange(range):
    r = range.split(':')
    #Many Python commands aren't inclusive at the end, but multiMutant is which is what this code is modeled after.
    return int(r[0]) - 1, int(r[1])

def createDir(d):
    if os.path.exists(d):
       os.system('rm -rf ' + d)
    os.system('mkdir ' + d)

def changeAt(string, index, target):
    return string[:index] + target + string[index + 1:]
    
    
#needs the following input pdbID, chainID, range of aminos, and wether or not em was used
#if em was used, include the -em flag at the end
def main():   
    if len(sys.argv) < 4:
        sys.exit("Please enter valid parameters")
        
    pdbID = sys.argv[1]
    chainID = sys.argv[2]
    r = getRange(sys.argv[3])
    filename = "D_"+pdbID+chainID+str(r[0] + 1)+":"+str(r[1])+"_out"
    if len(sys.argv) == 5:
        em = sys.argv[4]
    else:
        em = "-1"
    print("Running on file: " + filename)
    os.chdir(filename)
    
    allFilesExist = True

    #loops through all the files
    for root, dirs, files in os.walk(".", topdown=True):
        for name in dirs:
            os.chdir(name)
            #loops through to check files exist and arent empty
            for root, dirs, files in os.walk(".", topdown=True):
                numberCorrect = 0
                for name in files:
                    if(name.endswith('.fasta.txt') or name.endswith('.pdb')):
                        if(os.stat(name).st_size ==0):
                            print("File: "+ name+ " exists but is empty")
                        else:
                            numberCorrect = numberCorrect + 1
                    
                if(em != "-1"):
                    if(numberCorrect != 3):
                        allFilesExist = False
                else:
                    if(numberCorrect != 2):
                        allFilesExist = False
            os.chdir("..")
    print("==================")
    if(allFilesExist == False):
        print("File is missing") 
    else:
        print("All files were found and are not empty")


    
    
if __name__ == "__main__":
    main()
