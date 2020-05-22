#Double Mutation created by William Tian on May 21st, 2020.
#Uses multiMutant.sh (which in also depends upon ProMute) to mutate a given PDB ID twice.
#The results are consolidated into a file called "D_<PDBID><Chain><Range>_out"
#Redundant results are removed entirely, but the program still goes through the process of creating them.
#This is done through the use of a trie, and examining the FASTA sequences for each PDB
#The trie is created from the file trieHelper.py.
#doubleMutation.py is dependent on trieHelper for deleting redundant sequences.

import sys
import os
from trieHelper import *

#Macros
TEMP_F = "temp_doubleMutationHelper" #Name of the temporary directory where we store the singly mutated sequences.
PDB_T = trieHelper()

def echoPWD():
    print(os.popen('echo $PWD').read())

#Returns the path for the folder outputed by multiMutant.sh. It'll look like "<ID><CHAIN><RANGE>_out"
def getPath(argv):
    return argv[1] + argv[2] + argv[3] + '_out'

def createDir(t_f):
    if(os.path.exists(t_f)):
        os.system('rm -rf ' + t_f)
    os.system('mkdir ' + t_f)    

def getFASTA(fileName):
    with open(fileName) as f:
        return f.read()

#Deletes the leftover files from multiMutant. I thought they were annoying.
def cleanMultiMutant(argv):
    os.system('rm -rf ' + getPath(argv))

    #Commenting this out so you don't have to redownload the PDB file every time you want to call on the same sequence
    #os.system('rm promute/' + argv[1] + '.pdb')
    
#Mutates a given sequence with the given parameters.
#argv[0] is never used. argv[1] is the PDB ID.
#argv[2] is the Chain (Note: Case sensitive). argv[3] is the range (Note: Inclusive on both ends).
#Essentially just calls ./multiMutant.sh with the command line arguments passed
def callMultiMutant(argv):
    fPath = getPath(argv)
    os.system('rm -rf ' + fPath)

    print('Calling ./multiMutant.sh ' + argv[1] + ' ' + argv[2] + ' ' + argv[3])
        
    os.popen('./multiMutant.sh ' + argv[1] + ' ' + argv[2] + ' ' + argv[3]).read()  
    os.system('find ./' + fPath + ' -type d -empty -delete') #Removes redundant folders   

#MultiMutant.sh creates a bunch of PDB files across multiple folders in a new directory.
#This gathers all those files and puts them into a single temporary folder. Specified by t_f.
#This new folder will be discarded once the PDB files are mutated again, giving the original sequence two mutations.
def gatherPDBs(argv, t_f):
    fPath = getPath(argv)
    
    createDir(t_f)
    
    os.chdir('./' + fPath)
    for dirs in os.walk('.', topdown = False):
        if(dirs[0] != '.'):
            os.system('mv ' + dirs[0] + '/' + dirs[0] + '.pdb' + ' ../' + t_f)
    os.chdir('..')

#Moves the entire folder of sequences into one main folder.
def gatherDoubles(argv, t_f):
    fPath = getPath(argv)
    os.system('mv ' + fPath + '/* ' + t_f)
    os.system('rm -rf ' + fPath + ' ' + t_f + '/sequentialPipelineInvocation.sh')

#Uses PDB_T, a trie, to check if a FASTA sequence is found. If it is, we remove the folder, as it's redundant.
def removeRedundants(workingDir):
    print("Removing redundant sequences...")
    os.chdir(workingDir)

    i = 0
    j = 0
    for dirs in os.walk('.', topdown = False):
        if(dirs[0] != '.'):
            seq = getFASTA(dirs[0] + '/' + dirs[0] + '.fasta.txt').lower()
            if trieHelper.insertNode(PDB_T, seq) == True: #If True then we remove, since it's a redundant sequence
                os.system('rm -rf ' + dirs[0])
                i += 1
            j += 1
    os.chdir('..') #Changes back to the root directory for this file
    os.system('rm -rf ' + TEMP_F) #Removes the temporary folder which we used for our singly mutated sequences
    print("Removed %d redundant sequences out of %d total sequences." % (i, j))
    print("There are now %d sequences total." % (j - i))
    
# Grabs the temporary directory, and mutates everything in it again.
#Those results are then moved into the final directory, under the variable "dir"
def mutateDirectory(argv):
    dir = 'D_' + argv[1] +  argv[2] + argv[3] + '_out'
    createDir(dir)
    
    for file in os.listdir(TEMP_F):
        if file.endswith(".pdb"): #Add "and file.startswith()" to troubleshoot.
            os.system('mv ' + TEMP_F + '/' + file + ' promute') #Moves the current file in the loop into ./promute/
            file = file.replace('.pdb', '')
            temp_argv = [0, file, argv[2], argv[3]] #Create a new argv, with the PDB ID as the file we just obtained
            callMultiMutant(temp_argv)
            os.system('rm ./promute/' + file + '.pdb')

            gatherDoubles(temp_argv, dir)
    
def main():
    if(len(sys.argv) < 4):
        sys.exit("Please enter the correct command line arguments")
    callMultiMutant(sys.argv)
    gatherPDBs(sys.argv, TEMP_F)
    mutateDirectory(sys.argv)
    removeRedundants('D_' + getPath(sys.argv))
    cleanMultiMutant(sys.argv)
    
if __name__ == "__main__":
    main()
