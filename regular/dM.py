import sys
import os
import subprocess
import requests
import time

##
#Globals and Macros
##
D_DIR = ""
FASTA_SEQ = ""
aminoIndex = {
    0:'a', 1:'r', 2:'n', 3:'d', 4:'c', 5:'q', 6:'e', 7:'g', 8:'h', 9:'i',
    10:'l', 11:'k', 12:'m', 13:'f', 14:'p', 15:'s', 16:'t', 17:'w', 18:'y', 19:'v'}
#aminoIndex = {0:'k', 1:'m', 2:'n', 3:'d'}

PDB_DICT = {}
PDB_SINGLE_DICT = {}
REMOVALS = 0
DEV_NULL = open(os.devnull, 'w')
CHILD_PIDS = []

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

##
##
##
    
def initialize(pdbID, chainID, start, end):
    global D_DIR
    D_DIR = "D_%s.%s.%s_%s" % (pdbID, chainID, start, end)
    createDir(D_DIR) #Creates output file

    global FASTA_SEQ
    response = requests.get("https://www.rcsb.org/fasta/entry/" + pdbID).text.split()
    sequence = response[len(response) - 1][start - 1: end].lower()
    FASTA_SEQ = sequence
    PDB_DICT[FASTA_SEQ] = ""
       
def callProMute(pdbID, chainID, start, end, em = "no ", hphilic = "", hphobic = ""):
    os.chdir('promute')
    print("Gathering PDB files and mutating...\n")
    callProMuteHelper(FASTA_SEQ, pdbID, chainID, start, end, 1, em, hphilic, hphobic)
    os.chdir('..')

def callProMuteHelper(seq, pdbID, chainID, start, end, mutationNumber, em, hphilic, hphobic):
    #global PDB_SINGLE_DICT #Specifying global seems to have no affect...
    #global PDB_DICT
    global REMOVALS
    for residueNum in range(end-start+1):
        for target in aminoIndex:
            targetResidue = aminoIndex.get(target)
            newSeq = changeAt(seq, residueNum, targetResidue).lower()
            if not newSeq in PDB_SINGLE_DICT or mutationNumber == 1:
                if mutationNumber == 1:
                    parameters = "%s %s %d %s %s %s %s" % (pdbID, chainID, residueNum + start + 1, targetResidue, "no", "", "")
                else:
                    parameters = "%s %s %d %s %s %s %s" % (pdbID, chainID, residueNum + start + 1, targetResidue, em, hphilic, hphobic)
                command = "./proMute " + parameters
                newPdbID = ("%s.%s%d%s" % (pdbID, chainID, residueNum + start + 1, targetResidue)).upper()
                if mutationNumber == 1:
                    #print("\n-----SINGLE MUTATION-----")
                    #print(command)
                    #print(newPdbID)
                    PDB_SINGLE_DICT[newSeq] = pdbID
                    #os.system(command)
                    subprocess.call(command, stdout = DEV_NULL, shell = True)
                    callProMuteHelper(newSeq, newPdbID, chainID, start, end, 2, em, hphilic, hphobic)
                    #os.system('rm %s.fasta.txt %s.pdb' %(newPdbID, newPdbID))

                    if(newSeq in PDB_DICT):
                        d = PDB_DICT.pop(newSeq)
                        print("Removing a folder: %s_out" %(d))
                        REMOVALS += 1
                        #os.system('rm -rf %s*' % (d))
                        #os.system('rm -rf ../%s/%s_out' %(D_DIR, d))
                #elif mutationNumber == 2:
                elif not newSeq in PDB_DICT and mutationNumber == 2:
                    PDB_DICT[newSeq] = newPdbID
                    
                    newPID = os.fork()
                    if newPID == 0: #Child Process
                        #print("\n-----SECOND MUTATION-----")
                        #print(command)
                        #print(newPdbID)
                        #subprocess.call(command, shell = True)
                        subprocess.call(command, stdout = DEV_NULL, shell = True)
                        #createDir(newPdbID + '_out')
                        #print("Moving files into an %s_out" % (newPdbID))
                        #os.system('mv %s.fasta.txt %s.pdb %s_out' %(newPdbID, newPdbID, newPdbID))
                        #if em == "em ":
                            #os.system('mv %s_em.pdb %s_out' %(newPdbID, newPdbID))
                        #print("Moving %s_out to the ../%s" % (newPdbID, D_DIR))
                        #os.system('mv %s_out ../%s' %(newPdbID, D_DIR))
                        os._exit(0)
                    else:
                        CHILD_PIDS.append(newPID)

#Delete this
def movePDBs(em):
    os.chdir('promute')
    print("\nOrganizing and moving files over...") 
    for newPdbID in PDB_DICT.values():
        createDir(newPdbID + '_out')
        os.system('mv %s.fasta.txt %s.pdb %s_out' %(newPdbID, newPdbID, newPdbID))
        if(em == "em" or em == "srem"): #Unsure if it captures all possibilities
            try:
                os.system('mv %s_em.pdb %s_out' %(newPdbID, newPdbID))
            except:
                #This should ONLY happen when using srem, since not every file gets energy minimized.
                pass
        os.system('mv %s_out ../%s' %(newPdbID, D_DIR))
    os.chdir('..')

def cleanProMute(pdbID):
    os.chdir('promute')
    os.system('rm -rf %s*' % (pdbID))
    os.chdir('..')

#argv[1] = PDBID
#argv[2] = CHAINID
#argv[3] = RES#_START & RES#_END
#argv[4] = OPTIONAL
def main():
    if len(sys.argv) < 4:
        sys.exit("Please enter valid parameters")

    r = getRange(sys.argv[3])
    initialize(sys.argv[1], sys.argv[2], r[0] + 1, r[1])

    sys.argv[1] = sys.argv[1].upper()
    sys.argv[2] = sys.argv[2].upper()
    emFlag, hphilicFlag, hphobicFlag = "no", "", ""
    for i in range (4, len(sys.argv)):
        flag = sys.argv[i].lower()
        if flag == "em":
            emFlag = flag
        elif flag == "srem":
            emFlag = flag + " "
        elif flag == "hphilic":
            hphilicFlag = flag 
        elif flag == "hphobic":
            hphobicFlag = flag
    
    startTime = time.time()
    callProMute(sys.argv[1], sys.argv[2], r[0], r[1] - 1, emFlag, hphilicFlag, hphobicFlag)
    
    for i, child in enumerate(CHILD_PIDS):
        os.waitpid(child, 0)
    movePDBs(emFlag)
    cleanProMute(sys.argv[1])
    print("\nTime elapsed: %f minutes" % ((time.time() - startTime) / 60))
    print("Folder is %s" % (D_DIR))   
    
if __name__ == "__main__":
    main()
