import sys
import os
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
       
def callProMute(pdbID, chainID, start, end, opt):
    os.chdir('promute')
    print("Gathering PDB files and mutating...\n")
    callProMuteHelper(FASTA_SEQ, pdbID, chainID, start, end, opt, 1)
    os.chdir('..')

def callProMuteHelper(seq, pdbID, chainID, start, end, opt, mutationNumber):
    for residueNum in range(end-start+1):
        for target in aminoIndex:
            targetResidue = aminoIndex.get(target)
            newSeq = changeAt(seq, residueNum, targetResidue).lower()
            if not newSeq in PDB_SINGLE_DICT or mutationNumber == 1:
                parameters = "%s %s %d %s %s" % (pdbID, chainID, residueNum + start + 1, targetResidue, opt)
                command = "./proMute " + parameters
                newPdbID = ("%s.%s%d%s" % (pdbID, chainID, residueNum + start + 1, targetResidue)).upper()
                if mutationNumber == 1:
                    PDB_SINGLE_DICT[newSeq] = 1
                    os.system(command)
                    callProMuteHelper(newSeq, newPdbID, chainID, start, end, opt, 2)
                    os.system('rm %s.fasta.txt %s.pdb' %(newPdbID, newPdbID))

                    if(newSeq in PDB_DICT):
                        d = PDB_DICT.pop(newSeq)
                        os.system('rm -rf ../%s/%s_out' %(D_DIR, d))
                elif mutationNumber == 2:
                    PDB_DICT[newSeq] = newPdbID
                    os.system(command)
                    createDir(newPdbID + '_out')
                    os.system('mv %s.fasta.txt %s.pdb %s_out' %(pdbID, pdbID, pdbID))
                    os.system('mv %s_out ../%s' %(pdbID, D_DIR))
                    
def movePDBs():
    os.chdir('promute')
    print("\nOrganizing and moving files over...")
    for pdbID in PDB_DICT.values():
        createDir(pdbID + '_out')
        os.system('mv %s.fasta.txt %s.pdb %s_out' %(pdbID, pdbID, pdbID))
        os.system('mv %s_out ../%s' %(pdbID, D_DIR))
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

    startTime = time.time()
    if len(sys.argv) == 5:
        callProMute(sys.argv[1], sys.argv[2], r[0], r[1]-1, sys.argv[4])
    else:
        callProMute(sys.argv[1], sys.argv[2], r[0], r[1]-1, 'no')
        
    #movePDBs()

    print("\nTime elapsed: %f minutes" % ((time.time() - startTime) / 60))
    print("Folder is %s" % (D_DIR)) 
        
if __name__ == "__main__":
    main()
