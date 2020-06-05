import sys
import os
import subprocess
import threading
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

PDB_DICT = {}
PDB_SINGLE_DICT = {}
DEV_NULL = open(os.devnull, 'w')

THREADS = []
PROMUTES = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0} #We have 8 Promute folders. May change to 16 later...
MAXTHREADS = len(PROMUTES) + 1 #The 1 is because one of the threads is the original
SEMA = threading.BoundedSemaphore(len(PROMUTES))
MUTEX = threading.Lock()

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
       
def callProMute(pdbID, chainID, start, end, em = "no", hphilic = "", hphobic = ""):
    os.chdir('promute')
    print("Gathering PDB files and mutating...\n")
    callProMuteHelper(FASTA_SEQ, pdbID, chainID, start, end, 1, em, hphilic, hphobic)
    os.chdir('..')

#Recursive function that is very shallow but wide.
#The deepest it ever goes is 2 calls, but it makes 2 calls many times. Uses threads during the second call
def callProMuteHelper(seq, pdbID, chainID, start, end, mutationNumber, em, hphilic, hphobic):
    for residueNum in range(end-start+1):
        for target in aminoIndex:
            targetResidue = aminoIndex.get(target)
            newSeq = changeAt(seq, residueNum, targetResidue).lower()
            if not newSeq in PDB_SINGLE_DICT or mutationNumber == 1:
                if mutationNumber == 1:
                    parameters = "%s %s %d %s %s %s %s" % (pdbID, chainID, residueNum + start + 1, targetResidue, "no", "", "")
                else:
                    parameters = "%s %s %d %s %s %s %s" % (pdbID, chainID, residueNum + start + 1, targetResidue, "em", hphilic, hphobic)
                command = "./proMute " + parameters
                newPdbID = ("%s.%s%d%s" % (pdbID, chainID, residueNum + start + 1, targetResidue)).upper()
                if mutationNumber == 1:
                    PDB_SINGLE_DICT[newSeq] = ""
                    subprocess.call(command, stdout = DEV_NULL, shell = True)
                    callProMuteHelper(newSeq, newPdbID, chainID, start, end, 2, em, hphilic, hphobic)

                    if(newSeq in PDB_DICT):
                        d = PDB_DICT.pop(newSeq)
                #elif mutationNumber == 2:
                elif not newSeq in PDB_DICT and mutationNumber == 2:
                    PDB_DICT[newSeq] = newPdbID

                    command_wrapper = [command, pdbID, newPdbID]
                    t = threading.Thread(target=proMuteThreadWrapper, args = command_wrapper)
                    THREADS.append(t)
                    while True:
                        time.sleep(1)
                        if threading.active_count() < MAXTHREADS:
                           t.start() #We create a thread.
                           break

#Mutexes to place themselves in a dictionary.
#Then uses the "index" of that dictionary to call ProMute in one of the copied ProMute folders.
#Also uses a Semaphore because there's only so many ProMute copies. The Semaphore may be redundant.
def proMuteThreadWrapper(command, pdbID, newPdbID):
    threadData = threading.local()
    threadData.i = -1
    SEMA.acquire()
    while threadData.i == -1:
        MUTEX.acquire()
        for x, y in PROMUTES.items():
            if y == 0:
                PROMUTES[x] = 1
                threadData.i = x
                break
        MUTEX.release()

    os.system('cp ./%s.pdb ./promute_%s' % (pdbID, threadData.i))
    command = './promute_' + str(threadData.i) + '/' + command
    print(command)
    subprocess.call(command, stdout = DEV_NULL, shell = True)    

    MUTEX.acquire()
    PROMUTES[threadData.i] = 0
    MUTEX.release()
    SEMA.release()
    
def movePDBs(em):
    os.chdir('promute')
    print("\nOrganizing and moving files over...") 
    for newPdbID in PDB_DICT.values():
        createDir(newPdbID + '_out')
        os.system('mv %s.fasta.txt %s.pdb %s_em.pdb %s_out' %(newPdbID, newPdbID, newPdbID, newPdbID))
        
        ## Commented out this block of code to ensure that it always copies over the em.pdb
        
        #if(em == "em" or em == "srem"): #Unsure if it captures all possibilities. This may break if someone uses weird options or commands
            #try:
        #os.system('mv %s_em.pdb %s_out' %(newPdbID, newPdbID))
           # except:
                #This should ONLY happen when using srem, since not every file gets energy minimized.
            #    pass
        os.system('mv %s_out ../%s' %(newPdbID, D_DIR))
    os.chdir('..')

def cleanProMute(pdbID):
    os.chdir('promute')
    print("\nDeleting leftover files...")
    os.system('rm -rf %s*' % (pdbID))
    for i in range(MAXTHREADS - 1):      
        os.system('rm -rf /promute_%s/%s*' % (i, pdbID))
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
        if flag == "em" or flag == "-em":
            emFlag = flag
        elif flag == "srem" or flag == "-srem":
            emFlag = flag
        elif flag == "hphilic" or flag == "-hphilic":
            hphilicFlag = flag 
        elif flag == "hphobic" or flag == "-hphobic":
            hphobicFlag = flag
    
    startTime = time.time()
    callProMute(sys.argv[1], sys.argv[2], r[0], r[1] - 1, emFlag, hphilicFlag, hphobicFlag)
    
    for t in THREADS:
        t.join()
    
    movePDBs(emFlag)
    cleanProMute(sys.argv[1])
    print("\nTime elapsed: %f minutes" % ((time.time() - startTime) / 60))
    print("Folder is %s" % (D_DIR))   
    
if __name__ == "__main__":
    main()
