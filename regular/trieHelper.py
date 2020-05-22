#Trie Helper, this file was created by William Tian on May 22nd, 2020.
#doubleMutation.py is dependent on this, unless the code for redundant sequence checking is all commented out.
#A basic and bare Trie for the purpose of doubleMutation.py.
#Be advised when inserting into the trie. Wrong input will give unexpected results and and throw no error.

class trieNode:
    def __init__(self):
        self.children = [None] * 20
        
#Lower case because most of the FASTA sequence is lower case.
#The order is alphabetical from the amino acid's name, not their abbreviation.
map = {
    'a':0,
    'r':1,
    'n':2,
    'd':3,
    'c':4,
    'q':5,
    'e':6,
    'g':7,
    'h':8,
    'i':9,
    'l':10,
    'k':11,
    'm':12,
    'f':13,
    'p':14,
    's':15,
    't':16,
    'w':17,
    'y':18,
    'v':19
    }

class trieHelper:
    def __init__(self):
        self.root=self.getNode()

    def getNode(self):
        return trieNode()
    
    #Inserts the sequence into the Trie
    #Returns True if the sequence already exists (which means it's a redundant sequence)
    #Returns False if this sequence is unique
    def insertNode(self, seq):
        found = True
        currNode = self.root
        for c in seq:
            i = map.get(c)
            if not currNode.children[i]:
                currNode.children[i] = self.getNode()
                found = False
            currNode = currNode.children[i]
        return found
