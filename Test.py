import os
import regex
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors

string='atgcgagctagcggcgatagat'

allSites1=regex.finditer('gc', string, overlapped=True)
allSites2=regex.finditer('cg', string, overlapped=True)

sites1 = [] #List with sites coordinates
for i in allSites1:
    sites1.append(i.start())
print(sites1)

sites2 = [] #List with sites coordinates
for i in allSites2:
    sites2.append(i.start())
print(sites2)

print(sorted(sites1+sites2))

ar=[1,2,3,4,5,6,7,8]
print(ar[2:4])