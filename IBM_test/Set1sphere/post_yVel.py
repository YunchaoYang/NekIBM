#! /usr/bin/python3.6

import re
import os
import subprocess
from subprocess import call
from glob import glob

#get current folder and folder name
Local = os.getcwd()
print(Local + "\n")
folder=os.path.split(os.getcwd())
print (folder[-1])

nxrange     = range(16, 33, 8)
orderNrange = range(4,9,2)
filterrange = range(3,6,1)
ndhrange    = range(15,26,5)
level = 0

outputforce = open('SettlingPart_Force.dat', 'w')
with open('SettlingPart.dat', 'w') as output:
 #output.write("time, yVel, yPos, yForce \n")
 logfiles = glob(os.getcwd()+'/*log*') 
 
 latest_file = max(logfiles, key=os.path.getctime)
 print(latest_file)
 n = 0
 # readlastline
 try:
   filename=open(latest_file,'r')
   lines=filename.readlines()
   for line in lines:
     if line.startswith('Step'):
        line_items = line.split()
        ltime = line_items[3]
        outputforce.write(line)
        output.write(line)
     if line.startswith('Queen '):        
        outputforce.write(line)
     if line.startswith('Queen '):
        output.write(line)

 except IOError as e:
   print("Unable to open file") #Does not exist OR no read permissions
   
outputforce.close()   
