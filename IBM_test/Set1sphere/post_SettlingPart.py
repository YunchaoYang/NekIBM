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

f1 =  open('SettlingPart1.dat', 'w')
f2 =  open('SettlingPart2.dat', 'w')

f1.write('variables=t,istep,nid,stage,ipt,pid1,pid2,Vx,Vy,Vz,x,y,z,Fx,Fy,Fz \n')
f2.write('variables=t,istep,nid,stage,ipt,pid1,pid2,Vx,Vy,Vz,x,y,z,Fx,Fy,Fz \n')

 # readlastline
pid1 = '0' 
pid2 = '1'

tstep= {} # dictionary of {(time_step,time)}, key:time_step, word: time
latest_file="SettlingPart.dat"
# latest_file="SettlingPart_Force.dat"
try:
   filename=open(latest_file,'r')
   lines=filename.readlines()
   for i in range(len(lines)):
     line=lines[i]
     if line.startswith('Step'):
        line_items = line.split()
        timestep = line_items[1].replace(",","")
        ltime = line_items[3].replace(",","")        
        tstep[timestep] = ltime        
except IOError as e:
   print("Unable to open file") #Does not exist OR no read permissions

try:
   for i in range(len(lines)):
     line=lines[i]
     if line.startswith('Queen '):
        line_items = line.split()
        lstep = line_items[1] # time step
        stage = line_items[3] # stage 
        pid1_ = line_items[4] # pid1 
        pid2_ = line_items[5] # pid2
        
        if lstep in tstep:
           ltime = tstep[lstep]
        else:
           print("Time not found")
        #print(pid1_,pid2_,pid1,pid2)
        if stage == '3':
           if pid1_ == pid1 and pid2_ == pid2:
              f1.write(line.replace("Queen ",ltime))
           else:
              f2.write(line.replace("Queen ",ltime))
        
except IOError as e:
   print("Unable to open file") #Does not exist OR no read permissions
   
