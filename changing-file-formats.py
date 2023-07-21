#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 13:49:35 2023

@author: aditya
"""


###############################################################################
#        Reads from a .cssr file and appends molecules from half of the 
#        neighboring cells and writes it to .xyz format. 
###############################################################################


import numpy as np
import os


path = '/Users/aditya/Documents/pormetalomics/FAU/';
os.chdir(path)

# Reads first ths .cssr file to get the box dimensions 

filenameCSSR = "FAU.cssr";

h = open(filenameCSSR, 'r')
content = h.readlines()

print("Total number of lines : ", len(content));

linenumber = 0;
a = 0.0; b = 0.0; c = 0.0; natoms = 0;

for line in content: 
    
    if linenumber == 0:
        
       l = line.split();
       a = float(l[0]);
       b = float(l[1]);
       c = float(l[2]);
    elif linenumber == 2:
       l = line.split();
       natoms = int(l[0]);      
        
        
    linenumber = linenumber + 1;


print("Unit cell lengths are : " , a , b ,c )
print("Total number of atoms: ", natoms);

# Next we read the .xyz file
filenameXYZ = "FAU.xyz";
h = open(filenameXYZ, 'r')
content = h.readlines();

linenumber = 0;
atomName = [];
atomX = np.zeros(natoms);
atomY = np.zeros(natoms);
atomZ = np.zeros(natoms);

for line in content:    
    if linenumber >=2 :
        l = line.split();
        atomName.append(l[0]);
        atomX[linenumber-2] = float(l[1]);
        atomY[linenumber-2] = float(l[2]);
        atomZ[linenumber-2] = float(l[3]);
    linenumber = linenumber + 1;

atomNameNeigh = [];
atomXneigh = [];
atomXaux = [];
atomYneigh = [];
atomYaux = [];
atomZneigh = [];
atomZaux = [];


xt = [-1.0,0.0,1.0];
yt = [-1.0,0.0,1.0];
zt = [-1.0,0.0,1.0];

for ix in xt:
    for iy in yt:
        for iz in zt:
            
            for iatom in range(0,natoms):
                
                catomx = atomX[iatom] + a*ix;
                catomy = atomY[iatom] + b*iy;
                catomz = atomZ[iatom] + c*iz;
                catomName = atomName[iatom]; 
                
                if ( (catomx < 1.5*a) &  (catomx > -0.5*a) & 
                     (catomy < 1.5*b) & (catomy > -0.5*b) & 
                     (catomz < 1.5*c) & (catomz > -0.5*c) ) :
                    
                    atomNameNeigh.append(catomName);
                    atomXneigh.append(catomx);
                    atomYneigh.append(catomy);
                    atomZneigh.append(catomz);
                    
                    
natomsNeigh = len(atomXneigh);
outputFileName = "FAU-neigh.xyz";
mydatafile = open(outputFileName, 'w');

mydatafile.write(str(natomsNeigh) + "\n");
for i in range(0,natomsNeigh):
    
    k = atomNameNeigh[i];
    mydatafile.write(k + " " + str(atomXneigh[i]) + " " + str(atomYneigh[i]) + " " + str(atomZneigh[i]) + "\n")
    

                    