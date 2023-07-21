#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 13:08:53 2023

@author: aditya

# From the LAMMPS dumped file, this script generates a grid that has 
a) the distance grids
b) the stress grids that has been interpolated from the particle to the 
grid using the vtkShepherdMethod
c) This can be read into Paraview for segmentation analysis. 

"""

import vtk
import numpy as np
import os

# Let us first read the dump.stress. file from LAMMPS:
    
path  = "/Users/aditya/Documents/pormetalomics/testing-poreTDA/gold-structure";
os.chdir(path);

filename = "dump.stress.75000";

h = open(filename, 'r')

content = h.readlines()
    
linenumber = 0;
    
for line in content:
    
    if linenumber == 1:
        timestamp = int(line)
    elif linenumber == 3:
        natoms = int(line)
    elif linenumber == 5:
        l = line.split(' ');
        astart = float(l[0]);
        aend = float(l[1]);
    elif linenumber == 6:
        l = line.split(' ');
        bstart = float(l[0]);
        bend = float(l[1]);
    elif linenumber == 7:
        l = line.split(' ');
        cstart = float(l[0]);
        cend = float(l[1]);
    linenumber += 1
        
print(timestamp, natoms, astart, aend, bstart, bend, cstart, cend);

atomID = np.zeros(natoms, dtype = int);
atomType = np.zeros(natoms, dtype = int);
atomX = np.zeros(natoms);
atomY = np.zeros(natoms);
atomZ = np.zeros(natoms);
potentialEnergy = np.zeros(natoms);
stressXX = np.zeros(natoms);
stressYY = np.zeros(natoms);
stressZZ = np.zeros(natoms);
stressXY = np.zeros(natoms);
stressXZ = np.zeros(natoms);
stressYZ = np.zeros(natoms); 

linenumber = 0
    
for line in content:
        
    if linenumber > 8:
        
        lsplit = line.split(' ');
        atomID[linenumber-9] = int (lsplit[0])
        atomType[linenumber-9] = int (lsplit[1])
        atomX[linenumber-9] = float (lsplit[2])
        atomY[linenumber-9] = float (lsplit[3])
        atomZ[linenumber-9] = float (lsplit[4])
        potentialEnergy[linenumber-9] = float (lsplit[5]);
        stressXX[linenumber-9] = float (lsplit[5]);
        stressYY[linenumber-9] = float (lsplit[6]);
        stressZZ[linenumber-9] = float (lsplit[7]);
        stressXY[linenumber-9] = float (lsplit[8]);
        stressXZ[linenumber-9] = float (lsplit[9]);
        stressYZ[linenumber-9] = float (lsplit[10]);
        
    linenumber += 1
    

# We will also read the distance grid that is actually generated from zeoplusplus. 
# Load the host cube file
filenameDistanceGrid = "15.cube";

h = open(filenameDistanceGrid, 'r')
content = h.readlines()

print("Total number of lines = ", len(content))
linenumber = 0;

for line in content:
    
    if linenumber == 2:
        l = line.split('     ');
        natoms = int(l[0]);
    elif linenumber == 3:
        l = line.split('     ');
        na = int(l[0]);
        gridresA = float(l[1]);
    elif linenumber == 4:
        l = line.split('     ');
        nb = int(l[0]);
        gridresB = float(l[2]);
    elif linenumber == 5:
        l = line.split('     ');
        nc = int(l[0]);
        gridresC = float(l[3]);

    linenumber += 1
    
print("number of atoms = ", natoms)
print("Number of grid points: (na,nb,nc) = (", na, ",", nb, ",", nc,")")
print("Grid Resolution = (", gridresA, ",", gridresB, ",", gridresC, ")")


# Read the distance function 
cube = np.zeros((na,nb,nc));

linenumber = 0;
count2 = 0;
ix=0;iy=0;iz=0;

for line in content:
    
    if((linenumber > 6 + natoms)):
        lsplit = line.split('  ');
        for i in range(0,len(lsplit)):
            if (len(lsplit[i]) > 0) :
                cube[ix][iy][iz] = float(lsplit[i]);
                iz+=1;
        if (iz == nc):
            iz = 0;
            iy+=1
            if(iy == nb ):
                iy = 0;
                ix+=1;
    linenumber+=1
    

#%%  

points = vtk.vtkPoints();

for i in range(0,natoms) : 
    
    points.InsertPoint(i,atomX[i], atomY[i],atomZ[i]);
    
scalars = vtk.vtkFloatArray()
scalars.SetName("potentialEnergy")

for i in range(0,natoms) : 
    
    scalars.InsertValue(i, potentialEnergy[i]);

profile = vtk.vtkPolyData()
profile.SetPoints(points)
profile.GetPointData().SetScalars(scalars)

dim = 21

shepard1 = vtk.vtkShepardMethod()
shepard1.SetInputData(profile)
shepard1.SetModelBounds(astart,aend,bstart,bend, cstart,cend)
shepard1.SetSampleDimensions(dim,dim,dim)
shepard1.SetNullValue(0)
shepard1.SetMaximumDistance(1)
shepard1.SetPowerParameter(2)

timer = vtk.vtkExecutionTimer()
timer.SetFilter(shepard1)
shepard1.Update()
wallClock = timer.GetElapsedWallClockTime()
print ("Shephard (P=2):", wallClock)

xmlWriter = vtk.vtkXMLImageDataWriter()
outputFilename = "interpolate.vti";
xmlWriter.SetFileName(outputFilename);
xmlWriter.SetInputConnection(shepard1.GetOutputPort());
xmlWriter.Write();

