#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 13:08:53 2023

@author: aditya

# From the LAMMPS dumped file, this script generates a grid that has 
a) the distance grids
b) the stress grids that has been interpolated from the particle to the 
grid using the vtkShepherdMethod
c) The grid is a volume representation, and hence can be used in all the 
filters of paraview such as contouring, 

"""
import vtk
import numpy as np
import os
import sys

##############################################################################
# Parameters for the shepherd method:

# A) Max distance is the fraction of the diagonal of the box that is used to 
# interpolate. If maxDistance = 1.0, then all the points are used to interpolate
# Default value is 0.25. 

# B) PowerParameter is the inverse square weights for interpolation, which is ideal. 
    
maxDistance = 0.1;
powerParameter = 2;

##############################################################################

if (len(sys.argv) != 4):

    print(sys.argv)
    print("Enter three arguments, the dump.stress file, .cube file and a counter int (give 0 if only one file)");
    exit();

# Let us first read the dump.stress. file from LAMMPS:
path  = "/home/aditya/Documents/israel-structures/generation-of-distance-functions/stress-grids";
os.chdir(path);

filename = sys.argv[1];

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
     
print("Data read from", filename, " file");
print("File id: ", timestamp)
print("Number of atoms", natoms)
print("Cell length X : ", astart, aend)
print("Cell length Y : ", bstart, bend)
print("Cell length Z : ", cstart, cend)

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
        stressXX[linenumber-9] = float (lsplit[6]);
        stressYY[linenumber-9] = float (lsplit[7]);
        stressZZ[linenumber-9] = float (lsplit[8]);
        stressXY[linenumber-9] = float (lsplit[9]);
        stressXZ[linenumber-9] = float (lsplit[10]);
        stressYZ[linenumber-9] = float (lsplit[11]);
        
    linenumber += 1
    

# We will also read the distance grid that is actually generated from zeoplusplus. 
# Load the  cube file
filenameDistanceGrid = sys.argv[2];

h = open(filenameDistanceGrid, 'r')
content = h.readlines()

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

print("\n")    
print("Data read from the distance grid"); 
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
    

# Next, we create a polyData structure that can go as input to the 
# shepherd's Method. 

points = vtk.vtkPoints();

for i in range(0,natoms) : 
    
    points.InsertPoint(i,atomX[i], atomY[i],atomZ[i]);
    
PEarray = vtk.vtkFloatArray()
PEarray.SetName("potentialEnergy")
sigmaXXarray = vtk.vtkFloatArray();
sigmaXXarray.SetName("sigmaXX");
sigmaYYarray = vtk.vtkFloatArray();
sigmaYYarray.SetName("sigmaYY");
sigmaZZarray = vtk.vtkFloatArray();
sigmaZZarray.SetName("sigmaZZ");
sigmaXYarray = vtk.vtkFloatArray();
sigmaXYarray.SetName("sigmaXY");
sigmaXZarray = vtk.vtkFloatArray();
sigmaXZarray.SetName("sigmaXZ");
sigmaYZarray = vtk.vtkFloatArray();
sigmaYZarray.SetName("sigmaYZ");


for i in range(0,natoms) : 
    
    PEarray.InsertValue(i, potentialEnergy[i]);
    sigmaXXarray.InsertValue(i, stressXX[i]);
    sigmaYYarray.InsertValue(i, stressYY[i]);
    sigmaZZarray.InsertValue(i, stressZZ[i]);
    sigmaXYarray.InsertValue(i, stressXY[i]);
    sigmaXZarray.InsertValue(i, stressXZ[i]);
    sigmaYZarray.InsertValue(i, stressYZ[i]);



dimX = na; dimY = nb; dimZ = nc;

# PE interpolation
PEPoly = vtk.vtkPolyData()
PEPoly.SetPoints(points)
PEPoly.GetPointData().SetScalars(PEarray)

shepardPE = vtk.vtkShepardMethod()
shepardPE.SetInputData(PEPoly)
shepardPE.SetModelBounds(astart,aend,bstart,bend, cstart,cend)
shepardPE.SetSampleDimensions(dimX,dimY,dimZ)
shepardPE.SetNullValue(0)
shepardPE.SetMaximumDistance(maxDistance)
shepardPE.SetPowerParameter(powerParameter)

timer = vtk.vtkExecutionTimer()
timer.SetFilter(shepardPE)
shepardPE.Update()
wallClock = timer.GetElapsedWallClockTime()
print ("Shephard PE completed in:", wallClock)

# Sigma XX interpolation

sigmaXXPoly = vtk.vtkPolyData();
sigmaXXPoly.SetPoints(points);
sigmaXXPoly.GetPointData().SetScalars(sigmaXXarray);



shepardSigmaXX = vtk.vtkShepardMethod()
shepardSigmaXX.SetInputData(sigmaXXPoly)
shepardSigmaXX.SetModelBounds(astart,aend,bstart,bend, cstart,cend)
shepardSigmaXX.SetSampleDimensions(dimX,dimY,dimZ)
shepardSigmaXX.SetNullValue(0)
shepardSigmaXX.SetMaximumDistance(maxDistance)
shepardSigmaXX.SetPowerParameter(powerParameter)

timer = vtk.vtkExecutionTimer()
timer.SetFilter(shepardSigmaXX)
shepardSigmaXX.Update()
wallClock = timer.GetElapsedWallClockTime()
print ("Shephard sigmaXX completed in :", wallClock)

# Sigma YY interpolation
sigmaYYPoly = vtk.vtkPolyData();
sigmaYYPoly.SetPoints(points);
sigmaYYPoly.GetPointData().SetScalars(sigmaYYarray);


shepardSigmaYY = vtk.vtkShepardMethod()
shepardSigmaYY.SetInputData(sigmaYYPoly)
shepardSigmaYY.SetModelBounds(astart,aend,bstart,bend, cstart,cend)
shepardSigmaYY.SetSampleDimensions(dimX,dimY,dimZ)
shepardSigmaYY.SetNullValue(0)
shepardSigmaYY.SetMaximumDistance(maxDistance)
shepardSigmaYY.SetPowerParameter(powerParameter)

timer = vtk.vtkExecutionTimer()
timer.SetFilter(shepardSigmaYY)
shepardSigmaYY.Update()
wallClock = timer.GetElapsedWallClockTime()
print ("Shephard sigmaYY completed in :", wallClock)

# Sigma ZZ interpolation

sigmaZZPoly = vtk.vtkPolyData();
sigmaZZPoly.SetPoints(points);
sigmaZZPoly.GetPointData().SetScalars(sigmaZZarray);


shepardSigmaZZ = vtk.vtkShepardMethod()
shepardSigmaZZ.SetInputData(sigmaZZPoly)
shepardSigmaZZ.SetModelBounds(astart,aend,bstart,bend, cstart,cend)
shepardSigmaZZ.SetSampleDimensions(dimX,dimY,dimZ)
shepardSigmaZZ.SetNullValue(0)
shepardSigmaZZ.SetMaximumDistance(maxDistance)
shepardSigmaZZ.SetPowerParameter(powerParameter)

timer = vtk.vtkExecutionTimer()
timer.SetFilter(shepardSigmaZZ)
shepardSigmaZZ.Update()
wallClock = timer.GetElapsedWallClockTime()
print ("Shephard sigmaZZ completed in:", wallClock)

# Sigma XY interpolation

sigmaXYPoly = vtk.vtkPolyData();
sigmaXYPoly.SetPoints(points);
sigmaXYPoly.GetPointData().SetScalars(sigmaXYarray);

shepardSigmaXY = vtk.vtkShepardMethod()
shepardSigmaXY.SetInputData(sigmaXYPoly)
shepardSigmaXY.SetModelBounds(astart,aend,bstart,bend, cstart,cend)
shepardSigmaXY.SetSampleDimensions(dimX,dimY,dimZ)
shepardSigmaXY.SetNullValue(0)
shepardSigmaXY.SetMaximumDistance(maxDistance)
shepardSigmaXY.SetPowerParameter(powerParameter)

timer = vtk.vtkExecutionTimer()
timer.SetFilter(shepardSigmaXY)
shepardSigmaXY.Update()
wallClock = timer.GetElapsedWallClockTime()
print ("Shephard sigmaXY completed in :", wallClock)


# Sigma XZ interpolation

sigmaXZPoly = vtk.vtkPolyData();
sigmaXZPoly.SetPoints(points);
sigmaXZPoly.GetPointData().SetScalars(sigmaXZarray);

shepardSigmaXZ = vtk.vtkShepardMethod()
shepardSigmaXZ.SetInputData(sigmaXZPoly)
shepardSigmaXZ.SetModelBounds(astart,aend,bstart,bend, cstart,cend)
shepardSigmaXZ.SetSampleDimensions(dimX,dimY,dimZ)
shepardSigmaXZ.SetNullValue(0)
shepardSigmaXZ.SetMaximumDistance(maxDistance)
shepardSigmaXZ.SetPowerParameter(powerParameter)

timer = vtk.vtkExecutionTimer()
timer.SetFilter(shepardSigmaXZ)
shepardSigmaXZ.Update()
wallClock = timer.GetElapsedWallClockTime()
print ("Shephard sigmaXZ completed in :", wallClock)


# Sigma YZ interpolation

sigmaYZPoly = vtk.vtkPolyData();
sigmaYZPoly.SetPoints(points);
sigmaYZPoly.GetPointData().SetScalars(sigmaYZarray);

shepardSigmaYZ = vtk.vtkShepardMethod()
shepardSigmaYZ.SetInputData(sigmaYZPoly)
shepardSigmaYZ.SetModelBounds(astart,aend,bstart,bend, cstart,cend)
shepardSigmaYZ.SetSampleDimensions(dimX,dimY,dimZ)
shepardSigmaYZ.SetNullValue(0)
shepardSigmaYZ.SetMaximumDistance(maxDistance)
shepardSigmaYZ.SetPowerParameter(powerParameter)

timer = vtk.vtkExecutionTimer()
timer.SetFilter(shepardSigmaYZ)
shepardSigmaYZ.Update()
wallClock = timer.GetElapsedWallClockTime()
print ("Shephard sigmaYZ completed in :", wallClock)

# We create next the image grid to store the interpolated values.
grid = vtk.vtkImageData();
grid.SetExtent(0,dimX-1, 0, dimY-1, 0, dimZ-1);
grid.SetSpacing(gridresA,gridresB,gridresC);


distanceGrid = vtk.vtkDoubleArray();
distanceGrid.SetName("This is distance grid");
distanceGrid.SetNumberOfValues(grid.GetNumberOfPoints());

PEgrid = vtk.vtkDoubleArray();
PEgrid.SetName("PotentialEnergy");
PEgrid.SetNumberOfValues(grid.GetNumberOfPoints());

sigmaXXGrid = vtk.vtkDoubleArray();
sigmaXXGrid.SetName("sigmaXX");
sigmaXXGrid.SetNumberOfValues(grid.GetNumberOfPoints());

sigmaYYGrid = vtk.vtkDoubleArray();
sigmaYYGrid.SetName("sigmaYY");
sigmaYYGrid.SetNumberOfValues(grid.GetNumberOfPoints());

sigmaZZGrid = vtk.vtkDoubleArray();
sigmaZZGrid.SetName("sigmaZZ");
sigmaZZGrid.SetNumberOfValues(grid.GetNumberOfPoints());

sigmaXYGrid = vtk.vtkDoubleArray();
sigmaXYGrid.SetName("sigmaXY");
sigmaXYGrid.SetNumberOfValues(grid.GetNumberOfPoints());

sigmaXZGrid = vtk.vtkDoubleArray();
sigmaXZGrid.SetName("sigmaXZ");
sigmaXZGrid.SetNumberOfValues(grid.GetNumberOfPoints());

sigmaYZGrid = vtk.vtkDoubleArray();
sigmaYZGrid.SetName("sigmaYZ");
sigmaYZGrid.SetNumberOfValues(grid.GetNumberOfPoints());

for i in range(0,grid.GetNumberOfPoints()):

    distanceGrid.SetValue(i,0.0);
    PEgrid.SetValue(i,0.0);
    sigmaXXGrid.SetValue(i,0.0);
    sigmaYYGrid.SetValue(i,0.0);
    sigmaZZGrid.SetValue(i,0.0);
    sigmaXYGrid.SetValue(i,0.0);
    sigmaXZGrid.SetValue(i,0.0);
    sigmaYZGrid.SetValue(i,0.0);

grid.GetPointData().AddArray(distanceGrid);
grid.GetPointData().AddArray(PEgrid);
grid.GetPointData().AddArray(sigmaXXGrid);
grid.GetPointData().AddArray(sigmaYYGrid);
grid.GetPointData().AddArray(sigmaZZGrid);
grid.GetPointData().AddArray(sigmaXYGrid);
grid.GetPointData().AddArray(sigmaXZGrid);
grid.GetPointData().AddArray(sigmaYZGrid);


imageDataPE = shepardPE.GetOutput();
imageDataSXX = shepardSigmaXX.GetOutput();
imageDataSYY = shepardSigmaYY.GetOutput();
imageDataSZZ = shepardSigmaZZ.GetOutput();
imageDataSXY = shepardSigmaXY.GetOutput();
imageDataSXZ = shepardSigmaXZ.GetOutput();
imageDataSYZ = shepardSigmaYZ.GetOutput();


for i in range(0,dimX):

    for j in range(0,dimY):
        
        for k in range(0, dimZ):
            
            index = imageDataPE.GetScalarIndex(i,j,k);
            grid.GetPointData().GetArray("This is distance grid").SetVariantValue(index,cube[i][j][k]);
            grid.GetPointData().GetArray("PotentialEnergy").SetVariantValue(index,imageDataPE.GetScalarComponentAsDouble(i,j,k,0))
            grid.GetPointData().GetArray("sigmaXX").SetVariantValue(index,imageDataSXX.GetScalarComponentAsDouble(i,j,k,0))
            grid.GetPointData().GetArray("sigmaYY").SetVariantValue(index,imageDataSYY.GetScalarComponentAsDouble(i,j,k,0))
            grid.GetPointData().GetArray("sigmaZZ").SetVariantValue(index,imageDataSZZ.GetScalarComponentAsDouble(i,j,k,0))
            grid.GetPointData().GetArray("sigmaXY").SetVariantValue(index,imageDataSXY.GetScalarComponentAsDouble(i,j,k,0))
            grid.GetPointData().GetArray("sigmaXZ").SetVariantValue(index,imageDataSXZ.GetScalarComponentAsDouble(i,j,k,0))
            grid.GetPointData().GetArray("sigmaYZ").SetVariantValue(index,imageDataSYZ.GetScalarComponentAsDouble(i,j,k,0))


gridWriter = vtk.vtkXMLImageDataWriter();
filecounter = sys.argv[3];
outputFileName = "grid_" + filecounter + ".vti";
gridWriter.SetFileName(outputFileName);
gridWriter.SetInputData(grid);
gridWriter.Write();
