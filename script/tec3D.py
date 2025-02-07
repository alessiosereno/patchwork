#!/bin/python3

from pathlib import Path
import math
import numpy as np
import re, os
import sys
from pyiner import fini_get
from mesh_refine import refine_by_splitting

print(f"  Tecplot 3D mesh generation")
print()

ini_file = 'input.ini'
if not Path(ini_file).is_file():
    print('  ERROR : INI file not found -> mesh.tec generation interrupted')
    exit()

circular_sector = False
threedim = False

outpath = fini_get(ini_file, 'PATCHWORK', 'outpath')
if outpath:
    fileoutput = outpath + "/mesh.tec"

fs = fini_get(ini_file, 'PATCHWORK', 'fs', float)
if not fs:
    fs = 1.0

includefile = fini_get(ini_file, 'PATCHWORK', 'includefile')
if includefile:
    ini_file = includefile

fileinput = fini_get(ini_file, 'PATCHWORK', 'inputfile')
if fileinput:
    expanded_path = os.path.expanduser(fileinput)
    folder_path = os.path.dirname(expanded_path)
    if len(folder_path) == 0:
        folder_path = "."
    method = fini_get(ini_file, 'PATCHWORK', 'method')
    if method:
        fileinput = folder_path + "/mesh.dat"
        fileoutput = folder_path + "/mesh.tec"
    else:
        fileoutput = folder_path + "/mesh.tec"
elif outpath:
    fileinput = outpath + "/mesh.dat"
    fileoutput = outpath + "/mesh.tec"
else:
    fileinput = 'mesh.dat'
    fileoutput = 'mesh.tec'

# Read mesh type and input (2D,2Daxi,3D)
type = fini_get(ini_file,'PATCHWORK','type')
if type:
    if type=='2D':
        Nz_fix = 2
        width_max = fini_get(ini_file, 'PATCHWORK', 'width', float)
        if not width_max:
            width_max = 1e-3
        zmax = width_max
        stretch = 'null'
    elif type == '2Daxi':
        Nz_fix = 2
        theta_max = fini_get(ini_file, 'PATCHWORK', 'theta', float)
        if not theta_max:
            theta_max = 1*np.pi/180
        zmax = theta_max
        stretch = 'null'
        circular_sector = True
    elif type == '3D':
        Nz_fix = fini_get(ini_file, 'PATCHWORK', 'Nz')
        if Nz_fix:
            Nz_fix = int(Nz_fix)+1
            theta_max = fini_get(ini_file, 'PATCHWORK', 'theta', float)
            width_max = fini_get(ini_file, 'PATCHWORK', 'width', float)
            if theta_max:
                theta_max = theta_max*np.pi/180
                zmax = theta_max
                circular_sector = True
            elif width_max:
                zmax = width_max
            else:
                exit('  ERROR : theta or width parameter necessary!')
            stretch = fini_get(ini_file, 'PATCHWORK', 'stretch')
            if stretch:
                c = fini_get(ini_file, 'PATCHWORK', 'c')
                p = fini_get(ini_file, 'PATCHWORK', 'p')
            else:
                stretch = 'null'
        else:
            threedim = True
            print(f"  Look for a 3D mesh file")
else:
    stretch = 'null'
    threedim = True
    print(f"  Look for a 3D mesh file")

axis = fini_get(ini_file,'PATCHWORK','axis', float)

if not threedim:
    a = np.linspace(0,Nz_fix,num=Nz_fix)/Nz_fix
    if stretch == 'null':
        # No stretch
        zdir = zmax*a-0.5*zmax
    elif stretch == 'left':
        # Left
        zdir = zmax*np.tanh(c*a)/np.tanh(c)-0.5*zmax
    elif stretch == 'outer':
        # Outer
        zdir = zmax*(np.tanh(c*(a-p))+np.tanh(c*p))/(np.tanh(c*(1-p))+np.tanh(c*p))-0.5*zmax

if threedim:
    # Define the source and destination file names
    content_to_remove = "DT=(SINGLE SINGLE)"
    # Open the source file for readin
    with open(fileinput, "r") as source_file:
        # Open the destination file for writing
        with open(fileoutput, "w") as destination_file:
            # Iterate through each line in the source file
            for line in source_file:
                # Check if the line contains the content to remove
                if content_to_remove not in line:
                    # If it doesn't contain the content, write it to the destination file
                    destination_file.write(line)
    print("  Done!")
    exit()

# Read mesh.dat 2D
MESH = open(fileinput,"r")
MESH = MESH.readlines()

Nb = 0

Nx = []
Ny = []
Nz = []

xb = []
yb = []
zb = []

line = 3

while line <= len(MESH)-1:
    
    Nb   += 1
    line += 1
        
    # Extract Nx, Ny, and Nz

    # Define the regular expression pattern to match numbers
    pattern = r'\b\d+\b'
    numbers = re.findall(pattern, MESH[line])
    # Convert the extracted strings to integers
    numbers = [int(num) for num in numbers]

    Nx.append(numbers[0])
    Ny.append(numbers[1])
    Nz.append(numbers[2])
    
    # Skip useless lines
    line += 3

    # Extract X and Y coordinates

    x2d = []
    for j in range(Ny[Nb-1]):
        x1d = []
        for i in range(Nx[Nb-1]):
            x = float(MESH[line])
            x1d.append(x)
            line += 1
        x2d.append(x1d)
    
    xb.append(x2d)
    
    y2d = []
    for j in range(Ny[Nb-1]):
        y1d = []
        for i in range(Nx[Nb-1]):
            y = float(MESH[line])
            if axis and y <= axis:
                y = axis
            y1d.append(y)
            line += 1
        y2d.append(y1d)
    
    yb.append(y2d)

if not threedim:
    Nz[:] = [Nz_fix]*len(Nz)

# Build 3D mesh
print(f"  Building up the 3D grid ...")

Xf = []; Yf = []; Zf = []
XB = []; YB = []; ZB = []
for b in range(Nb):
    X3D = []
    Y3D = []
    Z3D = []
    for k in range(Nz[b]):
        zz = zdir[k]
        X2D = []
        Y2D = []
        Z2D = []
        for j in range(Ny[b]):
            X1D = []
            Y1D = []
            Z1D = []
            for i in range(Nx[b]):
                X1D.append(xb[b][j][i])
                if circular_sector:
                    Y1D.append(yb[b][j][i]*math.cos(zz))
                    Z1D.append(yb[b][j][i]*math.sin(zz))
                else:
                    Y1D.append(yb[b][j][i])
                    Z1D.append(zz)
            X2D.append(X1D)
            Y2D.append(Y1D)
            Z2D.append(Z1D)        
        X3D.append(X2D)
        Y3D.append(Y2D)
        Z3D.append(Z2D)
        zz += zmax/(Nz[b]-1)
    
    # Check for refinement request
    refine = fini_get(ini_file, 'PATCHWORK', 'refine', bool)
    if not refine:
        refine = False
    if refine:
        [Xf,Yf,Zf] = refine_by_splitting(X3D,Y3D,Z3D)
        XB.append(Xf)
        YB.append(Yf)
        ZB.append(Zf)
        Nx[b] = len(Xf[0][0])
        Ny[b] = len(Xf[0])
        Nz[b] = len(Xf)
    else:
        XB.append(X3D)
        YB.append(Y3D)
        ZB.append(Z3D)
    

# Write mesh.tec
MESH = open(fileoutput,"w")
MESH.writelines(' TITLE     = "PATCHWORK mesh" \n')
MESH.writelines(' VARIABLES = "X", "Y", "Z" \n')

for b in range(Nb):
    MESH.writelines(' ZONE T="BLOCCO '+"{0:>1}".format(str(b+1))+'" \n')
    MESH.writelines(' I='+"{0:>3}".format(str(Nx[b]))+", J="+"{0:>3}".format(str(Ny[b]))+", K="+"{0:>3}".format(str(Nz[b]))+", ZONETYPE=Ordered \n")
    MESH.writelines(' DATAPACKING=BLOCK \n')
    # MESH.writelines(' DT=(SINGLE SINGLE) \n')

    for k in range(Nz[b]):
        for j in range(Ny[b]):
            for i in range(Nx[b]):
                MESH.writelines(str(XB[b][k][j][i]/fs)+"\n")

    for k in range(Nz[b]):
        for j in range(Ny[b]):
            for i in range(Nx[b]):
                MESH.writelines(str(YB[b][k][j][i]/fs)+"\n")

    for k in range(Nz[b]):
        for j in range(Ny[b]):
            for i in range(Nx[b]):
                MESH.writelines(str(ZB[b][k][j][i]/fs)+"\n")

print(f"  Done!")


