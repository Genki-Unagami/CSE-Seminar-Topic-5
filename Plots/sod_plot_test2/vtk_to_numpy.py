#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 10:07:36 2020

@author: genki
"""

import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy


### Reference ###
exact_data = np.loadtxt("./t=0.150000[step2062].csv")

### Osher ###
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("./errors-15.vtu")
reader.Update()
data = reader.GetOutput()

points = data.GetPoints()
npts = points.GetNumberOfPoints()
x = vtk_to_numpy(points.GetData())
cell = (vtk_to_numpy(data.GetCells().GetConnectivityArray())).reshape(90000,4)
q = vtk_to_numpy(data.GetCellData().GetArray(0))

output_osher = np.zeros([300,4])
tmp = np.where(x[:,1]==0.5)[0]

for i in range(25):
    for j in range(13):
        if (j != 12):
            pointindex = 13 * i + j
            cellindex = np.where(cell[:,0] == tmp[pointindex])[0]
            loc = x[tmp[pointindex],0]
            
            ### Numerical Result ###
            rho = q[cellindex,0]
            u = q[cellindex,1] / rho
            p = (1.4-1)*(q[cellindex,4]-0.5*rho*(u**2))
            
            output_osher[12*i+j,0] = loc
            output_osher[12*i+j,1] = rho
            output_osher[12*i+j,2] = u
            output_osher[12*i+j,3] = p

#### Osher ###
#reader = vtk.vtkXMLUnstructuredGridReader()
#reader.SetFileName("./errors-12_osher.vtu")
#reader.Update()
#data = reader.GetOutput()
#
#points = data.GetPoints()
#npts = points.GetNumberOfPoints()
#x = vtk_to_numpy(points.GetData())
#cell = (vtk_to_numpy(data.GetCells().GetConnectivityArray())).reshape(90000,4)
#q = vtk_to_numpy(data.GetCellData().GetArray(0))
#
#output_osher = np.zeros([300,4])
#tmp = np.where(x[:,1]==0.5)[0]
#
#for i in range(25):
#    for j in range(13):
#        if (j != 12):
#            pointindex = 13 * i + j
#            cellindex = np.where(cell[:,0] == tmp[pointindex])[0]
#            loc = x[tmp[pointindex],0]
#            
#            ### Numerical Result ###
#            rho = q[cellindex,0]
#            u = q[cellindex,1] / rho
#            p = (1.4-1)*(q[cellindex,4]-0.5*rho*(u**2))
#            
#            output_osher[12*i+j,0] = loc
#            output_osher[12*i+j,1] = rho
#            output_osher[12*i+j,2] = u
#            output_osher[12*i+j,3] = p
            
            
            


fig1, ax1 = plt.subplots(figsize=(8, 6))
fig2, ax2 = plt.subplots(figsize=(8, 6))
fig3, ax3 = plt.subplots(figsize=(8, 6))

#ax1.scatter(output_roe[:,0], output_roe[:,1],facecolors='none',edgecolor='blue',label='roe')
ax1.scatter(output_osher[:,0], output_osher[:,1],facecolors='none',edgecolor='orange',label='osher')
ax1.plot(exact_data[:,0],exact_data[:,1],'k',label='reference')


#ax2.scatter(output_roe[:,0], output_roe[:,2],facecolors='none',edgecolor='blue',label='roe')
ax2.scatter(output_osher[:,0], output_osher[:,2],facecolors='none',edgecolor='orange',label='osher')
ax2.plot(exact_data[:,0],exact_data[:,4],'k',label='reference')

#ax3.scatter(output_roe[:,0], output_roe[:,3],facecolors='none',edgecolor='blue',label='roe')
ax3.scatter(output_osher[:,0], output_osher[:,3],facecolors='none',edgecolor='orange',label='osher')
ax3.plot(exact_data[:,0],exact_data[:,5],'k',label='reference')

ax1.set_xlabel(r'$x$', fontsize=18)
ax2.set_xlabel(r'$x$', fontsize=18)
ax3.set_xlabel(r'$x$', fontsize=18)
ax1.set_ylabel(r'$\rho$', fontsize=18)
ax2.set_ylabel(r'$u$', fontsize=18)
ax3.set_ylabel(r'$p$', fontsize=18)

ax1.legend(fontsize=14)
ax2.legend(fontsize=14)
ax3.legend(fontsize=14)
plt.show()