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
exact_data = np.loadtxt("t=1.800000[step4282].csv")

### Osher ###
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("entropy_osher.vtu")
reader.Update()
data = reader.GetOutput()

points = data.GetPoints()
npts = points.GetNumberOfPoints()
x = vtk_to_numpy(points.GetData())
cell = (vtk_to_numpy(data.GetCells().GetConnectivityArray())).reshape(91008,4)
q = vtk_to_numpy(data.GetCellData().GetArray(0))

output_osher = np.zeros([948,4])
tmp = np.where(x[:,1]==0.0105485)[0]

for i in range(79):
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
            
### Roe ###
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("entropy_roe.vtu")
reader.Update()
data = reader.GetOutput()

points = data.GetPoints()
npts = points.GetNumberOfPoints()
x = vtk_to_numpy(points.GetData())
cell = (vtk_to_numpy(data.GetCells().GetConnectivityArray())).reshape(91008,4)
q = vtk_to_numpy(data.GetCellData().GetArray(0))

output_roe = np.zeros([948,4])
tmp = np.where(x[:,1]==0.0105485)[0]

for i in range(79):
    for j in range(13):
        if (j != 12):
            pointindex = 13 * i + j
            cellindex = np.where(cell[:,0] == tmp[pointindex])[0]
            loc = x[tmp[pointindex],0]
            
            ### Numerical Result ###
            rho = q[cellindex,0]
            u = q[cellindex,1] / rho
            p = (1.4-1)*(q[cellindex,4]-0.5*rho*(u**2))
            
            output_roe[12*i+j,0] = loc
            output_roe[12*i+j,1] = rho
            output_roe[12*i+j,2] = u
            output_roe[12*i+j,3] = p

### Rusanov ###
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("entropy_rusanov.vtu")
reader.Update()
data = reader.GetOutput()

points = data.GetPoints()
npts = points.GetNumberOfPoints()
x = vtk_to_numpy(points.GetData())
cell = (vtk_to_numpy(data.GetCells().GetConnectivityArray())).reshape(91008,4)
q = vtk_to_numpy(data.GetCellData().GetArray(0))

output_rusanov = np.zeros([948,4])
tmp = np.where(x[:,1]==0.0105485)[0]

for i in range(79):
    for j in range(13):
        if (j != 12):
            pointindex = 13 * i + j
            cellindex = np.where(cell[:,0] == tmp[pointindex])[0]
            loc = x[tmp[pointindex],0]
            
            ### Numerical Result ###
            rho = q[cellindex,0]
            u = q[cellindex,1] / rho
            p = (1.4-1)*(q[cellindex,4]-0.5*rho*(u**2))
            
            output_rusanov[12*i+j,0] = loc
            output_rusanov[12*i+j,1] = rho
            output_rusanov[12*i+j,2] = u
            output_rusanov[12*i+j,3] = p
            
#exact_sorted = exact[np.argsort(exact[:,0])]

fig1, ax1 = plt.subplots(figsize=(8, 6))
fig2, ax2 = plt.subplots(figsize=(8, 6))
fig3, ax3 = plt.subplots(figsize=(8, 6))

fig4, ax4 = plt.subplots(figsize=(8, 6))
fig5, ax5 = plt.subplots(figsize=(8, 6))
fig6, ax6 = plt.subplots(figsize=(8, 6))

ax1.scatter(output_rusanov[:,0], output_rusanov[:,1],facecolors='none',edgecolor='red',label='rusanov')
ax1.scatter(output_roe[:,0], output_roe[:,1],facecolors='none',edgecolor='blue',label='roe')
ax4.scatter(output_rusanov[:,0], output_rusanov[:,1],facecolors='none',edgecolor='red',label='rusanov')
ax4.scatter(output_osher[:,0], output_osher[:,1],facecolors='none',edgecolor='orange',label='osher')
ax1.plot(exact_data[:,0],exact_data[:,1],'k',label='reference')
ax4.plot(exact_data[:,0],exact_data[:,1],'k',label='reference')
#ax1.plot(exact_sorted[:,0], exact_sorted[:,1],'k',label='exact')

ax2.scatter(output_rusanov[:,0], output_rusanov[:,2],facecolors='none',edgecolor='red',label='rusanov')
ax2.scatter(output_roe[:,0], output_roe[:,2],facecolors='none',edgecolor='blue',label='roe')
ax5.scatter(output_rusanov[:,0], output_rusanov[:,2],facecolors='none',edgecolor='red',label='rusanov')
ax5.scatter(output_osher[:,0], output_osher[:,2],facecolors='none',edgecolor='orange',label='osher')
ax2.plot(exact_data[:,0],exact_data[:,4],'k',label='reference')
ax5.plot(exact_data[:,0],exact_data[:,4],'k',label='reference')
#ax2.plot(exact_sorted[:,0], exact_sorted[:,2],'k',label='exact')
#
ax3.scatter(output_rusanov[:,0], output_rusanov[:,3],facecolors='none',edgecolor='red',label='rusanov')
ax3.scatter(output_roe[:,0], output_roe[:,3],facecolors='none',edgecolor='blue',label='roe')
ax6.scatter(output_rusanov[:,0], output_rusanov[:,3],facecolors='none',edgecolor='red',label='rusanov')
ax6.scatter(output_osher[:,0], output_osher[:,3],facecolors='none',edgecolor='orange',label='osher')
ax3.plot(exact_data[:,0],exact_data[:,5],'k',label='reference')
ax6.plot(exact_data[:,0],exact_data[:,5],'k',label='reference')
#ax3.plot(exact_sorted[:,0], exact_sorted[:,3],'k',label='exact')

ax1.set_xlabel(r'$x$', fontsize=18)
ax2.set_xlabel(r'$x$', fontsize=18)
ax3.set_xlabel(r'$x$', fontsize=18)
ax1.set_ylabel(r'$\rho$', fontsize=18)
ax2.set_ylabel(r'$u$', fontsize=18)
ax3.set_ylabel(r'$p$', fontsize=18)
ax4.set_xlabel(r'$x$', fontsize=18)
ax5.set_xlabel(r'$x$', fontsize=18)
ax6.set_xlabel(r'$x$', fontsize=18)
ax4.set_ylabel(r'$\rho$', fontsize=18)
ax5.set_ylabel(r'$u$', fontsize=18)
ax6.set_ylabel(r'$p$', fontsize=18)

ax1.legend(loc=3, fontsize=14)
ax2.legend(loc=3, fontsize=14)
ax3.legend(loc=3, fontsize=14)
ax4.legend(loc=3, fontsize=14)
ax5.legend(loc=3, fontsize=14)
ax6.legend(loc=3, fontsize=14)
plt.show()