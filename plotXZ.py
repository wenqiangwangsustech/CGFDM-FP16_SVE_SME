#!/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import sys


it = 600
var = 'Vx' #Vx Vz
plotError = 0


Keys = ['A', 'B', 'C', 'D', 'E']


version = ["FP32-CGFDM", "FP16-CGFDM", "FP32-SVE", "FP16-SVE"]
standardVersion = version[0]
plotVersion = version[3]


ModelName = "Gaussian-shaped Topography Model"

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
DT = params["DT"]
DH = params["DH"]



standardFileName = "output_" + standardVersion
fileName = [] 
for ver in version:
	fileName.append("output_" + ver + "/" + "%s_%d"%(var, it) )

#coordinate file
fileNameX = standardFileName + "/coordX"
fileNameY = standardFileName + "/coordY"
fileNameZ = standardFileName + "/coordZ"

FAST_AXIS = params["FAST_AXIS"]

sliceX = params["sliceX"] - grid.frontNX
sliceY = params["sliceY"] - grid.frontNY
sliceZ = params["sliceZ"] - grid.frontNZ

for mpiSliceX in range( grid.PX ):
	if sliceX[mpiSliceX] >= 0 and sliceX[mpiSliceX] < grid.nx[mpiSliceX]:
		break

for mpiSliceY in range( grid.PY ):
	if sliceY[mpiSliceY] >= 0 and sliceY[mpiSliceY] < grid.ny[mpiSliceY]:
		break

for mpiSliceZ in range( grid.PZ ):
	if sliceZ[mpiSliceZ] >= 0 and sliceZ[mpiSliceZ] < grid.nz[mpiSliceZ]:
		break

if FAST_AXIS == 'Z':
	dataX = np.zeros( [grid.NX, grid.NZ] )
	dataZ = np.zeros( [grid.NX, grid.NZ] )
else:
	dataX = np.zeros( [grid.NZ, grid.NX] )
	dataZ = np.zeros( [grid.NZ, grid.NX] )

versionData = { }

for ver in range(len(version)):
	if FAST_AXIS == 'Z':
		data  = np.zeros( [grid.NX, grid.NZ] )
	else:
		data  = np.zeros( [grid.NZ, grid.NX] )
	mpiY = mpiSliceY
	for mpiZ in range( grid.PZ ):
		for mpiX in range( grid.PX ):
			fileX = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
			fileZ = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameZ, mpiX, mpiY, mpiZ ), "rb" )
			file  = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileName[ver], mpiX, mpiY, mpiZ ), "rb" )
			#print(fileName[ver])
			nx = grid.nx[mpiX]
			nz = grid.nz[mpiZ]
			#print( "nx = %d, nz = %d" % ( nx, nz ) )
			datax = np.fromfile( fileX, dtype='float32', count = nx * nz )
			#print( np.shape( datax ) )
			dataz = np.fromfile( fileZ, dtype='float32', count = nx * nz )
			data_  = np.fromfile( file, dtype='float32', count = nx * nz )
			#print("max = %f\n"%np.abs(np.max(data_)))
			I  = grid.frontNX[mpiX]
			I_ = grid.frontNX[mpiX] + nx
			K  = grid.frontNZ[mpiZ]
			K_ = grid.frontNZ[mpiZ] + nz
	
			if FAST_AXIS == 'Z':
				dataX[I:I_, K:K_] = np.reshape( datax, ( nx, nz ) )
				dataZ[I:I_, K:K_] = np.reshape( dataz, ( nx, nz ) )
				data [I:I_, K:K_] = np.reshape( data_, ( nx, nz ) )
			else:
				dataX[K:K_, I:I_] = np.reshape( datax, ( nz, nx ) )
				dataZ[K:K_, I:I_] = np.reshape( dataz, ( nz, nx ) )
				data [K:K_, I:I_] = np.reshape( data_, ( nz, nx ) )
	print(version[ver])
	versionData[version[ver]] = data


plt.figure(2, dpi = 300)
unit = 1000 # 1km = 1000m
dataX = dataX/unit
dataZ = dataZ/unit
vm = np.max( np.abs( versionData[standardVersion] ) )
print( "vm = %f" % vm   )

gain = 2


if plotError != 1:
	data = versionData[plotVersion]
	plt.pcolormesh( dataX, dataZ, data, vmax = vm/gain, vmin = -vm/gain, cmap = "seismic" )

if plotError == 1:
	data = np.abs(versionData[plotVersion] - versionData[standardVersion])/vm * 100
	tvm = np.max( data )
	print("tvm = %f" % tvm )
	plt.pcolormesh( dataX, dataZ, data, vmax = tvm, vmin = 0.0, cmap = "Reds" )


plt.plot( dataX[:, -1], dataZ[:, -1],  'k', linewidth = 2  )
plt.plot( dataX[:,  0], dataZ[:,  0],  'k', linewidth = 2  )
plt.plot( dataX[-1, :], dataZ[-1, :],  'k', linewidth = 2  )
plt.plot( dataX[ 0, :], dataZ[ 0, :],  'k', linewidth = 2  )

Xmin = np.min( dataX )
Xmax = np.max( dataX )
Zmin = np.min( dataZ )
Zmax = np.max( dataZ )
xticksRange = [ -15, -10, -5, 0, 5, 10, 15 ]
zticksRange = [ -20, -15, -10, -5, 0,  3 ] 

plt.xticks(xticksRange)
plt.yticks(zticksRange)

plt.tick_params(axis='both',which='major')

plt.ylabel( "Vertical Distance-Z (km)" )
plt.xlabel( "Horizontal Distance-X (km)" )

plt.axis( "image" )

if var == 'Vx':
	varLow = "$v_x$"
if var == 'Vy':
	varLow = "$v_y$"
if var == 'Vz':
	varLow = "$v_z$"


sourceX = params['sourceX']
sourceY = params['sourceY']
sourceZ = params['sourceZ']

centerX = params['centerX']
centerY = params['centerY']



stationFile = open( "station.json" )
stationDir = json.load( stationFile )

station = stationDir["station(point)"] 


srcX = ( sourceX - centerX ) * DH 
srcY = ( sourceY - centerY ) * DH 
srcZ = - ( grid.NZ - 1 - sourceZ ) * DH 

print("srcX = %f, srcY = %f, srcZ = %f"%(srcX, srcY, srcZ))


recvX = []
recvY = []
recvZ = []

for key_it in range(len(Keys)):
	value = station[Keys[key_it]]
	recvX.append( ( value[0] - centerX ) * DH  )
	recvY.append( ( value[1] - centerY ) * DH  )
	recvZ.append( dataZ[-1, value[0]]*unit - (grid.NZ - 1 - value[2]) * DH )
	

lenRecv = len( recvX )

print( recvX )
print( recvZ )


if plotError == 1:
	cb = plt.colorbar( format = '%.2f', shrink = 0.65 )
	cb.ax.set_title( "%")
	cbRange = np.linspace( 0, tvm, 5 )
	cb.set_ticks( cbRange )
	plt.title( "%s: Relative Errors at time = %.2fs" % ( plotVersion, it * DT ), fontsize = 10 )
	#plt.gca( ).set_ylim( [Zmin, Zmax + 0.5 ] )
	plt.savefig( "png/ReError_%s_%s_%d.pdf"%(plotVersion, var, it), bbox_inches='tight' )
else:
	cb = plt.colorbar( format = '%.1e', shrink = 0.65 )
	cb.ax.set_title( "m/s" )
	plt.plot( srcX / unit, srcZ / unit, 'r*', markersize=10 )
	for iRecv in range( lenRecv ):
		#print( Keys[iRecv] + ":" + "x = %f, z = %f" % (recvX[iRecv] / unit, recvZ[iRecv] / unit )  )
		if Keys[iRecv] == 'C':
			plt.text( recvX[iRecv] / unit + 0.5, recvZ[iRecv] / unit-0.2, Keys[iRecv] )
			plt.plot( recvX[iRecv] / unit, recvZ[iRecv] / unit,  'g' + 'v', markersize=5 )
		if Keys[iRecv] != 'C':
			plt.text( recvX[iRecv] / unit -0.4, recvZ[iRecv] / unit + 0.5, Keys[iRecv] )
			plt.plot( recvX[iRecv] / unit, recvZ[iRecv] / unit,  'g' + 'v', markersize=5 )

	cbRange = np.linspace( -vm/gain, vm/gain, 5 )
	cb.set_ticks( cbRange )
	plt.title( "%s: %s at time = %.2fs" % ( plotVersion, varLow, it * DT ), fontsize = 10 )
	#plt.gca( ).set_ylim( [Zmin, Zmax + 0.5 ] )
	plt.savefig( "png/%s_%s_%d.pdf"%(plotVersion, var, it), bbox_inches='tight' )


