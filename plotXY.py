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


FreeSurf =1

it = 800
var = 'Vz' #Vy Vz

version = ["FP32-CGFDM","FP16-SVE"] # 
standardVersion = version[0]
plotVersion = version[1]

plotError = 1

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
DT = params["DT"]


standardFileName = "output_" + standardVersion
fileName = [] 
for ver in version:
	if FreeSurf == 1:
		fileName.append("output_" + ver + "/FreeSurf" + "%s_%d"%(var, it) )
	if FreeSurf != 1:
		fileName.append("output_" + ver + "/" + "%s_%d"%(var, it) )


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
	dataX = np.zeros( [grid.NX, grid.NY] )
	dataY = np.zeros( [grid.NX, grid.NY] )
else:
	dataX = np.zeros( [grid.NY, grid.NX] )
	dataY = np.zeros( [grid.NY, grid.NX] )


versionData = { }

for ver in range(len(version)):
	if FAST_AXIS == 'Z':
		data  = np.zeros( [grid.NX, grid.NY] )
	else:
		data  = np.zeros( [grid.NY, grid.NX] )
	for mpiY in range( grid.PY ):
		for mpiX in range( grid.PX ):
			mpiZ = mpiSliceZ
			XFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
			YFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
			if FreeSurf:
				mpiZ = grid.PZ - 1
			else:
				mpiZ = mpiSliceZ
			print( fileName[ver] )
			File  = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileName[ver], mpiX, mpiY, mpiZ ), "rb" )
			
			ny = grid.ny[mpiY]
			nx = grid.nx[mpiX]
	
			print( "ny = %d, nx = %d" % ( nx, ny ) )
			datax = np.fromfile( XFile, dtype='float32', count = ny * nx )
			datay = np.fromfile( YFile, dtype='float32', count = ny * nx )
			data_ = np.fromfile(  File, dtype='float32', count = ny * nx )
	
			J  = grid.frontNY[mpiY]
			J_ = grid.frontNY[mpiY] + ny
			I  = grid.frontNX[mpiX]
			I_ = grid.frontNX[mpiX] + nx
	
			if FAST_AXIS == 'Z':
				dataX[I:I_, J:J_] = np.reshape( datax, (nx, ny) )
				dataY[I:I_, J:J_] = np.reshape( datay, (nx, ny) )
				data [I:I_, J:J_] = np.reshape( data_, (nx, ny) )
			else:
				dataX[J:J_, I:I_] = np.reshape( datax, (ny, nx) )
				dataY[J:J_, I:I_] = np.reshape( datay, (ny, nx) )
				data [J:J_, I:I_] = np.reshape( data_, (ny, nx) )

	vm = np.max( np.abs( data ) )
	print("%s: %f"%(version[ver], vm))
	versionData[version[ver]] = data
	vm = np.max( np.abs( versionData[version[ver]] ) )
	print("%s: %f"%(version[ver], vm))


v1 = np.max( np.abs( versionData[plotVersion] ) )
v2 = np.max( np.abs( versionData[standardVersion] ) )
print("%s: %f; %s %f"%(plotVersion, v1,  standardVersion, v2))

plt.figure( 2 )
unit = 1000 # 1km = 1000m
dataX = dataX/unit
dataY = dataY/unit
vm = np.max( np.abs( versionData[standardVersion] ) )
print( "vm = %f" % vm   )

gain = 2


simple = 2
dataX = dataX[::simple,::simple]
dataY = dataY[::simple,::simple]

if plotError != 1:
	data = versionData[plotVersion]
	data = data[::simple,::simple]
	plt.pcolormesh( dataX, dataY, data, vmax = vm/gain, vmin = -vm/gain, cmap = "seismic" )

if plotError == 1:
	data = np.abs(versionData[plotVersion] - versionData[standardVersion])/vm * 100
	tvm = np.max( data )
	print("tvm = %f" % tvm )
	data = data[::simple,::simple]
	plt.pcolormesh( dataX, dataY, data, vmax = tvm, vmin = 0.0, cmap = "Reds" )



Xmin = np.min( dataX )
Xmax = np.max( dataX )
Ymin = np.min( dataY )
Ymax = np.max( dataY )



xticksRange = [ -120, -60, 0, 60, 120 ]
yticksRange = [ -120, -60, 0, 60, 120 ] 

plt.xticks(xticksRange)
plt.yticks(yticksRange)

plt.tick_params(axis='both',which='major')

plt.ylabel( "South-North (km)" )
plt.xlabel( "West-East (km)" )


plt.axis( "image" )

if var == 'Vx':
	varLow = "$v_x$"
if var == 'Vy':
	varLow = "$v_y$"
if var == 'Vz':
	varLow = "$v_z$"


if plotError == 1:
	cb = plt.colorbar( format = '%.2f', shrink = 0.9 )
	cb.ax.set_title( "%")
	cbRange = np.linspace( 0, tvm, 5 )
	cb.set_ticks( cbRange )
	plt.title( "%s: Relative Errors at time = %.2fs" % ( plotVersion, it * DT ), fontsize = 10 )
	#plt.gca( ).set_ylim( [Zmin, Zmax + 0.5 ] )
	plt.savefig( "png/ReError_%s_%s_%d.pdf"%(plotVersion, var, it), bbox_inches='tight' )
else:
	cb = plt.colorbar( format = '%.1e', shrink = 0.9 )
	cb.ax.set_title( "m/s" )
	cbRange = np.linspace( -vm/gain, vm/gain, 5 )
	cb.set_ticks( cbRange )
	plt.title( "%s: %s at time = %.2fs" % ( plotVersion, varLow, it * DT ), fontsize = 10 )
	#plt.gca( ).set_ylim( [Zmin, Zmax + 0.5 ] )
	plt.savefig( "png/%s_%s_%d.pdf"%(plotVersion, var, it), bbox_inches='tight' )



