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


it = 1000
var = 'Vx' #Vx Vz
mod = 'FP16+SVE'



if len( sys.argv ) > 1:
	it = int( sys.argv[1] )
	var = str( sys.argv[2] )




jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

sample = 1

outputPath = params["out"]
fileNameX = params["out"] + "/coordX"
fileNameY = params["out"] + "/coordY"
fileNameZ = params["out"] + "/coordZ"


varname = "./png/%s_%d_%s" % ( var, it, mod )

fileName = params["out"] + "/%s_%d"%( var, it )
#fileName = params["out"] + "/Vs"
#fileName = "output" +"/%s_%d"%( var, it )
print( "Draw " + fileName  )


FAST_AXIS = params["FAST_AXIS"]

DT = params["DT"]
'''
for Z in range( grid.PZ ):
	for Y in range( grid.PY ):
		for X in range( grid.PX ):
			print( "nx = %d, ny = %d, nz = %d\n" % ( grid.nx[X], grid.ny[Y], grid.nz[Z] ) )
'''

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
	data  = np.zeros( [grid.NX, grid.NZ] )
else:
	dataX = np.zeros( [grid.NZ, grid.NX] )
	dataZ = np.zeros( [grid.NZ, grid.NX] )
	data  = np.zeros( [grid.NZ, grid.NX] )


mpiY = mpiSliceY
for mpiZ in range( grid.PZ ):
	for mpiX in range( grid.PX ):
		fileX = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		fileZ = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameZ, mpiX, mpiY, mpiZ ), "rb" )
		file  = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		nx = grid.nx[mpiX]
		nz = grid.nz[mpiZ]
		print( "nx = %d, nz = %d" % ( nx, nz ) )
		datax = np.fromfile( fileX, dtype='float32', count = nx * nz )
		#print( np.shape( datax ) )
		dataz = np.fromfile( fileZ, dtype='float32', count = nx * nz )
		data_  = np.fromfile( file, dtype='float32', count = nx * nz )
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

plt.figure(2, dpi = 300)
dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1000 # 1km = 1000m
dataX = dataX/unit
dataZ = dataZ/unit
#plt.pcolormesh( dataX, dataZ, data cmap = "jet" )


#np.save( "./numpyData/%s_%s_%d.npy"%(mod, var, it), data)
dataFP32 = np.load("./numpyData/%s_%s_%d.npy"%(mod, var, it))
data = np.abs( dataFP32 - data )/2
vm = np.max( np.abs( data ) ) 
print( "vm = %f" % vm   )
#plt.pcolormesh( dataX, dataZ, data, vmax = vm, vmin = -vm, cmap = "seismic" )
plt.pcolormesh( dataX, dataZ, data, vmax = vm, vmin = 0, cmap = "Reds" )


#plt.pcolormesh( dataX, dataZ, data, cmap = "jet" )

plt.plot( dataX[:, -1], dataZ[:, -1],  'k', linewidth = 2  )
plt.plot( dataX[:,  0], dataZ[:,  0],  'k', linewidth = 2  )
plt.plot( dataX[-1, :], dataZ[-1, :],  'k', linewidth = 2  )
plt.plot( dataX[ 0, :], dataZ[ 0, :],  'k', linewidth = 2  )



Xmin = np.min( dataX )
Xmax = np.max( dataX )
Zmin = np.min( dataZ )
Zmax = np.max( dataZ )
xticksRange = [ -10, -5, 0, 5, 10 ]
 
zticksRange = [ -20, -16, -12, -8, -4, 0, 3 ] #np.linspace( -10, 2.5, 5 )

plt.xticks(xticksRange)
plt.yticks(zticksRange)



plt.tick_params(axis='both',which='major')

plt.ylabel( "Vertical Distance-Z (km)" )
plt.xlabel( "Horizontal Distance-X (km)" )


#plt.pcolormesh( data[::sample, ::sample] // unit, vmin = -vm / 2, vmax = vm, cmap = "jet" )
#plt.pcolormesh( data[5:grid.NZ -5:sample, 5:grid.NX - 5:sample] // unit, cmap = "seismic" )
#plt.pcolormesh( data, vmax = vm, vmin = -vm, cmap = "jet" )
#plt.imshow( data, cmap = "jet", origin= "lower" )
cb = plt.colorbar( format = '%.1e', shrink = 0.8 )
#plt.colorbar( orientation = "horizontal" )
plt.axis( "image" )
#plt.plot( dataZ[grid.NZ - 1, :] )

#print( grid )
if var == 'Vx':
	varLow = "$v_x$"
if var == 'Vy':
	varLow = "$v_y$"
if var == 'Vz':
	varLow = "$v_z$"

cb.ax.set_title( "%s (m/s)" % varLow )
cbRange = np.linspace( 0, vm, 5 )
#cbRange = np.linspace( -vm, vm, 5 )
cb.set_ticks( cbRange )
plt.title( "%s Error: time = %.2fs" % ( mod, it * DT ), fontsize = 14 )
#plt.title( "Absolute Errors (x10): time = %.2fs" % ( it * DT ), fontsize = 14 )
plt.savefig( varname + "Error.png" )



