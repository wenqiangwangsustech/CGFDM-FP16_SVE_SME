#!/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys



def pltFormat( plot, var, grid, mpiSliceY, it, fileNameX, fileNameZ, fileName, dataX, dataZ, data ):
	plot.clf( )
	print( "it = %d" % it )
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
	#plt.figure(1)
	vm = np.max( np.abs( data ) ) / 2
	plot.pcolormesh( dataX, dataZ, data, vmax = vm, vmin = -vm, cmap = "seismic" )
	cb = plot.colorbar( )
	cb.ax.tick_params( labelsize = 25 )
	plot.axis( "image" )
	plot.title( var + "%d" % it, fontsize = 25 )
	

it = 1600
var = 'Vx' #Vx Vz

itStart = 0
itEnd = 1220



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


varname = "./png/%s_%d" % ( var, it )

#fileName = params["out"] + "/Vs"
fileName = "output" +"/%s_%d"%( var, it )
print( "Draw " + fileName  )


FAST_AXIS = params["FAST_AXIS"]


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




TMAX = params["TMAX"] 
DT = params["DT"]
IT_SKIP = params["IT_SKIP"]

NT = int( TMAX / DT )


fig = plt.figure( 1, figsize = ( 16, 16 ) )
def update( it ):
	print( "iterator: %d" % it)
	fileName = params["out"] + "/%s_%d"%( var, it )
	pltFormat( plt, var, grid, mpiSliceY, it, fileNameX, fileNameZ, fileName, dataX, dataZ, data )
	return fig

anim = animation.FuncAnimation( fig, update, frames = range( itStart, itEnd, IT_SKIP ), interval = 50 )
anim.save( varname + ".gif", writer = 'pillow', fps = 100 )


