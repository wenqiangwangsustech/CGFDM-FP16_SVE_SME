#!/usr/bin/env python

'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import  sys
from scipy.io import savemat

FreeSurf = 1

it = 3000

var = 'Vy' #Vy Vz

PGVh = 1

if len( sys.argv ) > 1:
	it = int( sys.argv[1] )



PGV = 0
if PGVh == 1:
	PGV = 0
if PGVh == 0:
	PGV = 1

#jsonsFile = open( "params.json" )
#jsonsFile = open( "paramsDir/paramsCGFDM3D.json" )
#jsonsFile = open( "paramsDir/paramsCGFDM3D-CJMVS.json" )
#jsonsFile = open( "paramsDir/paramsCGFDM3D-LSRK.json" )
jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
FAST_AXIS = params["FAST_AXIS"]

sample = 1

outputPath = params["out"]
fileNameX = params["out"] + "/lon"
fileNameY = params["out"] + "/lat"


if FreeSurf == 1:
	#fileName = params["out"] + "/FreeSurf%s_%d" % ( var, it )
    fileName = params["out"] + "/PGV"
else:
	fileName = params["out"] + "/%s_%d" % ( var, it )


PGVSIZE = 2


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
	dataX = np.zeros( [grid.NX, grid.NY] )
	dataY = np.zeros( [grid.NX, grid.NY] )
	data  = np.zeros( [grid.NX, grid.NY] )
	Intensity = np.zeros( [grid.NX, grid.NY] )
else:
	dataX = np.zeros( [grid.NY, grid.NX] )
	dataY = np.zeros( [grid.NY, grid.NX] )
	data  = np.zeros( [grid.NY, grid.NX] )
	Intensity = np.zeros( [grid.NY, grid.NX] )


for mpiY in range( grid.PY ):
	for mpiX in range( grid.PX ):
		mpiZ = mpiSliceZ#grid.PZ - 1
		XFile = open( "%s_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		YFile = open( "%s_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
		if FreeSurf:
			mpiZ = grid.PZ - 1
		else:
			mpiZ = mpiSliceZ
		File  = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		
		ny = grid.ny[mpiY]
		nx = grid.nx[mpiX]

		print( "ny = %d, nx = %d" % ( nx, ny ) )
		datax = np.fromfile( XFile, dtype='float32', count = ny * nx )
		datay = np.fromfile( YFile, dtype='float32', count = ny * nx )
		data_ = np.fromfile(  File, dtype='float32', count = ny * nx * PGVSIZE )
		dataTmp = np.reshape( data_, ( ny * nx, PGVSIZE ) )
		J  = grid.frontNY[mpiY]
		J_ = grid.frontNY[mpiY] + ny
		I  = grid.frontNX[mpiX]
		I_ = grid.frontNX[mpiX] + nx

		if FAST_AXIS == 'Z':
			dataX[I:I_, J:J_] = np.reshape( datax, (nx, ny) )
			dataY[I:I_, J:J_] = np.reshape( datay, (nx, ny) )
			data [I:I_, J:J_] = np.reshape( dataTmp[:, PGVh], (nx, ny) )
		else:
			dataX[J:J_, I:I_] = np.reshape( datax, (ny, nx) )
			dataY[J:J_, I:I_] = np.reshape( datay, (ny, nx) )
			data [J:J_, I:I_] = np.reshape( dataTmp[:, PGVh], (ny, nx) )

PGV = data

for j in range( grid.NY ):
	for i in range( grid.NX ):
		if  PGV[j,i] >= 1.41:
			Intensity[j,i] = 11
		if  PGV[j,i] >= 0.72 and  PGV[j,i] < 1.41:
			Intensity[j,i] = 10
		if  PGV[j,i] >= 0.36 and  PGV[j,i] < 0.72:
			Intensity[j,i] = 9
		if  PGV[j,i] >= 0.19 and  PGV[j,i] < 0.36:
			Intensity[j,i] = 8
		if  PGV[j,i] >= 0.10 and  PGV[j,i] < 0.19:
			Intensity[j,i] = 7
		if  PGV[j,i] >= 0.05 and  PGV[j,i] < 0.1:
			Intensity[j,i] = 6
		if  PGV[j,i] >= 0.02 and  PGV[j,i] < 0.05:
			Intensity[j,i] = 5
		if  PGV[j,i] >= 0.01 and  PGV[j,i] < 0.02:
			Intensity[j,i] = 4
		if  PGV[j,i] >= 0.005 and  PGV[j,i] < 0.01:
			Intensity[j,i] = 3
		if  PGV[j,i] >= 0.001 and  PGV[j,i] < 0.005:
			Intensity[j,i] = 2
		if  PGV[j,i] < 0.001:
			Intensity[j,i] = 1



nPML = params["nPML"]
NX = grid.NX
NY = grid.NY
PGVh = PGV[nPML:NY - nPML, nPML:NX-nPML]
np.save( "./numpyData/PGVh.npy", PGVh )
logPGVh = np.log( PGV[nPML:NY - nPML, nPML:NX-nPML] )
np.save( "./numpyData/logPGVh.npy", logPGVh )

Intensity = Intensity[nPML:NY - nPML, nPML:NX-nPML]
lon = dataX[nPML:NY - nPML, nPML:NX-nPML]
lat = dataY[nPML:NY - nPML, nPML:NX-nPML]

print(np.max(Intensity))

np.save( "./numpyData/Intensity.npy", Intensity )
np.save( "./numpyData/lon.npy", lon )
np.save( "./numpyData/lat.npy", lat )
    
savemat('./numpyData/lon.mat', {'lon':lon})
savemat('./numpyData/lat.mat', {'lat':lat})
savemat('./numpyData/PGVh.mat', {'PGVh':PGVh})
savemat('./numpyData/logPGVh.mat', {'logPGVh':logPGVh})

#Intensity = 3.00 * np.log10(PGV) + 9.77

plt.figure( 2 )

dpi = 300
#vm = np.max( np.abs( PGV ) ) / 2
plt.pcolormesh( lon, lat, Intensity, cmap = "seismic" ) #origin= "lower" )
plt.colorbar()
#plt.pcolormesh( dataX, dataY, data, cmap = "jet" ) #origin= "lower" )

plt.axis( "image" )
plt.title( 'Intensity' )

plt.show( )
plt.savefig( fileName + ".png" )
#plt.plot( dataY[grid.NZ - 1, :] )

#print( grid )


