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



version = ["FP32-CGFDM","FP16-SVE"] # 
standardVersion = version[0]
plotVersion = version[0]

plotError = 0

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
FAST_AXIS = params["FAST_AXIS"]


PGVh = 0
PGV = 1
if PGVh == 1:
	PGV = 0
if PGVh == 0:
	PGV = 1



standardFileName = "output_" + standardVersion
fileName = [] 
for ver in version:
	fileName.append("output_" + ver + "/PGV" )


fileNameX = params["out"] + "/lon"
fileNameY = params["out"] + "/lat"
fileNameTerrain = params["out"] + "/terrain"


PGVSIZE = 2

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
	Terr  = np.zeros( [grid.NX, grid.NY] )
	Intensity = np.zeros( [grid.NX, grid.NY] )
else:
	dataX = np.zeros( [grid.NY, grid.NX] )
	dataY = np.zeros( [grid.NY, grid.NX] )
	Terr  = np.zeros( [grid.NY, grid.NX] )
	Intensity = np.zeros( [grid.NY, grid.NX] )


versionData = { }

for ver in range(len(version)):
	if FAST_AXIS == 'Z':
		data  = np.zeros( [grid.NX, grid.NY] )
	else:
		data  = np.zeros( [grid.NY, grid.NX] )
	for mpiY in range( grid.PY ):
		for mpiX in range( grid.PX ):
			mpiZ = mpiSliceZ
			XFile = open( "%s_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
			YFile = open( "%s_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
			mpiZ = grid.PZ - 1
			File  = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileName[ver], mpiX, mpiY, mpiZ ), "rb" )
			TerFile = open( "%s_mpi_%d_%d_%d.bin" % ( fileNameTerrain, mpiX, mpiY, mpiZ ), "rb" )
			
			ny = grid.ny[mpiY]
			nx = grid.nx[mpiX]

			print( "ny = %d, nx = %d" % ( nx, ny ) )
			datax = np.fromfile( XFile, dtype='float32', count = ny * nx )
			datay = np.fromfile( YFile, dtype='float32', count = ny * nx )
			data_ = np.fromfile(  File, dtype='float32', count = ny * nx * PGVSIZE )
			terr  = np.fromfile(  TerFile, dtype='float32', count = ny * nx )
			dataTmp = np.reshape( data_, ( ny * nx, PGVSIZE ) )
			J  = grid.frontNY[mpiY]
			J_ = grid.frontNY[mpiY] + ny
			I  = grid.frontNX[mpiX]
			I_ = grid.frontNX[mpiX] + nx

			if FAST_AXIS == 'Z':
				dataX[I:I_, J:J_] = np.reshape( datax, (nx, ny) )
				dataY[I:I_, J:J_] = np.reshape( datay, (nx, ny) )
				Terr [I:I_, J:J_] = np.reshape( terr,  (nx, ny) )
				data [I:I_, J:J_] = np.reshape( dataTmp[:, PGVh], (nx, ny) )
			else:
				dataX[J:J_, I:I_] = np.reshape( datax, (ny, nx) )
				dataY[J:J_, I:I_] = np.reshape( datay, (ny, nx) )
				Terr [J:J_, I:I_] = np.reshape( terr,  (ny, nx) )
				data [J:J_, I:I_] = np.reshape( dataTmp[:, PGVh], (ny, nx) )
	print(version[ver])
	versionData[version[ver]] = data

PGV = versionData[plotVersion]

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
np.save( "./numpyData/%s_PGVh.npy"%plotVersion, PGVh )
logPGVh = np.log( PGV[nPML:NY - nPML, nPML:NX-nPML] )
np.save( "./numpyData/%s_logPGVh.npy"%plotVersion, logPGVh )

Intensity = Intensity[nPML:NY - nPML, nPML:NX-nPML]
lon = dataX[nPML:NY - nPML, nPML:NX-nPML]
lat = dataY[nPML:NY - nPML, nPML:NX-nPML]

Terr = Terr[nPML:NY - nPML, nPML:NX-nPML]
print(np.max(Intensity))

np.save( "./numpyData/%s_terrain.npy"%plotVersion, Terr)
np.save( "./numpyData/%s_Intensity.npy"%plotVersion, Intensity )
np.save( "./numpyData/%s_lon.npy"%plotVersion, lon )
np.save( "./numpyData/%s_lat.npy"%plotVersion, lat )
    
savemat('./numpyData/%s_terrain.mat'%plotVersion, {'terrain':Terr})
savemat('./numpyData/%s_lon.mat'%plotVersion, {'lon':lon})
savemat('./numpyData/%s_lat.mat'%plotVersion, {'lat':lat})
savemat('./numpyData/%s_PGVh.mat'%plotVersion, {'PGVh':PGVh})
savemat('./numpyData/%s_logPGVh.mat'%plotVersion, {'logPGVh':logPGVh})


plt.figure( 2 )

plt.pcolormesh( lon, lat, Intensity, cmap = "seismic" ) #origin= "lower" )
plt.colorbar()

plt.axis( "image" )
plt.title( '%s Intensity'%plotVersion)

plt.show( )
plt.savefig( "png/"+ plotVersion + "_PGV.png" )

