#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 17, 2021
20:49
'''
import os
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
FAST_AXIS = params["FAST_AXIS"]

index = []

globalX = []
globalY = []
globalZ = []

for iZ in range( grid.PZ ):
	for iY in range( grid.PY ):
		for iX in range( grid.PX ):
			fileName = "%s/srcIndex_%d_%d_%d.bin" % ( params['out'], iX, iY, iZ )
			fileNameP = "%s/source_mpi_%d_%d_%d.bin" % ( params['out'], iX, iY, iZ )

			if os.path.exists( fileNameP ) is False:
				#print( fileName + " doesnot exist" )
				continue
			print( fileName  )
			sourceFile = open( fileName, 'rb' )

			sourceFileP = open( fileNameP, 'rb' )
			p = np.fromfile( sourceFileP, dtype='int32', count = 1 )[0]

			_nx_ = grid._nx_[iX]
			_ny_ = grid._ny_[iY]
			_nz_ = grid._nz_[iZ]

			frontNX = grid.frontNX[iX]
			frontNY = grid.frontNY[iY]
			frontNZ = grid.frontNZ[iZ]


			tmpIndex = np.fromfile( sourceFile, dtype='int64', count = p )
			if FAST_AXIS == 'Z':
				tmpX = tmpIndex // ( _nz_ * _ny_ ) + frontNX
				tmpY = tmpIndex % ( _nz_ * _ny_ ) // _nz_ + frontNY
				tmpZ = tmpIndex % ( _nz_ * _ny_ ) % _nz_ + frontNZ

			else:
				tmpZ = tmpIndex // ( _nx_ * _ny_ ) + frontNZ
				tmpY = tmpIndex % ( _nx_ * _ny_ ) // _nx_ + frontNY
				tmpX = tmpIndex % ( _nx_ * _ny_ ) % _nx_ + frontNX
			
			print(  "p = %d" % p  )
			for i in range( p ):
				globalZ.append( tmpZ[i] )
				globalY.append( tmpY[i] )
				globalX.append( tmpX[i] )

			index.append( tmpIndex )

			sourceFile.close( )
			sourceFileP.close( )


plt.scatter( globalY, globalZ )
#plt.axis( "image" )
dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1000 # 1km = 1000m
#fig = plt.figure( 1 )
#ax = Axes3D( fig )
#ax.scatter( globalX, globalY, globalZ )
plt.savefig( "srcIndex.png" )
