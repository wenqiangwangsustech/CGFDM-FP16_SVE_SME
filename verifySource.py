#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 17, 2021
13:14
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

pointNum = [ ]
nt = 0
dt = 0.0

index = []
Mxx = []
Myy = []
Mzz = []
Mxy = []
Mxz = []
Myz = []

globalX = []
globalY = []
globalZ = []

for iZ in range( grid.PZ ):
	for iY in range( grid.PY ):
		for iX in range( grid.PX ):
			fileName = "%s/source_mpi_%d_%d_%d.bin" % ( params['out'], iX, iY, iZ )
			#print( fileName  )
			if os.path.exists( fileName ) is False:
				#print( fileName + " doesnot exist" )
				continue
			sourceFile = open( fileName, 'rb' )
			print( fileName )
			p = np.fromfile( sourceFile, dtype='int32', count = 1 )[0]
			pointNum.append( p )
			nt = np.fromfile( sourceFile, dtype='int32', count = 1 )[0]
			dt = np.fromfile( sourceFile, dtype='float32', count = 1 )[0]
			
			_nx_ = grid._nx_[iX]
			_ny_ = grid._ny_[iY]
			_nz_ = grid._nz_[iZ]

			frontNX = grid.frontNX[iX]
			frontNY = grid.frontNY[iY]
			frontNZ = grid.frontNZ[iZ]

			for i in range( p ):
				tmpIndex = np.fromfile( sourceFile, dtype='int64', count = 1 )[0]

				if FAST_AXIS == 'Z':
					tmpX = tmpIndex // ( _nz_ * _ny_ ) + frontNX
					tmpY = tmpIndex % ( _nz_ * _ny_ ) // _nz_ + frontNY
					tmpZ = tmpIndex % ( _nz_ * _ny_ ) % _nz_ + frontNZ

				else:
					tmpZ = tmpIndex // ( _nx_ * _ny_ ) + frontNZ
					tmpY = tmpIndex % ( _nx_ * _ny_ ) // _nx_ + frontNY
					tmpX = tmpIndex % ( _nx_ * _ny_ ) % _nx_ + frontNX

				globalZ.append( tmpZ )
				globalY.append( tmpY )
				globalX.append( tmpX )

				index.append( tmpIndex )
				
				Mxx.append( np.fromfile( sourceFile, dtype='float32', count = nt ) )
				Myy.append( np.fromfile( sourceFile, dtype='float32', count = nt ) )
				Mzz.append( np.fromfile( sourceFile, dtype='float32', count = nt ) )
				Mxy.append( np.fromfile( sourceFile, dtype='float32', count = nt ) )
				Mxz.append( np.fromfile( sourceFile, dtype='float32', count = nt ) )
				Myz.append( np.fromfile( sourceFile, dtype='float32', count = nt ) )

			sourceFile.close( )


#plt.scatter( globalY, globalZ )
#plt.axis( "image" )
#dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
#unit = 1000 # 1km = 1000m
fig = plt.figure( )
ax = Axes3D( fig )
ax.scatter( globalX, globalY, globalZ, c = globalZ )

#plt.savefig( "srcIndex.png" )
#plt.plot( Myy )
plt.savefig( "Mxx.png" )
