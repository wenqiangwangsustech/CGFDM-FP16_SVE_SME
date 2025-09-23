#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fileX = open( "output/source_coord_X.bin", "rb" )
fileY = open( "output/source_coord_Y.bin", "rb" )
fileZ = open( "output/source_coord_Z.bin", "rb" )


x = np.fromfile( fileX, dtype = "float32" )
y = np.fromfile( fileY, dtype = "float32" )
z = np.fromfile( fileZ, dtype = "float32" )


fileX.close( )
fileY.close( )
fileZ.close( )

dpi = 300
fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1000 # 1km = 1000m
ax = Axes3D( fig )
ax.scatter( x, y, z )
ax.axis( 'image' )

