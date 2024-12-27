#!/usr/bin/env python

"""
Created on 2022-02-24
23:59
@author: Wenqiang Wang 
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import struct
import sys
import os
from pyproj import Proj
import json
from mpl_toolkits.mplot3d import Axes3D
from pyscripts.GRID import GRID


stationDic = { }

file = open( 'stationLonLat.txt', 'r' )
LonLatTxt = file.readlines( )
file.close( )

stationNum = len( LonLatTxt )
stationName = []
Lon = np.zeros( stationNum  )
Lat = np.zeros( stationNum  )
for i in range( stationNum ):
	staStr = LonLatTxt[i].split( )
	stationName.append( staStr[0] )
	Lon[i] = float( staStr[1] )
	Lat[i] = float( staStr[2] )

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

DH = params["DH"]
centerX = params["centerX"]
centerY = params["centerY"]

nPML = params['nPML']
XRange1 = nPML + 1
XRange2 = params['NX'] - nPML

YRange1 = nPML + 1
YRange2 = params['NY'] - nPML


NX = params['NX']
NY = params['NY']
NZ = params['NZ']



latC = params["centerLatitude"]
lonC = params["centerLongitude"]
proj = Proj( proj='aeqd', lat_0 = latC, lon_0 = lonC, ellps = "WGS84" )
NZ = params['NZ']

XCoord = np.zeros( stationNum )
YCoord = np.zeros( stationNum )
XCoord, YCoord = proj( Lon, Lat )

X = np.zeros( stationNum )
Y = np.zeros( stationNum )
n = 0


stationDicTmp = { }
for i in range( stationNum ):
	X[i] = round( XCoord[i] / DH ) + centerX
	Y[i] = round( YCoord[i] / DH ) + centerY
	if X[i] >= XRange1 and X[i] < XRange2 and Y[i] >= YRange1 and Y[i] < YRange2:
		I = int( X[i] ); J = int( Y[i] ); K = int( NZ - 1 );
		index = I + J * NX + K * NX * NY
		stationDicTmp[index] = stationName[i]


stationDic = { }


for key in stationDicTmp.keys( ):
    K = key // ( NX * NY )
    J = key % ( NX * NY ) // NX
    I = key % NX
    stationDic[stationDicTmp[key]] = [I, J, K]



stationParam = {

	"point"   : 1,
	"lon_lat" : 0,

	"station(point)" : stationDic
}

print( stationParam  )

stationJsonFile = open( "station.json", 'w' )

json_str = json.dumps( stationParam, indent = 4  )
stationJsonFile.write( json_str )
stationJsonFile.close( )


#json_str = json.dumps( stationParam )
#print( json_str )


