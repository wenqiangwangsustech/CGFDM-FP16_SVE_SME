#!/usr/bin/python
'''
Wenqiang Wang@SUSTech 

Date: 2021/12/21

'''


import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import json as js
import struct
from pyproj import Proj
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import CloughTocher2DInterpolator

def gaussFunc( t, t0, a ):
    source = np.exp(-( t - t0 )**2 / ( a ** 2 ) ) / ( np.sqrt( np.pi ) * a )
    return source



NT = 4000
dt = 0.02
strike = 164
dip    = 46
rake   = 122.0 * np.ones([NT])

lat = 35.83
lon = 102.81
Z   = -5000


Mw = 6.2

M0 = ((Mw + 6.06) * 3 * 0.5)**10

t = np.linspace(0, dt * NT, NT)

half_duration = 5.2
t0 = 4 * half_duration
rate = gaussFunc( t, t0, half_duration )
rate = rate / np.max(rate)

area = M0 / (3e10)


jsonsFile = open( "params.json" )
params = js.load( jsonsFile )

sourceFileName = params["sourceDir"] + "/point.bin"
sourceFile = open( sourceFileName, "wb" )

NPTS = 1
value = struct.pack( "i", NPTS )
sourceFile.write( value )

value = struct.pack( "i", NT )
sourceFile.write( value )

value = struct.pack( "f", dt )
sourceFile.write( value )

value = struct.pack( "f",  lon)
sourceFile.write( value )
value = struct.pack( "f",  lat)
sourceFile.write( value )
value = struct.pack( "f",  Z )
sourceFile.write( value )
value = struct.pack( "f",  area)
sourceFile.write( value )
value = struct.pack( "f",  strike)
sourceFile.write( value )
value = struct.pack( "f",  dip)
sourceFile.write( value )
tvalue = struct.pack( "f" * NT,  *( rake) )
sourceFile.write( tvalue )
tvalue = struct.pack( "f" * NT,  *( rate) )
sourceFile.write( tvalue )

sourceFile.close( )


