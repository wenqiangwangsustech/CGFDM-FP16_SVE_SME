#!/usr/bin/env python
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Wed 27 Jul 2022 08:51:10 PM CST
################################


import numpy as np
import matplotlib.pyplot as plt
import os
import json

exeFile = [
'FP32-CGFDM', 
'FP16-SVE']

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )



for iExe in exeFile:
	params['out'] = 'output_%s' % iExe
	
	print( params )
	paramsJsonFile = open( "params.json_%s" % iExe, 'w' )
	
	json_str = json.dumps( params, indent = 4 )
	paramsJsonFile.write( json_str )
	paramsJsonFile.close( )
