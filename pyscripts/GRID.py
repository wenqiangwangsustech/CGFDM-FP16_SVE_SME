#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np

class GRID:
	def __init__( self, params, HALO = 3 ):
		resX = 0 	
		resY = 0 
		resZ = 0 	
			
		self.PX = params["PX"] 
		self.PY = params["PY"] 
		self.PZ = params["PZ"] 

		self._NX_ = params["NX"] + 2 * HALO 
		self._NY_ = params["NY"] + 2 * HALO 
		self._NZ_ = params["NZ"] + 2 * HALO 

		self._NX = params["NX"] + HALO 
		self._NY = params["NY"] + HALO 
		self._NZ = params["NZ"] + HALO 

		self.NX = params["NX"] 
		self.NY = params["NY"] 
		self.NZ = params["NZ"] 

		resX = self.NX % self.PX 
		resY = self.NY % self.PY 
		resZ = self.NZ % self.PZ

		self.nx = np.zeros( self.PX, dtype = 'int32' ) + self.NX // self.PX 
		self.ny = np.zeros( self.PY, dtype = 'int32' ) + self.NY // self.PY 
		self.nz = np.zeros( self.PZ, dtype = 'int32' ) + self.NZ // self.PZ 
		
		self.frontNX = np.zeros( self.PX, dtype = 'int32' )
		self.frontNY = np.zeros( self.PY, dtype = 'int32' )
		self.frontNZ = np.zeros( self.PZ, dtype = 'int32' )
		

		self.NonUniform = params["Non-Uniform"]

		for mpiX in range( self.PX ):
			if ( mpiX < resX ):
				self.nx[mpiX] += 1 
				self.frontNX[mpiX] = mpiX * self.nx[mpiX] 
			else:
				self.frontNX[mpiX] = resX * ( self.nx[mpiX] + 1 ) + ( mpiX - resX ) * self.nx[mpiX] 
			if self.NonUniform == 1:
				self.init_heter_grid_X( params, mpiX )

		for mpiY in range( self.PY ):
			if ( mpiY < resY ):
				self.ny[mpiY] += 1 
				self.frontNY[mpiY] = mpiY * self.ny[mpiY] 
			else:
				self.frontNY[mpiY] = resY * ( self.ny[mpiY] + 1 ) + ( mpiY - resY ) * self.ny[mpiY] 
			if self.NonUniform == 1:
				self.init_heter_grid_Y( params, mpiY )

		for mpiZ in range( self.PZ ):
			if ( mpiZ < resZ ):
				self.nz[mpiZ] += 1 
				self.frontNZ[mpiZ] = mpiZ * self.nz[mpiZ] 
			else:
				self.frontNZ[mpiZ] = resZ * ( self.nz[mpiZ] + 1 ) + ( mpiZ - resZ ) * self.nz[mpiZ] 
			if self.NonUniform == 1:
				self.init_heter_grid_Z( params, mpiZ )

		self._frontNX = self.frontNX + HALO 
		self._frontNY = self.frontNY + HALO 
		self._frontNZ = self.frontNZ + HALO 

		self._nx = self.nx + HALO 
		self._ny = self.ny + HALO 
		self._nz = self.nz + HALO 

		self._nx_ = self.nx + 2 * HALO 
		self._ny_ = self.ny + 2 * HALO 
		self._nz_ = self.nz + 2 * HALO 

		self.originalX = params["centerX"] 
		self.originalY = params["centerY"] 

		self._originalX = self.originalX + HALO 
		self._originalY = self.originalY + HALO 

		self.halo = HALO 

		self.DH = params["DH"] 

	def init_heter_grid_X( self, params, mpiX ):
		if self.PX >= 3:
			if (mpiX == 0) or (mpiX == self.PX - 1):
				if mpiX == 0:
					self.nx[mpiX] = params["nxl"]
					self.frontNX[mpiX] = 0;
				if mpiX == self.PX - 1:
					self.nx[mpiX] = params["nxr"]
					self.frontNX[mpiX] = self.NX - self.nx[mpiX];
			else:
				ResNX = params["NX"] - params["nxl"] - params["nxr"]
				self.nx[mpiX] = ResNX // (params["PX"]-2)
				resX = ResNX % (params["PX"]-2)
				if mpiX - 1 < resX:
					self.nx[mpiX] +=1
					self.frontNX[mpiX] = (mpiX-1) * self.nx[mpiX] + params["nxl"]
				else:
					self.frontNX[mpiX] = resX * ( self.nx[mpiX] + 1 ) + ( mpiX - 1 - resX ) * self.nx[mpiX] + params["nxl"]

	def init_heter_grid_Y( self, params, mpiY ):
		if self.PY >= 3:
			if (mpiY == 0) or (mpiY == self.PY - 1):
				if mpiY == 0:
					self.ny[mpiY] = params["nyl"]
					self.frontNY[mpiY] = 0;
				if mpiY == self.PY - 1:
					self.ny[mpiY] = params["nyr"]
					self.frontNY[mpiY] = self.NY - self.ny[mpiY];
			else:
				ResNY = params["NY"] - params["nyl"] - params["nyr"]
				self.ny[mpiY] = ResNY // (params["PY"]-2)
				resY = ResNY % (params["PY"]-2)
				if mpiY - 1 < resY:
					self.ny[mpiY] +=1
					self.frontNY[mpiY] = (mpiY-1) * self.ny[mpiY] + params["nyl"]
				else:
					self.frontNY[mpiY] = resY * ( self.ny[mpiY] + 1 ) + ( mpiY - 1 - resY ) * self.ny[mpiY] + params["nyl"]

	def init_heter_grid_Z( self, params, mpiZ ):
		if self.PZ >= 3:
			if (mpiZ == 0) or (mpiZ == self.PZ - 1):
				if mpiZ == 0:
					self.nz[mpiZ] = params["nzl"]
					self.frontNZ[mpiZ] = 0;
				if mpiZ == self.PZ - 1:
					self.nz[mpiZ] = params["nzr"]
					self.frontNZ[mpiZ] = self.NZ - self.nz[mpiZ];
			else:
				ResNZ = params["NZ"] - params["nzl"] - params["nzr"]
				self.nz[mpiZ] = ResNZ // (params["PZ"]-2)
				resZ = ResNZ % (params["PZ"]-2)
				if mpiZ - 1 < resZ:
					self.nz[mpiZ] +=1
					self.frontNZ[mpiZ] = (mpiZ-1) * self.nz[mpiZ] + params["nzl"]
				else:
					self.frontNZ[mpiZ] = resZ * ( self.nz[mpiZ] + 1 ) + ( mpiZ - 1 - resZ ) * self.nz[mpiZ] + params["nzl"]


def main( ):
	jsonsFile = open( "params.json" )
	params = json.load( jsonsFile )
	grid = GRID( params )
	print( grid )


if __name__ == '__main__':
	main()




