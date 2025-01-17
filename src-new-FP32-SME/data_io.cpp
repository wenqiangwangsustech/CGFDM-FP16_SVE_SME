/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:data_io.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-06
*   Discription:
*
================================================================*/
#include "header.h"

void locateSlice( PARAMS params, GRID grid, SLICE * slice )
{
	int sliceX = params.sliceX;
	int sliceY = params.sliceY;
	int sliceZ = params.sliceZ;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;
	
 	slice->X = sliceX - frontNX + HALO;
 	slice->Y = sliceY - frontNY + HALO;
 	slice->Z = sliceZ - frontNZ + HALO;
	
	//printf( "slice.X = %d, slice.Y = %d, slice.Z = %d\n", slice->X, slice->Y, slice->Z );
	
}

void locateFreeSurfSlice( GRID grid, SLICE * slice )
{
	int sliceX = -1;
	int sliceY = -1;
	int sliceZ = grid.NZ - 1;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;
	
 	slice->X = sliceX - frontNX + HALO;
 	slice->Y = sliceY - frontNY + HALO;
 	slice->Z = sliceZ - frontNZ + HALO;
	
	//printf( "slice.X = %d, slice.Y = %d, slice.Z = %d\n", slice->X, slice->Y, slice->Z );
	
}


void allocDataout( GRID grid, char XYZ, float ** dataout )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	long long num = 0;
	
	switch( XYZ )
	{
		case 'X':
			num = ny * nz;
			break;
		case 'Y':
			num = nx * nz;
			break;
		case 'Z':
			num = nx * ny;
			break;
	}
		
	float * pData = NULL;
	long long size = sizeof( float ) * num;

	CHECK( Malloc( ( void ** )&pData, size ) );
	CHECK( Memset(  pData, 0, size ) ); 
	
	*dataout = pData;
}

void freeDataout( float * dataout  )
{

	Free( dataout  );
}


__GLOBAL__
void pack_iodata_x( void * datain, int VAR, int VARSIZE, float * dataout, int nx, int ny, int nz, int I, int CVS, int FP_TYPE )
{

	//printf( "datain = %p\n", datain  );
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	int _nz_ = nz + 2 * HALO;

	double c = 1.0;

	float * datain1;
	FLOAT * datain2;

	if ( 1 == FP_TYPE ) datain1 = ( float * ) datain;
	if ( 2 == FP_TYPE )	datain2 = ( FLOAT * ) datain;


#ifdef GPU_CUDA
	int j0 = threadIdx.x + blockIdx.x * blockDim.x;
	int k0 = threadIdx.y + blockIdx.y * blockDim.y;
#else
	int j0 = 0;
	int k0 = 0;
#endif
	
	int i = I;
	int j = j0 + HALO;
	int k = k0 + HALO;

	long long index, pos;
	
	if ( CVS == 1 )
		c = 1.0 / Cv;
	if ( CVS == 2 )
		c = 1.0 / Cs;

	CALCULATE2D( j0, k0, 0, ny, 0, nz )
		j = j0 + HALO;
		k = k0 + HALO;
		index = INDEX( i, j, k ) * VARSIZE + VAR;	
		pos = Index2D( j0, k0, ny, nz );
		if ( 1 == FP_TYPE ) dataout[pos] = datain1[index] * c;
		if ( 2 == FP_TYPE )	dataout[pos] = ( float )datain2[index] * c;
		//if ( 1 == FP_TYPE )
		//printf( "1:dataout = %f\n", dataout[pos] );
		//printf( "FP_TYPE = %d\n", FP_TYPE );
	END_CALCULATE2D( )

}


__GLOBAL__
void pack_iodata_y( void * datain, int VAR, int VARSIZE, float * dataout, int nx, int ny, int nz, int J, int CVS, int FP_TYPE )
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	int _nz_ = nz + 2 * HALO;

	double c = 1.0;

	float * datain1;
	FLOAT * datain2;

	if ( 1 == FP_TYPE ) datain1 = ( float * ) datain;
	if ( 2 == FP_TYPE )	datain2 = ( FLOAT * ) datain;

#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int k0 = threadIdx.y + blockIdx.y * blockDim.y;
#else
	int i0 = 0;
	int k0 = 0;
#endif
	
	int i = i0 + HALO;
	int j = J;
	int k = k0 + HALO;

	long long index, pos;
	
	if ( CVS == 1 )
		c = 1.0 / Cv;
	if ( CVS == 2 )
		c = 1.0 / Cs;

	CALCULATE2D( i0, k0, 0, nx, 0, nz )
		i = i0 + HALO;
		k = k0 + HALO;
		index = INDEX( i, j, k ) * VARSIZE + VAR;	
		pos = Index2D( i0, k0, nx, nz );
		if ( 1 == FP_TYPE ) dataout[pos] = datain1[index] * c;
		if ( 2 == FP_TYPE )	dataout[pos] = ( float )datain2[index] * c;
		//printf( "2:dataout = %f\n", dataout[pos] );
	END_CALCULATE2D( )

}


__GLOBAL__
void pack_iodata_z( void * datain, int VAR, int VARSIZE, float * dataout, int nx, int ny, int nz, int K, int CVS, int FP_TYPE )
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	int _nz_ = nz + 2 * HALO;

	double c = 1.0;

	float * datain1;
	FLOAT * datain2;

	if ( 1 == FP_TYPE ) datain1 = ( float * ) datain;
	if ( 2 == FP_TYPE )	datain2 = ( FLOAT * ) datain;

#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
#else
	int i0 = 0;
	int j0 = 0;
#endif
	
	int i = i0 + HALO;
	int j = j0 + HALO;
	int k = K;

	long long index, pos;

	if ( CVS == 1 )
		c = 1.0 / Cv;
	if ( CVS == 2 )
		c = 1.0 / Cs;

	
	CALCULATE2D( i0, j0, 0, nx, 0, ny )
		i = i0 + HALO;
		j = j0 + HALO;
		index = INDEX( i, j, k ) * VARSIZE + VAR;	
		pos = Index2D( i0, j0, nx, ny );
		if ( 1 == FP_TYPE ) dataout[pos] = datain1[index] * c;
		if ( 2 == FP_TYPE )	dataout[pos] = ( float )datain2[index] * c;
		//printf( "3:dataout = %f\n", dataout[pos] );
	END_CALCULATE2D( )


}

void allocDataout_cpu(  GRID grid, char XYZ, float ** dataout )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	long long num = 0;
	
	switch( XYZ )
	{
		case 'X':
			num = ny * nz;
			break;
		case 'Y':
			num = nx * nz;
			break;
		case 'Z':
			num = nx * ny;
			break;
	}
		
	long long size = sizeof( float ) * num;

	*dataout = (float * )malloc( size );
	Memset(  *dataout, 0, size );
	
}

void freeDataout_cpu( float * dataout  )
{
	Free( dataout  );
}



void allocSliceData( GRID grid, SLICE slice, SLICE_DATA * sliceData, SLICE_DATA * sliceDataCpu )
{
	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	if ( slice.X >= HALO && slice.X < _nx )
		allocDataout( grid, 'X', &( sliceData->x ) );
	if ( slice.Y >= HALO && slice.Y < _ny )
		allocDataout( grid, 'Y', &( sliceData->y ) );
	if ( slice.Z >= HALO && slice.Z < _nz )
		allocDataout( grid, 'Z', &( sliceData->z ) );

#ifdef GPU_CUDA
	if ( slice.X >= HALO && slice.X < _nx )
		allocDataout_cpu( grid, 'X', &( sliceDataCpu->x ) );
	if ( slice.Y >= HALO && slice.Y < _ny )
		allocDataout_cpu( grid, 'Y', &( sliceDataCpu->y ) );
	if ( slice.Z >= HALO && slice.Z < _nz )
		allocDataout_cpu( grid, 'Z', &( sliceDataCpu->z ) );
#else
	if ( slice.X >= HALO && slice.X < _nx )
		sliceDataCpu->x = sliceData->x;
	if ( slice.Y >= HALO && slice.Y < _ny )
		sliceDataCpu->y = sliceData->y;
	if ( slice.Z >= HALO && slice.Z < _nz )
		sliceDataCpu->z = sliceData->z;
#endif

}


void freeSliceData( GRID grid, SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu )
{
	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	if ( slice.X >= HALO && slice.X < _nx )
		freeDataout( sliceData.x );
	if ( slice.Y >= HALO && slice.Y < _ny )
		freeDataout( sliceData.y );
	if ( slice.Z >= HALO && slice.Z < _nz )
		freeDataout( sliceData.z );

#ifdef GPU_CUDA
	if ( slice.X >= HALO && slice.X < _nx )
		freeDataout_cpu( sliceDataCpu.x );
	if ( slice.Y >= HALO && slice.Y < _ny )
		freeDataout_cpu( sliceDataCpu.y );
	if ( slice.Z >= HALO && slice.Z < _nz )
		freeDataout_cpu( sliceDataCpu.z );
#endif
}

void data2D_output_bin( GRID grid, SLICE slice, 
						MPI_COORD thisMPICoord, 
						void * datain, int VAR, int VARSIZE, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu,	
						const char * name, int CVS, int FP_TYPE )
{

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;


	if ( slice.X >= HALO && slice.X < _nx )
	{

#ifdef GPU_CUDA
		dim3 threads( 32, 16, 1 );
		dim3 blocks;
		blocks.x = ( ny + threads.x - 1 ) / threads.x;
		blocks.y = ( nz + threads.y - 1 ) / threads.y;
		blocks.z = 1;
		pack_iodata_x
		<<< blocks, threads >>>
		( datain, VAR, VARSIZE, sliceData.x, nx, ny, nz, slice.X, CVS, FP_TYPE );
		long long size = sizeof( float ) * ny * nz;
		CHECK( cudaMemcpy( sliceDataCpu.x, sliceData.x, size, cudaMemcpyDeviceToHost ) );

#else
		pack_iodata_x( datain, VAR, VARSIZE, sliceData.x, nx, ny, nz, slice.X, CVS, FP_TYPE );

		sliceDataCpu.x = sliceData.x;

#endif

		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_X_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 

		fwrite( sliceDataCpu.x, sizeof( float ), ny * nz, fp );
		
		fclose( fp );
	}

	if ( slice.Y >= HALO && slice.Y < _ny )
	{

#ifdef GPU_CUDA
		dim3 threads( 32, 16, 1 );
		dim3 blocks;
		blocks.x = ( nx + threads.x - 1 ) / threads.x;
		blocks.y = ( nz + threads.y - 1 ) / threads.y;
		blocks.z = 1;
		pack_iodata_y
		<<< blocks, threads >>>
		( datain, VAR, VARSIZE, sliceData.y, nx, ny, nz, slice.Y, CVS, FP_TYPE );
		long long size = sizeof( float ) * nx * nz;
		CHECK( cudaMemcpy( sliceDataCpu.y, sliceData.y, size, cudaMemcpyDeviceToHost ) );
#else
		pack_iodata_y( datain, VAR, VARSIZE, sliceData.y, nx, ny, nz, slice.Y, CVS, FP_TYPE );
		sliceDataCpu.y = sliceData.y;
#endif

		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_Y_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 

		fwrite( sliceDataCpu.y, sizeof( float ), nx * nz, fp );
		
		fclose( fp );
	}

	if ( slice.Z >= HALO && slice.Z < _nz )
	{

#ifdef GPU_CUDA
		dim3 threads( 32, 16, 1 );
		dim3 blocks;
		blocks.x = ( nx + threads.x - 1 ) / threads.x;
		blocks.y = ( ny + threads.y - 1 ) / threads.y;
		blocks.z = 1;
		pack_iodata_z
		<<< blocks, threads >>>
		( datain, VAR, VARSIZE, sliceData.z, nx, ny, nz, slice.Z, CVS, FP_TYPE );
		long long size = sizeof( float ) * nx * ny;
		CHECK( cudaMemcpy( sliceDataCpu.z, sliceData.z, size, cudaMemcpyDeviceToHost ) );
#else
		pack_iodata_z( datain, VAR, VARSIZE, sliceData.z, nx, ny, nz, slice.Z, CVS, FP_TYPE );
		sliceDataCpu.z = sliceData.z;
#endif

		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_Z_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 
		
		fwrite( sliceDataCpu.z, sizeof( float ), nx * ny, fp );

		fclose( fp );
	}

}


void data2D_XYZ_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, FLOAT * wave, SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu, char var, int it )
{

	//printf( "wave = %p\n", wave );

	int FP_TYPE = 2;
	long long num = grid._nx_ * grid._ny_ * grid._nz_;
	switch ( var )
	{
		case 'V':
			{

				char VxFileName[128], VyFileName[128], VzFileName[128];

				sprintf( VxFileName, "%s/Vx_%d", params.OUT, it );
				sprintf( VyFileName, "%s/Vy_%d", params.OUT, it );
				sprintf( VzFileName, "%s/Vz_%d", params.OUT, it );

				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 0/*x*/, 1, sliceData, sliceDataCpu, VxFileName, 1 , FP_TYPE );
				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 1/*x*/, 1, sliceData, sliceDataCpu, VyFileName, 1 , FP_TYPE );
				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 2/*x*/, 1, sliceData, sliceDataCpu, VzFileName, 1 , FP_TYPE );
			}	
			break;	
		case 'T':
			{

				char TxxFileName[128], TyyFileName[128], TzzFileName[128];
				char TxyFileName[128], TxzFileName[128], TyzFileName[128];

				sprintf( TxxFileName, "%s/Txx_%d", params.OUT, it );
				sprintf( TyyFileName, "%s/Tyy_%d", params.OUT, it );
				sprintf( TzzFileName, "%s/Tzz_%d", params.OUT, it );
                                             
				sprintf( TxyFileName, "%s/Txy_%d", params.OUT, it );
				sprintf( TxzFileName, "%s/Txz_%d", params.OUT, it );
				sprintf( TyzFileName, "%s/Tyz_%d", params.OUT, it );

				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 3/*x*/, 1, sliceData, sliceDataCpu, TxxFileName, 2 , FP_TYPE );
				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 4/*x*/, 1, sliceData, sliceDataCpu, TyyFileName, 2 , FP_TYPE );
				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 5/*x*/, 1, sliceData, sliceDataCpu, TzzFileName, 2 , FP_TYPE );

				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 6/*x*/, 1, sliceData, sliceDataCpu, TxyFileName, 2 , FP_TYPE );
				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 7/*x*/, 1, sliceData, sliceDataCpu, TxzFileName, 2 , FP_TYPE );
				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 8/*x*/, 1, sliceData, sliceDataCpu, TyzFileName, 2 , FP_TYPE );
			}	
			break;	

		case 'F':
			{
				char VxFileName[128], VyFileName[128], VzFileName[128];
				sprintf( VxFileName, "%s/FreeSurfVx_%d", params.OUT, it );
				sprintf( VyFileName, "%s/FreeSurfVy_%d", params.OUT, it );
				sprintf( VzFileName, "%s/FreeSurfVz_%d", params.OUT, it );

				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 0/*x*/, 1, sliceData, sliceDataCpu, VxFileName, 1 , FP_TYPE );
				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 1/*x*/, 1, sliceData, sliceDataCpu, VyFileName, 1 , FP_TYPE );
				data2D_output_bin( grid, slice, thisMPICoord, wave, num * 2/*x*/, 1, sliceData, sliceDataCpu, VzFileName, 1 , FP_TYPE );
			}

	}

}

void data2D_Model_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * coord, float * medium, SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu )
{
	int FP_TYPE = 1;

	char XName[256], YName[256], ZName[256];

	{
		sprintf( XName, "%s/coordX", params.OUT );
		sprintf( YName, "%s/coordY", params.OUT );
		sprintf( ZName, "%s/coordZ", params.OUT );

		data2D_output_bin( grid, slice, thisMPICoord, coord, 0/*x*/, CSIZE, sliceData, sliceDataCpu, XName, 0 , FP_TYPE );
		data2D_output_bin( grid, slice, thisMPICoord, coord, 1/*y*/, CSIZE, sliceData, sliceDataCpu, YName, 0 , FP_TYPE );
		data2D_output_bin( grid, slice, thisMPICoord, coord, 2/*z*/, CSIZE, sliceData, sliceDataCpu, ZName, 0 , FP_TYPE );

		memset( XName, 0, 256 );
		memset( YName, 0, 256 );
		memset( ZName, 0, 256 );
	}
	{
		sprintf( XName, "%s/Vs",  params.OUT );
		sprintf( YName, "%s/Vp",  params.OUT );
		sprintf( ZName, "%s/rho", params.OUT );

		data2D_output_bin( grid, slice, thisMPICoord, medium, 0/*Vs */, MSIZE, sliceData, sliceDataCpu, XName, 0 , FP_TYPE );
		data2D_output_bin( grid, slice, thisMPICoord, medium, 1/*Vp */, MSIZE, sliceData, sliceDataCpu, YName, 0 , FP_TYPE );
		data2D_output_bin( grid, slice, thisMPICoord, medium, 2/*Rho*/, MSIZE, sliceData, sliceDataCpu, ZName, 0 , FP_TYPE );
	}

}

