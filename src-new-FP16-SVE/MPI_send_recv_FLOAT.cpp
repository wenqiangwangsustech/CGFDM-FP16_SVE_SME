/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:MPI_send_recv.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-17
*   Discription:
*
================================================================*/
#include "header.h"

typedef void (*FLOAT_PACK_UNPACK_FUNC)( FLOAT * wave, FLOAT * thisSend, int xStartHalo, int _nx_, int _ny_, int _nz_ );


__GLOBAL__
void FLOAT_packX( FLOAT * wave, FLOAT * thisSend, int xStartHalo, int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;
	
	long long num = _nx_ * _ny_ * _nz_;
	CALCULATE3D( i0, j0, k0, 0, HALO, 0, _ny_, 0, _nz_ )
		i = i0 + xStartHalo;
		j = j0;
		k = k0;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, HALO, _ny_, _nz_ ) * WSIZE;
		thisSend[pos + 0] = wave[index + num * 0];
		thisSend[pos + 1] = wave[index + num * 1];
		thisSend[pos + 2] = wave[index + num * 2];
		thisSend[pos + 3] = wave[index + num * 3];
		thisSend[pos + 4] = wave[index + num * 4];
		thisSend[pos + 5] = wave[index + num * 5];
		thisSend[pos + 6] = wave[index + num * 6];
		thisSend[pos + 7] = wave[index + num * 7];
		thisSend[pos + 8] = wave[index + num * 8];
    END_CALCULATE3D( )
}


__GLOBAL__
void FLOAT_unpackX( FLOAT * wave, FLOAT * thisSend, int xStartHalo, int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;
	
	long long num = _nx_ * _ny_ * _nz_;
	CALCULATE3D( i0, j0, k0, 0, HALO, 0, _ny_, 0, _nz_ )
		i = i0 + xStartHalo;
		j = j0;
		k = k0;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, HALO, _ny_, _nz_ ) * WSIZE;
		wave[index + num * 0] = thisSend[pos + 0];
		wave[index + num * 1] = thisSend[pos + 1];
		wave[index + num * 2] = thisSend[pos + 2];
		wave[index + num * 3] = thisSend[pos + 3];
		wave[index + num * 4] = thisSend[pos + 4];
		wave[index + num * 5] = thisSend[pos + 5];
		wave[index + num * 6] = thisSend[pos + 6];
		wave[index + num * 7] = thisSend[pos + 7];
		wave[index + num * 8] = thisSend[pos + 8];
    END_CALCULATE3D( )
}

void FLOAT_PackUnpackX( FLOAT * wave, FLOAT * thisSendRecv, 
	int xStartHalo, int _nx_, int _ny_, int _nz_, FLOAT_PACK_UNPACK_FUNC pack_unpack_func )
{

	int X0 = HALO;

#ifdef GPU_CUDA
	dim3 threads( 4, 8, 16);
	dim3 blocks;
	blocks.x = ( HALO + threads.x - 1 ) / threads.x;
	blocks.y = ( _ny_ + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ + threads.z - 1 ) / threads.z;
	pack_unpack_func<<< blocks, threads >>>
	( wave, thisSendRecv, xStartHalo, _nx_, _ny_, _nz_ );
	CHECK( cudaDeviceSynchronize( ) );
#else
	pack_unpack_func
	( wave, thisSendRecv, xStartHalo, _nx_, _ny_, _nz_ );
#endif
}


__GLOBAL__
void FLOAT_packY( FLOAT * wave, FLOAT * thisSend, int yStartHalo, int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;

	long long num = _nx_ * _ny_ * _nz_;
	CALCULATE3D( i0, j0, k0, 0, _nx_, 0, HALO, 0, _nz_ )
		i = i0;
		j = j0 + yStartHalo;
		k = k0;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, _nx_, HALO, _nz_ ) * WSIZE;
		thisSend[pos + 0] = wave[index + num * 0];
		thisSend[pos + 1] = wave[index + num * 1];
		thisSend[pos + 2] = wave[index + num * 2];
		thisSend[pos + 3] = wave[index + num * 3];
		thisSend[pos + 4] = wave[index + num * 4];
		thisSend[pos + 5] = wave[index + num * 5];
		thisSend[pos + 6] = wave[index + num * 6];
		thisSend[pos + 7] = wave[index + num * 7];
		thisSend[pos + 8] = wave[index + num * 8];
    END_CALCULATE3D( )

}

__GLOBAL__
void FLOAT_unpackY( FLOAT * wave, FLOAT * thisSend, int yStartHalo, int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;

	long long num = _nx_ * _ny_ * _nz_;
	CALCULATE3D( i0, j0, k0, 0, _nx_, 0, HALO, 0, _nz_ )
		i = i0;
		j = j0 + yStartHalo;
		k = k0;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, _nx_, HALO, _nz_ ) * WSIZE;
		wave[index + num * 0] = thisSend[pos + 0];
		wave[index + num * 1] = thisSend[pos + 1];
		wave[index + num * 2] = thisSend[pos + 2];
		wave[index + num * 3] = thisSend[pos + 3];
		wave[index + num * 4] = thisSend[pos + 4];
		wave[index + num * 5] = thisSend[pos + 5];
		wave[index + num * 6] = thisSend[pos + 6];
		wave[index + num * 7] = thisSend[pos + 7];
		wave[index + num * 8] = thisSend[pos + 8];
    END_CALCULATE3D( )

}
void FLOAT_PackUnpackY( FLOAT * wave, FLOAT * thisSendRecv, 
	int yStartHalo, int _nx_, int _ny_, int _nz_, FLOAT_PACK_UNPACK_FUNC pack_unpack_func )
{
	int Y0 = HALO;

#ifdef GPU_CUDA
	dim3 threads( 8, 4, 16);
	dim3 blocks;
	blocks.x = ( _nx_ + threads.x - 1 ) / threads.x;
	blocks.y = ( HALO + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ + threads.z - 1 ) / threads.z;
	pack_unpack_func<<< blocks, threads >>>
	( wave, thisSendRecv, yStartHalo, _nx_, _ny_, _nz_ );
	CHECK( cudaDeviceSynchronize( ) );
#else
	pack_unpack_func
	( wave, thisSendRecv, yStartHalo, _nx_, _ny_, _nz_ );
#endif
}


__GLOBAL__
void FLOAT_packZ( FLOAT * wave, FLOAT * thisSend, int zStartHalo, int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;
	
	long long num = _nx_ * _ny_ * _nz_;
	CALCULATE3D( i0, j0, k0, 0, _nx_, 0, _ny_, 0, HALO )
		i = i0;
		j = j0;
		k = k0 + zStartHalo;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, _nx_, _ny_, HALO ) * WSIZE;
		thisSend[pos + 0] = wave[index + num * 0];
		thisSend[pos + 1] = wave[index + num * 1];
		thisSend[pos + 2] = wave[index + num * 2];
		thisSend[pos + 3] = wave[index + num * 3];
		thisSend[pos + 4] = wave[index + num * 4];
		thisSend[pos + 5] = wave[index + num * 5];
		thisSend[pos + 6] = wave[index + num * 6];
		thisSend[pos + 7] = wave[index + num * 7];
		thisSend[pos + 8] = wave[index + num * 8];
	END_CALCULATE3D( )

}

__GLOBAL__
void FLOAT_unpackZ( FLOAT * wave, FLOAT * thisSend, int zStartHalo, int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;
	
	long long num = _nx_ * _ny_ * _nz_;
	CALCULATE3D( i0, j0, k0, 0, _nx_, 0, _ny_, 0, HALO )
		i = i0;
		j = j0;
		k = k0 + zStartHalo;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, _nx_, _ny_, HALO ) * WSIZE;
		wave[index + num * 0] = thisSend[pos + 0];
		wave[index + num * 1] = thisSend[pos + 1];
		wave[index + num * 2] = thisSend[pos + 2];
		wave[index + num * 3] = thisSend[pos + 3];
		wave[index + num * 4] = thisSend[pos + 4];
		wave[index + num * 5] = thisSend[pos + 5];
		wave[index + num * 6] = thisSend[pos + 6];
		wave[index + num * 7] = thisSend[pos + 7];
		wave[index + num * 8] = thisSend[pos + 8];
	END_CALCULATE3D( )

}
void FLOAT_PackUnpackZ( FLOAT * wave, FLOAT * thisSendRecv, 
	int zStartHalo, int _nx_, int _ny_, int _nz_, FLOAT_PACK_UNPACK_FUNC pack_unpack_func )
{
#ifdef GPU_CUDA
	dim3 threads( 16, 8, 4);
	dim3 blocks;
	blocks.x = ( _nx_ + threads.x - 1 ) / threads.x;
	blocks.y = ( _ny_ + threads.y - 1 ) / threads.y;
	blocks.z = ( HALO + threads.z - 1 ) / threads.z;
	pack_unpack_func<<< blocks, threads >>>
	( wave, thisSendRecv, zStartHalo, _nx_, _ny_, _nz_ );
	CHECK( cudaDeviceSynchronize( ) );
#else
	pack_unpack_func
	( wave, thisSendRecv, zStartHalo, _nx_, _ny_, _nz_ );
#endif
}


void FLOAT_alloc_send_recv( GRID grid, FLOAT ** send, FLOAT ** recv, char XYZ )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	long long num = 0;
	
	switch( XYZ )
	{
		case 'X':
			num = _ny_ * _nz_* HALO * WSIZE;
			break;       
		case 'Y':         
			num = _nx_ * _nz_* HALO * WSIZE;
			break;       
		case 'Z':         
			num = _nx_ * _ny_* HALO * WSIZE;
			break;
	}
		
	FLOAT * pSend = NULL;
	FLOAT * pRecv = NULL;
	long long size = sizeof( FLOAT ) * num;

	CHECK( Malloc( ( void ** )&pSend, size ) );
	CHECK( Memset(  pSend, 0, size ) ); 
	
	CHECK( Malloc( ( void ** )&pRecv, size ) );
	CHECK( Memset(  pRecv, 0, size ) ); 

	*send = pSend;
	*recv = pRecv;

}


void FLOAT_allocSendRecv( GRID grid, MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA_FLOAT * sr)
{
	memset( sr, 0, sizeof( SEND_RECV_DATA_FLOAT ) );
	
	/*if ( mpiNeighbor.X1 > 0 )*/	FLOAT_alloc_send_recv( grid, &( sr->thisXSend1 ), &( sr->thisXRecv1 ), 'X' );
	/*if ( mpiNeighbor.Y1 > 0 )*/	FLOAT_alloc_send_recv( grid, &( sr->thisYSend1 ), &( sr->thisYRecv1 ), 'Y' );
	/*if ( mpiNeighbor.Z1 > 0 )*/	FLOAT_alloc_send_recv( grid, &( sr->thisZSend1 ), &( sr->thisZRecv1 ), 'Z' ); 

	/*if ( mpiNeighbor.X2 > 0 )*/	FLOAT_alloc_send_recv( grid, &( sr->thisXSend2 ), &( sr->thisXRecv2 ), 'X' );
	/*if ( mpiNeighbor.Y2 > 0 )*/	FLOAT_alloc_send_recv( grid, &( sr->thisYSend2 ), &( sr->thisYRecv2 ), 'Y' );
	/*if ( mpiNeighbor.Z2 > 0 )*/	FLOAT_alloc_send_recv( grid, &( sr->thisZSend2 ), &( sr->thisZRecv2 ), 'Z' );
	
}

void FLOAT_freeSendRecv( MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA_FLOAT sr )
{
	/*if ( mpiNeighbor.X1 > 0 )*/	{ Free( sr.thisXSend1);	Free( sr.thisXRecv1 ); };
	/*if ( mpiNeighbor.Y1 > 0 )*/	{ Free( sr.thisYSend1);	Free( sr.thisYRecv1 ); };
	/*if ( mpiNeighbor.Z1 > 0 )*/	{ Free( sr.thisZSend1);	Free( sr.thisZRecv1 ); }; 
                                 
	/*if ( mpiNeighbor.X2 > 0 )*/	{ Free( sr.thisXSend2);	Free( sr.thisXRecv2 ); };
	/*if ( mpiNeighbor.Y2 > 0 )*/	{ Free( sr.thisYSend2);	Free( sr.thisYRecv2 ); };
	/*if ( mpiNeighbor.Z2 > 0 )*/	{ Free( sr.thisZSend2);	Free( sr.thisZRecv2 ); };
	
}


void FLOAT_mpiSendRecv( MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, GRID grid, FLOAT * wave, SEND_RECV_DATA_FLOAT sr, int VARSIZE )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	long long num = 0;

	FLOAT * thisXSend1 = sr.thisXSend1;
	FLOAT * thisXRecv1 = sr.thisXRecv1;
	FLOAT * thisYSend1 = sr.thisYSend1;
	FLOAT * thisYRecv1 = sr.thisYRecv1;
	FLOAT * thisZSend1 = sr.thisZSend1;
	FLOAT * thisZRecv1 = sr.thisZRecv1;
                                 
	FLOAT * thisXSend2 = sr.thisXSend2;
	FLOAT * thisXRecv2 = sr.thisXRecv2;
	FLOAT * thisYSend2 = sr.thisYSend2;
	FLOAT * thisYRecv2 = sr.thisYRecv2;
	FLOAT * thisZSend2 = sr.thisZSend2;
	FLOAT * thisZRecv2 = sr.thisZRecv2;

	int xStartHalo, yStartHalo, zStartHalo;

	MPI_Status stat;

//x direction data exchange
	xStartHalo = nx;
	if ( mpiNeighbor.X2 >= 0 ) FLOAT_PackUnpackX( wave, thisXSend2, xStartHalo, _nx_, _ny_, _nz_, FLOAT_packX );

	num = HALO * _ny_ * _nz_ * VARSIZE * sizeof( FLOAT );
	MPI_Sendrecv( sr.thisXSend2, num, MPI_CHAR, mpiNeighbor.X2, 101,
				  sr.thisXRecv1, num, MPI_CHAR, mpiNeighbor.X1, 101,
				  comm_cart, &stat );

	xStartHalo = 0;
	if ( mpiNeighbor.X1 >= 0 ) FLOAT_PackUnpackX( wave, thisXRecv1, xStartHalo, _nx_, _ny_, _nz_, FLOAT_unpackX );

	
	xStartHalo = HALO;
	if ( mpiNeighbor.X1 >= 0 ) FLOAT_PackUnpackX( wave, thisXSend1, xStartHalo, _nx_, _ny_, _nz_, FLOAT_packX );

	num = HALO * _ny_ * _nz_ * VARSIZE * sizeof( FLOAT );
	MPI_Sendrecv( sr.thisXSend1, num, MPI_CHAR, mpiNeighbor.X1, 102,
				  sr.thisXRecv2, num, MPI_CHAR, mpiNeighbor.X2, 102,
				  comm_cart, &stat );

	xStartHalo = _nx;
	if ( mpiNeighbor.X2 >= 0 ) FLOAT_PackUnpackX( wave, thisXRecv2, xStartHalo, _nx_, _ny_, _nz_, FLOAT_unpackX );

//y direction data exchange
	yStartHalo = ny;
	if ( mpiNeighbor.Y2 >= 0 ) FLOAT_PackUnpackY( wave, thisYSend2, yStartHalo, _nx_, _ny_, _nz_, FLOAT_packY );

	num = HALO * _nx_ * _nz_ * VARSIZE * sizeof( FLOAT );
	MPI_Sendrecv( sr.thisYSend2, num, MPI_CHAR, mpiNeighbor.Y2, 103,
				  sr.thisYRecv1, num, MPI_CHAR, mpiNeighbor.Y1, 103,
				  comm_cart, &stat );

	yStartHalo = 0;
	if ( mpiNeighbor.Y1 >= 0 ) FLOAT_PackUnpackY( wave, thisYRecv1, yStartHalo, _nx_, _ny_, _nz_, FLOAT_unpackY );

	
	yStartHalo = HALO;
	if ( mpiNeighbor.Y1 >= 0 ) FLOAT_PackUnpackY( wave, thisYSend1, yStartHalo, _nx_, _ny_, _nz_, FLOAT_packY );

	num = HALO * _nx_ * _nz_ * VARSIZE * sizeof( FLOAT );
	MPI_Sendrecv( sr.thisYSend1, num, MPI_CHAR, mpiNeighbor.Y1, 104,
				  sr.thisYRecv2, num, MPI_CHAR, mpiNeighbor.Y2, 104,
				  comm_cart, &stat );

	yStartHalo = _ny;
	if ( mpiNeighbor.Y2 >= 0 ) FLOAT_PackUnpackY( wave, thisYRecv2, yStartHalo, _nx_, _ny_, _nz_, FLOAT_unpackY );

//z direction data exchange
	zStartHalo = nz;
	if ( mpiNeighbor.Z2 >= 0 ) FLOAT_PackUnpackZ( wave, thisZSend2, zStartHalo, _nx_, _ny_, _nz_, FLOAT_packZ );

	num = HALO * _nx_ * _ny_ * VARSIZE * sizeof( FLOAT );
	MPI_Sendrecv( sr.thisZSend2, num, MPI_CHAR, mpiNeighbor.Z2, 105,
				  sr.thisZRecv1, num, MPI_CHAR, mpiNeighbor.Z1, 105,
				  comm_cart, &stat );

	zStartHalo = 0;
	if ( mpiNeighbor.Z1 >= 0 ) FLOAT_PackUnpackZ( wave, thisZRecv1, zStartHalo, _nx_, _ny_, _nz_, FLOAT_unpackZ );

	
	zStartHalo = HALO;
	if ( mpiNeighbor.Z1 >= 0 ) FLOAT_PackUnpackZ( wave, thisZSend1, zStartHalo, _nx_, _ny_, _nz_, FLOAT_packZ );

	num = HALO * _nx_ * _ny_ * VARSIZE * sizeof( FLOAT );
	MPI_Sendrecv( sr.thisZSend1, num, MPI_CHAR, mpiNeighbor.Z1, 106,
				  sr.thisZRecv2, num, MPI_CHAR, mpiNeighbor.Z2, 106,
				  comm_cart, &stat );

	zStartHalo = _nz;
	if ( mpiNeighbor.Z2 >= 0 ) FLOAT_PackUnpackZ( wave, thisZRecv2, zStartHalo, _nx_, _ny_, _nz_, FLOAT_unpackZ );

}

