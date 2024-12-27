/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:matrix.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-02
*   Discription:
*
================================================================*/
#include "header.h"

void allocContravariant( GRID grid, float ** contravariant )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 

		
	long long size = sizeof( float ) * num * CONSIZE;

	CHECK( Malloc( ( void ** )contravariant, size ) );
	if ( *contravariant == NULL )
	{
		printf( "can't allocate contravariant memory!\n" );
		MPI_Abort( MPI_COMM_WORLD, 10001 );
	}
	CHECK( Memset( *contravariant, 0, size ) ); 

}

void freeContravariant( float * contravariant )
{	
	Free( contravariant );
}

void allocJac( GRID grid, float ** Jac )
{

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	long long num = _nx_ * _ny_ * _nz_; 

	float * pJac = NULL;
	long long size = sizeof( float ) * num;

	Malloc( ( void ** )&pJac, size );
	Memset(  pJac, 0, size ); 

	*Jac = pJac;

}


void freeJac( float * Jac )
{
	Free( Jac );
}


//When change the fast axises:
/*
 * =============================================
 *             BE careful!!!!!!!!!!!!
 * =============================================
 */

__GLOBAL__
void solve_con_jac( float * con, float * coord, float * Jac, int _nx_, int _ny_, int _nz_, float rDH )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x + HALO; 
	int j = threadIdx.y + blockIdx.y * blockDim.y + HALO;
	int k = threadIdx.z + blockIdx.z * blockDim.z + HALO;
#else
	int i = HALO;
	int j = HALO;
	int k = HALO;
#endif
	long long index = 0;

	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;

	float JacobiInv = 0.0f;

	float x_xi = 0.0f, x_et = 0.0f, x_zt = 0.0f;
	float y_xi = 0.0f, y_et = 0.0f, y_zt = 0.0f;
	float z_xi = 0.0f, z_et = 0.0f, z_zt = 0.0f;

	float xi_x = 0.0f, et_x = 0.0f, zt_x = 0.0f;
	float xi_y = 0.0f, et_y = 0.0f, zt_y = 0.0f;
	float xi_z = 0.0f, et_z = 0.0f, zt_z = 0.0f;

	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )
		index = INDEX( i, j, k );
		x_xi = 0.5f * ( LC( coord/*x*/, 0, CSIZE, 1, xi ) + LC( coord/*x*/, 0, CSIZE, -1, xi ) );	
		x_et = 0.5f * ( LC( coord/*x*/, 0, CSIZE, 1, et ) + LC( coord/*x*/, 0, CSIZE, -1, et ) );	
		x_zt = 0.5f * ( LC( coord/*x*/, 0, CSIZE, 1, zt ) + LC( coord/*x*/, 0, CSIZE, -1, zt ) );

		y_xi = 0.5f * ( LC( coord/*y*/, 1, CSIZE, 1, xi ) + LC( coord/*y*/, 1, CSIZE, -1, xi ) );	
		y_et = 0.5f * ( LC( coord/*y*/, 1, CSIZE, 1, et ) + LC( coord/*y*/, 1, CSIZE, -1, et ) );	
		y_zt = 0.5f * ( LC( coord/*y*/, 1, CSIZE, 1, zt ) + LC( coord/*y*/, 1, CSIZE, -1, zt ) );

		/*
		if ( index == INDEX( _nx/2 - 19, _ny / 2 - 31, _nz - 1 ) )
		{
			printf( "z-1 = %3.10e\n", coord[INDEX_et( i, j, k, -1 )*CSIZE+2]  );
			printf( "z+0 = %3.10e\n", coord[INDEX_et( i, j, k, +0 )*CSIZE+2]  );
			printf( "z+1 = %3.10e\n", coord[INDEX_et( i, j, k, +1 )*CSIZE+2]  );
			printf( "z+2 = %3.10e\n", coord[INDEX_et( i, j, k, +2 )*CSIZE+2]  );
			printf( "z+3 = %3.10e\n", coord[INDEX_et( i, j, k, +3 )*CSIZE+2]  );

			printf( "z-2 = %3.10e\n", coord[INDEX_et( i, j, k, -2 )*CSIZE+2]  );
			printf( "z-3 = %3.10e\n", coord[INDEX_et( i, j, k, -3 )*CSIZE+2]  );

		}
		*/

		
		z_xi = 0.5f * ( LC( coord/*z*/, 2, CSIZE, 1, xi ) + LC( coord/*z*/, 2, CSIZE, -1, xi ) );	
		z_et = 0.5f * ( LC( coord/*z*/, 2, CSIZE, 1, et ) + LC( coord/*z*/, 2, CSIZE, -1, et ) );	
		z_zt = 0.5f * ( LC( coord/*z*/, 2, CSIZE, 1, zt ) + LC( coord/*z*/, 2, CSIZE, -1, zt ) );

		Jac[index] = x_xi * y_et * z_zt + x_et * y_zt * z_xi + x_zt * y_xi * z_et - x_zt * y_et * z_xi - x_et * y_xi * z_zt - x_xi * z_et * y_zt;
		JacobiInv = 1.0f / Jac[index];


		xi_x = ( y_et * z_zt - y_zt * z_et ) * JacobiInv;	
		xi_y = ( x_zt * z_et - x_et * z_zt ) * JacobiInv;	
		xi_z = ( x_et * y_zt - x_zt * y_et ) * JacobiInv;

		et_x = ( y_zt * z_xi - y_xi * z_zt ) * JacobiInv;	
		et_y = ( x_xi * z_zt - x_zt * z_xi ) * JacobiInv;	
		et_z = ( x_zt * y_xi - x_xi * y_zt ) * JacobiInv;

		zt_x = ( y_xi * z_et - y_et * z_xi ) * JacobiInv;	
		zt_y = ( x_et * z_xi - x_xi * z_et ) * JacobiInv;	
		zt_z = ( x_xi * y_et - x_et * y_xi ) * JacobiInv;

		con[index*CONSIZE+0] = xi_x; con[index*CONSIZE+1] = xi_y; con[index*CONSIZE+2] = xi_z;
		con[index*CONSIZE+3] = et_x; con[index*CONSIZE+4] = et_y; con[index*CONSIZE+5] = et_z;
		con[index*CONSIZE+6] = zt_x; con[index*CONSIZE+7] = zt_y; con[index*CONSIZE+8] = zt_z;
		//if ( index == INDEX( _nx/2 - 10, _ny / 2 - 10, _nz - 1 ) )
		/*
		if ( index == INDEX( _nx/2 - 19, _ny / 2 - 31, _nz - 1 ) )
		{
			printf( "x = %3.10e\n", coord[index*CSIZE+0]  );
			printf( "y = %3.10e\n", coord[index*CSIZE+1]  );
			printf( "z = %3.10e\n", coord[index*CSIZE+2]  );

			printf( "Jac = %3.10e\n", JacobiInv );

			printf("x_xi = %3.10e\n", x_xi );
			printf("y_xi = %3.10e\n", y_xi );
			printf("z_xi = %3.10e\n", z_xi );
			printf("x_et = %3.10e\n", x_et );
			printf("y_et = %3.10e\n", y_et );
			printf("z_et = %3.10e\n", z_et );
			printf("x_zt = %3.10e\n", x_zt );
			printf("y_zt = %3.10e\n", y_zt );
			printf("z_zt = %3.10e\n", z_zt );
		}
		*/

	END_CALCULATE3D( )
}

//When change the fast axises:
/*
 * =============================================
 *             BE careful!!!!!!!!!!!!
 * =============================================
 */
void allocMat_rDZ( GRID grid, Mat_rDZ * mat_rDZ )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	
	long long num = _nx_ * _ny_; 
		
	float * pSurf = NULL;
	long long size = sizeof( float ) * num * 2 * MATSIZE;

	Malloc( ( void ** )&pSurf, size );
	if ( pSurf == NULL )
	{
		printf( "can't allocate Mat_rDZ memory!\n" );
		MPI_Abort( MPI_COMM_WORLD, 10001 );
	}
	Memset(  pSurf, 0, size ); 
	
	mat_rDZ->_rDZ_DX = pSurf + 0 * num; 
	mat_rDZ->_rDZ_DY = pSurf + MATSIZE * num; 

}

void freeMat_rDZ( Mat_rDZ mat_rDZ )
{
	Free( mat_rDZ._rDZ_DX );
}


__GLOBAL__
void solve_coordinate_on_free_surface( 
	float * con, float * coord, float * Jac,
	float * medium, Mat_rDZ mat_rDZ, 
	int _nx_, int _ny_, int _nz_ ) 
{
	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;
	

#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x + HALO; 
	int j = threadIdx.y + blockIdx.y * blockDim.y + HALO;
	int k = threadIdx.z + blockIdx.z * blockDim.z + _nz - 1;
#else
	int i = HALO;
	int j = HALO;
	int k = _nz - 1;
#endif

	long long index = 0;
	long long pos = 0;

	float DZ[9];
	float DX[9];
	float DY[9];
	float rDZ[9]; //the inverse matrix of DZ
	float DZ_det;
	float chi = 0.0f;
	float lambda = 0.0f;
	float mu = 0.0f;
	float lam_mu = 0.0f;
	
	float * _rDZ_DX = mat_rDZ._rDZ_DX;
	float * _rDZ_DY = mat_rDZ._rDZ_DY;



	int indexOnSurf = 0;

	float xi_x = 0.0f, et_x = 0.0f, zt_x = 0.0f;
	float xi_y = 0.0f, et_y = 0.0f, zt_y = 0.0f;
	float xi_z = 0.0f, et_z = 0.0f, zt_z = 0.0f;	


	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, _nz - 1, _nz_ )
		index = INDEX( i, j, k );

		if ( k == _nz - 1 )
		{	
			xi_x = con[index*CONSIZE+0];  xi_y = con[index*CONSIZE+1];  xi_z = con[index*CONSIZE+2];
			et_x = con[index*CONSIZE+3];  et_y = con[index*CONSIZE+4];  et_z = con[index*CONSIZE+5];
			zt_x = con[index*CONSIZE+6];  zt_y = con[index*CONSIZE+7];  zt_z = con[index*CONSIZE+8];

			indexOnSurf = Index2D( i, j, _nx_, _ny_ );
			mu 		= medium[index*MSIZE+0];
			lambda 	= medium[index*MSIZE+1];
			chi = 2.0f * lambda + mu;//希腊字母
			lam_mu = lambda + mu;
			//printf( "mu = %f, lambda = %f\n", mu, lambda );
			/****************
			---				  ---
			| DZ[0] DZ[1] DZ[2] |
			| DZ[3] DZ[4] DZ[5] |
			| DZ[6] DZ[7] DZ[8] |
			---				  ---
			*****************/
			DZ[0] = chi * zt_x * zt_x + mu * ( zt_y * zt_y + zt_z * zt_z );	
			DZ[1] = lam_mu * zt_x * zt_y; 									
			DZ[2] = lam_mu * zt_x * zt_z;
			
			DZ[3] = DZ[1];													
			DZ[4] = chi * zt_y * zt_y + mu * ( zt_x * zt_x + zt_z * zt_z );	
			DZ[5] = lam_mu * zt_y * zt_z;
			
			DZ[6] = DZ[2];													
			DZ[7] = DZ[5];													
			DZ[8] = chi * zt_z * zt_z + mu * ( zt_x * zt_x + zt_y * zt_y );

			DZ_det 	= DZ[0] * DZ[4] * DZ[8] 
					+ DZ[1] * DZ[5] * DZ[6] 
					+ DZ[2] * DZ[7] * DZ[3] 
					- DZ[2] * DZ[4] * DZ[6] 
					- DZ[1] * DZ[3] * DZ[8] 
					- DZ[0] * DZ[7] * DZ[5];

			DX[0] = chi * zt_x * xi_x + mu * ( zt_y * xi_y + zt_z * xi_z );		
			DX[1] = lambda * zt_x * xi_y + mu * zt_y * xi_x;					
			DX[2] = lambda * zt_x * xi_z + mu * zt_z * xi_x;
			
			DX[3] = lambda * zt_y * xi_x + mu * zt_x * xi_y;					
			DX[4] = chi * zt_y * xi_y + mu * ( zt_x * xi_x + zt_z * xi_z );		
			DX[5] = lambda * zt_y * xi_z + mu * zt_z * xi_y;
			
			DX[6] = lambda * zt_z * xi_x + mu * zt_x * xi_z;					
			DX[7] = lambda * zt_z * xi_y + mu * zt_y * xi_z;					
			DX[8] = chi * zt_z * xi_z + mu * ( zt_x * xi_x + zt_y * xi_y );

			DY[0] = chi * zt_x * et_x + mu * ( zt_y * et_y + zt_z * et_z );		
			DY[1] = lambda * zt_x * et_y + mu * zt_y * et_x;					
			DY[2] = lambda * zt_x * et_z + mu * zt_z * et_x;
			
			DY[3] = lambda * zt_y * et_x + mu * zt_x * et_y;					
			DY[4] = chi * zt_y * et_y + mu * ( zt_x * et_x + zt_z * et_z );		
			DY[5] = lambda * zt_y * et_z + mu * zt_z * et_y;
			
			DY[6] = lambda * zt_z * et_x + mu * zt_x * et_z;					
			DY[7] = lambda * zt_z * et_y + mu * zt_y * et_z;					
			DY[8] = chi * zt_z * et_z + mu * ( zt_x * et_x + zt_y * et_y );

			rDZ[0] = (   DZ[4] * DZ[8] - DZ[5] * DZ[7] ) / DZ_det; 
			rDZ[1] = ( - DZ[3] * DZ[8] + DZ[5] * DZ[6] ) / DZ_det; 
			rDZ[2] = (   DZ[3] * DZ[7] - DZ[4] * DZ[6] ) / DZ_det; 
			rDZ[3] = ( - DZ[1] * DZ[8] + DZ[2] * DZ[7] ) / DZ_det; 
			rDZ[4] = (   DZ[0] * DZ[8] - DZ[2] * DZ[6] ) / DZ_det; 
			rDZ[5] = ( - DZ[0] * DZ[7] + DZ[1] * DZ[6] ) / DZ_det; 
			rDZ[6] = (   DZ[1] * DZ[5] - DZ[2] * DZ[4] ) / DZ_det; 
			rDZ[7] = ( - DZ[0] * DZ[5] + DZ[2] * DZ[3] ) / DZ_det; 
			rDZ[8] = (   DZ[0] * DZ[4] - DZ[1] * DZ[3] ) / DZ_det;
			
			_rDZ_DX[indexOnSurf*MATSIZE+0] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DX[0], DX[3], DX[6] );//M11
			_rDZ_DX[indexOnSurf*MATSIZE+1] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DX[1], DX[4], DX[7] );//M12
			_rDZ_DX[indexOnSurf*MATSIZE+2] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DX[2], DX[5], DX[8] );//M13
                               //   
			_rDZ_DX[indexOnSurf*MATSIZE+3] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DX[0], DX[3], DX[6] );//M21
			_rDZ_DX[indexOnSurf*MATSIZE+4] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DX[1], DX[4], DX[7] );//M22
			_rDZ_DX[indexOnSurf*MATSIZE+5] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DX[2], DX[5], DX[8] );//M23
                               //   
			_rDZ_DX[indexOnSurf*MATSIZE+6] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DX[0], DX[3], DX[6] );//M31
			_rDZ_DX[indexOnSurf*MATSIZE+7] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DX[1], DX[4], DX[7] );//M32
			_rDZ_DX[indexOnSurf*MATSIZE+8] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DX[2], DX[5], DX[8] );//M33
                               //   
			_rDZ_DY[indexOnSurf*MATSIZE+0] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DY[0], DY[3], DY[6] );//M11
			_rDZ_DY[indexOnSurf*MATSIZE+1] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DY[1], DY[4], DY[7] );//M12
			_rDZ_DY[indexOnSurf*MATSIZE+2] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DY[2], DY[5], DY[8] );//M13
                               //   
			_rDZ_DY[indexOnSurf*MATSIZE+3] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DY[0], DY[3], DY[6] );//M21
			_rDZ_DY[indexOnSurf*MATSIZE+4] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DY[1], DY[4], DY[7] );//M22
			_rDZ_DY[indexOnSurf*MATSIZE+5] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DY[2], DY[5], DY[8] );//M23
			                   //   
			_rDZ_DY[indexOnSurf*MATSIZE+6] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DY[0], DY[3], DY[6] );//M31
			_rDZ_DY[indexOnSurf*MATSIZE+7] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DY[1], DY[4], DY[7] );//M32
			_rDZ_DY[indexOnSurf*MATSIZE+8] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DY[2], DY[5], DY[8] );//M33

			/*
			if ( index == INDEX( _nx/2 - 19, _ny / 2 - 31, _nz - 1 ) )
			//if ( index == INDEX( _nx/2 - 10, _ny / 2 - 10, _nz - 1 ) )
			{
				printf("xi_x = %3.10e\n", xi_x );
				printf("xi_y = %3.10e\n", xi_y );
				printf("xi_z = %3.10e\n", xi_z );
				printf("et_x = %3.10e\n", et_x );
				printf("et_y = %3.10e\n", et_y );
				printf("et_z = %3.10e\n", et_z );
				printf("zt_x = %3.10e\n", zt_x );
				printf("zt_y = %3.10e\n", zt_y );
				printf("zt_z = %3.10e\n", zt_z );
				                              
				printf("_rDZ_DX0 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+0]);
				printf("_rDZ_DX1 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+1]);
				printf("_rDZ_DX2 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+2]);
				                                                //   
				printf("_rDZ_DX3 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+3]);
				printf("_rDZ_DX4 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+4]);
				printf("_rDZ_DX5 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+5]);
				                                                //   
				printf("_rDZ_DX6 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+6]);
				printf("_rDZ_DX7 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+7]);
				printf("_rDZ_DX8 = %3.10e\n", _rDZ_DX[indexOnSurf*MATSIZE+8]);
				                                                 //   
				printf("_rDZ_DY0 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+0]);
				printf("_rDZ_DY1 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+1]);
				printf("_rDZ_DY2 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+2]);
				                                                 //   
				printf("_rDZ_DY3 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+3]);
				printf("_rDZ_DY4 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+4]);
				printf("_rDZ_DY5 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+5]);
				                                                 //   
				printf("_rDZ_DY6 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+6]);
				printf("_rDZ_DY7 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+7]);
				printf("_rDZ_DY8 = %3.10e\n", _rDZ_DY[indexOnSurf*MATSIZE+8]);
			}
			*/
		}
		else
		{
			pos = INDEX( i, j, 2 * ( _nz - 1 ) - k );
			Jac[index] = Jac[pos];
			
			con[index*CONSIZE+0] = con[pos*CONSIZE+0]; 
			con[index*CONSIZE+1] = con[pos*CONSIZE+1];
			con[index*CONSIZE+2] = con[pos*CONSIZE+2];
			con[index*CONSIZE+3] = con[pos*CONSIZE+3]; 
			con[index*CONSIZE+4] = con[pos*CONSIZE+4]; 
			con[index*CONSIZE+5] = con[pos*CONSIZE+5];
			con[index*CONSIZE+6] = con[pos*CONSIZE+6]; 
			con[index*CONSIZE+7] = con[pos*CONSIZE+7]; 
			con[index*CONSIZE+8] = con[pos*CONSIZE+8];

		}

	END_CALCULATE3D( )

}


void solveContravariantJac( MPI_Comm comm_cart, GRID grid, float * con, float * coord, float * Jac, float * medium, Mat_rDZ mat_rDZ )
{
	float DH = grid.DH;
	float rDH = 1.0 / DH;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

#ifdef GPU_CUDA
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	dim3 threads( 32, 4, 4);
	dim3 blocks;
	blocks.x = ( nx + threads.x - 1 ) / threads.x;
	blocks.y = ( ny + threads.y - 1 ) / threads.y;
	blocks.z = ( nz + threads.z - 1 ) / threads.z;

	solve_con_jac
	<<< blocks, threads >>>
	( con, coord, Jac, _nx_, _ny_, _nz_, rDH );
#else
	solve_con_jac( con, coord, Jac, _nx_, _ny_, _nz_, rDH );
#endif
	

#ifdef FREE_SURFACE
	
#ifdef GPU_CUDA
	dim3 threadSurf( 32, 16, 1);
	dim3 blockSurf;
	blockSurf.x = ( nx + threadSurf.x - 1 ) / threadSurf.x;
	blockSurf.y = ( ny + threadSurf.y - 1 ) / threadSurf.y;
	blockSurf.z = HALO + 1;

	solve_coordinate_on_free_surface 
	<<< blockSurf, threadSurf >>>
	(	con, coord, Jac,
		medium, mat_rDZ, 
		_nx_, _ny_, _nz_ ); 
#else
	solve_coordinate_on_free_surface 
	(	con, coord, Jac,
		medium, mat_rDZ, 
		_nx_, _ny_, _nz_ ); 
	
#endif //GPU_CUDA

#endif //FREE_SURFACE

}
