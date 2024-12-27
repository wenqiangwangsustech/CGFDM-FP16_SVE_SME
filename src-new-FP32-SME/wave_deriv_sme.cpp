/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:wave.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-03
*   Discription:
*
================================================================*/

#include "header.h"
#include <arm_sve.h>
#include <arm_sme.h>


#ifdef PML
#define TIMES_PML_BETA_X * pml_beta_x
#define TIMES_PML_BETA_Y * pml_beta_y
#define TIMES_PML_BETA_Z * pml_beta_z
#else
#define TIMES_PML_BETA_X 
#define TIMES_PML_BETA_Y 
#define TIMES_PML_BETA_Z 
#endif


void allocWave( GRID grid, WAVE * wave )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 
		
	FLOAT * pWave = NULL;
	long long size = sizeof( FLOAT ) * num * WSIZE * 4;

	CHECK( Malloc( ( void ** )&pWave, size ) );

	if ( pWave== NULL )
	{
		printf( "can't allocate Wave memory!\n" );
		MPI_Abort( MPI_COMM_WORLD, 10001 );
	}
	CHECK( Memset( pWave, 0, size ) ); 

	wave->h_W = pWave + 0 * WSIZE * num;
	wave->  W = pWave + 1 * WSIZE * num;
	wave->t_W = pWave + 2 * WSIZE * num;
	wave->m_W = pWave + 3 * WSIZE * num;
}


void freeWave( WAVE wave )
{
	Free( wave.h_W );
}



__GLOBAL__		
__arm_new("za")
void wave_deriv( FLOAT * h_W, FLOAT * W, FLOAT * CJM, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int _nx_, int _ny_, int _nz_, float rDH, int FB1, int FB2, int FB3, float DT )	
__arm_streaming 
{																	

	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;

#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x + HALO; 
	int j = threadIdx.y + blockIdx.y * blockDim.y + HALO;
	int k = threadIdx.z + blockIdx.z * blockDim.z + HALO;
#else
	int i = HALO;
	int j = HALO;
	int k = HALO;
#endif
	long long index;

#ifdef PML
	svfloat32_t pml_beta_x;
	float pml_beta_y;
	float pml_beta_z;
#endif
																	
	svfloat32_t mu;												
	svfloat32_t lambda;											
	svfloat32_t buoyancy;											
																	
	svfloat32_t xi_x;   svfloat32_t xi_y; 	svfloat32_t xi_z;		
	svfloat32_t et_x;   svfloat32_t et_y; 	svfloat32_t et_z;		
	svfloat32_t zt_x;   svfloat32_t zt_y; 	svfloat32_t zt_z;		
																	
	svfloat32_t Txx_xi; svfloat32_t Tyy_xi; svfloat32_t Txy_xi; 	
	svfloat32_t Txx_et; svfloat32_t Tyy_et; svfloat32_t Txy_et; 	
	svfloat32_t Txx_zt; svfloat32_t Tyy_zt; svfloat32_t Txy_zt; 	
	svfloat32_t Txz_xi; svfloat32_t Tyz_xi; svfloat32_t Tzz_xi;	
	svfloat32_t Txz_et; svfloat32_t Tyz_et; svfloat32_t Tzz_et;	
	svfloat32_t Txz_zt; svfloat32_t Tyz_zt; svfloat32_t Tzz_zt;	
	svfloat32_t Vx_xi ; svfloat32_t Vx_et ; svfloat32_t Vx_zt ;	
	svfloat32_t Vy_xi ; svfloat32_t Vy_et ; svfloat32_t Vy_zt ;	
	svfloat32_t Vz_xi ; svfloat32_t Vz_et ; svfloat32_t Vz_zt ;

	svfloat32_t T0;
	svfloat32_t T1;
	svfloat32_t T2;
	svfloat32_t T3;
	svfloat32_t T4;
	svfloat32_t T5;
	svfloat32_t T6;
	svfloat32_t T7;
	svfloat32_t T8;

	int len = svcntw();

	svfloat32_t Z_1, Z0, Z1, Z2, Z3, Z4;
	svfloat32_t TMP;

	svbool_t pg;
	svbool_t pgIdx;

	svfloat32_t coef;
	svfloat32_t lam2mu;
	



	long long num = _nx_ * _ny_ * _nz_; 
	CALCULATE3D_SVE( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )  
		pg = svwhilelt_b32( i, _nx );
		svzero_za();
		index = INDEX( i, j, k ); 									

		L_SVE( Vx_xi, W, 0, FB1, xi ); 				
		L_SVE( Vy_xi, W, 1, FB1, xi ); 				
		L_SVE( Vz_xi, W, 2, FB1, xi ); 				
		L_SVE(Txx_xi, W, 3, FB1, xi ); 				
		L_SVE(Tyy_xi, W, 4, FB1, xi ); 				
		L_SVE(Tzz_xi, W, 5, FB1, xi ); 				
		L_SVE(Txy_xi, W, 6, FB1, xi ); 				
		L_SVE(Txz_xi, W, 7, FB1, xi ); 				
		L_SVE(Tyz_xi, W, 8, FB1, xi ); 				

		L_SVE( Vx_et, W, 0, FB2, et );				
		L_SVE( Vy_et, W, 1, FB2, et );				
		L_SVE( Vz_et, W, 2, FB2, et );				
		L_SVE(Txx_et, W, 3, FB2, et );				
		L_SVE(Tyy_et, W, 4, FB2, et );				
		L_SVE(Tzz_et, W, 5, FB2, et );				
		L_SVE(Txy_et, W, 6, FB2, et );				
		L_SVE(Txz_et, W, 7, FB2, et );				
		L_SVE(Tyz_et, W, 8, FB2, et );				

   		L_SVE( Vx_zt, W, 0, FB3, zt );				
   		L_SVE( Vy_zt, W, 1, FB3, zt );				
   		L_SVE( Vz_zt, W, 2, FB3, zt );				
  		L_SVE(Txx_zt, W, 3, FB3, zt );				
  		L_SVE(Tyy_zt, W, 4, FB3, zt );				
  		L_SVE(Tzz_zt, W, 5, FB3, zt );				
  		L_SVE(Txy_zt, W, 6, FB3, zt );				
  		L_SVE(Txz_zt, W, 7, FB3, zt );				
  		L_SVE(Tyz_zt, W, 8, FB3, zt );

#ifdef PML
		pml_beta_x = svld1( pg, &(pml_beta.x[i]) );
		pml_beta_y = pml_beta.y[j];
		pml_beta_z = pml_beta.z[k];
#endif


#ifdef PML
		 Vx_xi = svmul_m( pg,  Vx_xi, pml_beta_x );
		 Vy_xi = svmul_m( pg,  Vy_xi, pml_beta_x );
		 Vz_xi = svmul_m( pg,  Vz_xi, pml_beta_x );
		 Vx_et = svmul_m( pg,  Vx_et, pml_beta_y );
		 Vy_et = svmul_m( pg,  Vy_et, pml_beta_y );
		 Vz_et = svmul_m( pg,  Vz_et, pml_beta_y );
		 Vx_zt = svmul_m( pg,  Vx_zt, pml_beta_z );
		 Vy_zt = svmul_m( pg,  Vy_zt, pml_beta_z );
		 Vz_zt = svmul_m( pg,  Vz_zt, pml_beta_z );

		Txx_xi = svmul_m( pg, Txx_xi, pml_beta_x );
		Tyy_xi = svmul_m( pg, Tyy_xi, pml_beta_x );
		Tzz_xi = svmul_m( pg, Tzz_xi, pml_beta_x );
		Txy_xi = svmul_m( pg, Txy_xi, pml_beta_x );
		Txz_xi = svmul_m( pg, Txz_xi, pml_beta_x );
		Tyz_xi = svmul_m( pg, Tyz_xi, pml_beta_x );

		Txx_et = svmul_m( pg, Txx_et, pml_beta_y );
		Tyy_et = svmul_m( pg, Tyy_et, pml_beta_y );
		Tzz_et = svmul_m( pg, Tzz_et, pml_beta_y );
		Txy_et = svmul_m( pg, Txy_et, pml_beta_y );
		Txz_et = svmul_m( pg, Txz_et, pml_beta_y );
		Tyz_et = svmul_m( pg, Tyz_et, pml_beta_y );

		Txx_zt = svmul_m( pg, Txx_zt, pml_beta_z );
		Tyy_zt = svmul_m( pg, Tyy_zt, pml_beta_z );
		Tzz_zt = svmul_m( pg, Tzz_zt, pml_beta_z );
		Txy_zt = svmul_m( pg, Txy_zt, pml_beta_z );
		Txz_zt = svmul_m( pg, Txz_zt, pml_beta_z );
		Tyz_zt = svmul_m( pg, Tyz_zt, pml_beta_z );
#endif
		xi_x = svld1( pg, &(CJM[index+num*0])); 
		xi_y = svld1( pg, &(CJM[index+num*1])); 
		xi_z = svld1( pg, &(CJM[index+num*2]));	
		et_x = svld1( pg, &(CJM[index+num*3])); 
		et_y = svld1( pg, &(CJM[index+num*4])); 
		et_z = svld1( pg, &(CJM[index+num*5]));	
		zt_x = svld1( pg, &(CJM[index+num*6])); 
		zt_y = svld1( pg, &(CJM[index+num*7])); 
		zt_z = svld1( pg, &(CJM[index+num*8]));	

		mu 		 = svld1( pg, &(CJM[index+num*10]));
		lambda 	 = svld1( pg, &(CJM[index+num*11]));	
		buoyancy = svld1( pg, &(CJM[index+num*12]));
		buoyancy = svmul_m( pg, buoyancy, Crho);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 0); //dVx/dt
		T0 = svmul_m(pg, Txx_xi, xi_x); 
		T1 = svmul_m(pg, Txy_xi, xi_y); 
		T2 = svmul_m(pg, Txz_xi, xi_z); 
		T3 = svmul_m(pg, Txx_et, et_x); 
		T4 = svmul_m(pg, Txy_et, et_y); 
		T5 = svmul_m(pg, Txz_et, et_z); 
		T6 = svmul_m(pg, Txx_zt, zt_x); 
		T7 = svmul_m(pg, Txy_zt, zt_y); 
		T8 = svmul_m(pg, Txz_zt, zt_z); 
		
		svmopa_za32_m(3, pgIdx, pg, buoyancy, T0);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T1);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T2);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T3);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T4);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T5);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T6);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T7);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T8);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 1); //dVy/dt
		T0 = svmul_m(pg, Txy_xi, xi_x);
		T1 = svmul_m(pg, Tyy_xi, xi_y);
		T2 = svmul_m(pg, Tyz_xi, xi_z);
		T3 = svmul_m(pg, Txy_et, et_x);
		T4 = svmul_m(pg, Tyy_et, et_y);
		T5 = svmul_m(pg, Tyz_et, et_z);
		T6 = svmul_m(pg, Txy_zt, zt_x);
		T7 = svmul_m(pg, Tyy_zt, zt_y);
		T8 = svmul_m(pg, Tyz_zt, zt_z);

		svmopa_za32_m(3, pgIdx, pg, buoyancy, T0);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T1);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T2);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T3);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T4);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T5);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T6);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T7);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T8);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 2); //dVz/dt  
		T0 = svmul_m(pg, Txz_xi, xi_x); 
		T1 = svmul_m(pg, Tyz_xi, xi_y); 
		T2 = svmul_m(pg, Tzz_xi, xi_z); 
		T3 = svmul_m(pg, Txz_et, et_x); 
		T4 = svmul_m(pg, Tyz_et, et_y); 
		T5 = svmul_m(pg, Tzz_et, et_z); 
		T6 = svmul_m(pg, Txz_zt, zt_x); 
		T7 = svmul_m(pg, Tyz_zt, zt_y); 
		T8 = svmul_m(pg, Tzz_zt, zt_z); 

		svmopa_za32_m(3, pgIdx, pg, buoyancy, T0);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T1);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T2);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T3);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T4);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T5);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T6);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T7);
        svmopa_za32_m(3, pgIdx, pg, buoyancy, T8);


		//Tile 0:xi 1:et 2:zt 3:d /dt
		lam2mu = svmul_m(pg, mu, 2.0f);
		lam2mu = svadd_m(pg, lam2mu, lambda);
		svfloat32_t T0 = svmul_m(pg, Vx_xi, xi_x);
		svfloat32_t T1 = svmul_m(pg, Vx_et, et_x);
		svfloat32_t T2 = svmul_m(pg, Vx_zt, zt_x);
		svfloat32_t T3 = svmul_m(pg, Vy_xi, xi_y);
		svfloat32_t T4 = svmul_m(pg, Vy_et, et_y);
		svfloat32_t T5 = svmul_m(pg, Vy_zt, zt_y);
		svfloat32_t T6 = svmul_m(pg, Vz_xi, xi_z);
		svfloat32_t T7 = svmul_m(pg, Vz_et, et_z);
		svfloat32_t T8 = svmul_m(pg, Vz_zt, zt_z);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 3); //dTxx/dt  
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T0);
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T1);
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T2);
		svmopa_za32_m(3, pgIdx, pg, lambda, T3);
		svmopa_za32_m(3, pgIdx, pg, lambda, T4);
		svmopa_za32_m(3, pgIdx, pg, lambda, T5);
		svmopa_za32_m(3, pgIdx, pg, lambda, T6);
		svmopa_za32_m(3, pgIdx, pg, lambda, T7);
		svmopa_za32_m(3, pgIdx, pg, lambda, T8);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 4); //dTyy/dt
		svmopa_za32_m(3, pgIdx, pg, lambda, T0);
		svmopa_za32_m(3, pgIdx, pg, lambda, T1);
		svmopa_za32_m(3, pgIdx, pg, lambda, T2);
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T3);
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T4);
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T5);
		svmopa_za32_m(3, pgIdx, pg, lambda, T6);
		svmopa_za32_m(3, pgIdx, pg, lambda, T7);
		svmopa_za32_m(3, pgIdx, pg, lambda, T8);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 5); //dTzz/dt
		svmopa_za32_m(3, pgIdx, pg, lambda, T0);
		svmopa_za32_m(3, pgIdx, pg, lambda, T1);
		svmopa_za32_m(3, pgIdx, pg, lambda, T2);
		svmopa_za32_m(3, pgIdx, pg, lambda, T3);
		svmopa_za32_m(3, pgIdx, pg, lambda, T4);
		svmopa_za32_m(3, pgIdx, pg, lambda, T5);
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T6);
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T7);
		svmopa_za32_m(3, pgIdx, pg, lam2mu, T8);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 6); //dTxx/dt  
		T0 = svmul_m(pg, Vx_xi, xi_y);
		T1 = svmul_m(pg, Vy_xi, xi_x);
		T2 = svmul_m(pg, Vx_et, et_y);
		T3 = svmul_m(pg, Vy_et, et_x);
		T4 = svmul_m(pg, Vx_zt, zt_y);
		T5 = svmul_m(pg, Vy_zt, zt_x);

		svmopa_za32_m(3, pgIdx, pg, mu, T0);
		svmopa_za32_m(3, pgIdx, pg, mu, T1);
		svmopa_za32_m(3, pgIdx, pg, mu, T2);
		svmopa_za32_m(3, pgIdx, pg, mu, T3);
		svmopa_za32_m(3, pgIdx, pg, mu, T4);
		svmopa_za32_m(3, pgIdx, pg, mu, T5);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 7); //dTyy/dt
		T0 = svmul_m(pg, Vx_xi, xi_z);
		T1 = svmul_m(pg, Vz_xi, xi_x);
		T2 = svmul_m(pg, Vx_et, et_z);
		T3 = svmul_m(pg, Vz_et, et_x);
		T4 = svmul_m(pg, Vx_zt, zt_z);
		T5 = svmul_m(pg, Vz_zt, zt_x);

		svmopa_za32_m(3, pgIdx, pg, mu, T0);
		svmopa_za32_m(3, pgIdx, pg, mu, T1);
		svmopa_za32_m(3, pgIdx, pg, mu, T2);
		svmopa_za32_m(3, pgIdx, pg, mu, T3);
		svmopa_za32_m(3, pgIdx, pg, mu, T4);
		svmopa_za32_m(3, pgIdx, pg, mu, T5);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 8); //dTzz/dt
		T0 = svmul_m(pg, Vy_xi, xi_z);
		T1 = svmul_m(pg, Vz_xi, xi_y);
		T2 = svmul_m(pg, Vy_et, et_z);
		T3 = svmul_m(pg, Vz_et, et_y);
		T4 = svmul_m(pg, Vy_zt, zt_z);
		T5 = svmul_m(pg, Vz_zt, zt_y);

		svmopa_za32_m(3, pgIdx, pg, mu, T0);
		svmopa_za32_m(3, pgIdx, pg, mu, T1);
		svmopa_za32_m(3, pgIdx, pg, mu, T2);
		svmopa_za32_m(3, pgIdx, pg, mu, T3);
		svmopa_za32_m(3, pgIdx, pg, mu, T4);
		svmopa_za32_m(3, pgIdx, pg, mu, T5);


		T0 = svread_hor_za32_f32_m(T0, pg, 3, 0);
		T1 = svread_hor_za32_f32_m(T1, pg, 3, 1);
		T2 = svread_hor_za32_f32_m(T2, pg, 3, 2);
		T3 = svread_hor_za32_f32_m(T3, pg, 3, 3);
		T4 = svread_hor_za32_f32_m(T4, pg, 3, 4);
		T5 = svread_hor_za32_f32_m(T5, pg, 3, 5);
		T6 = svread_hor_za32_f32_m(T6, pg, 3, 6);
		T7 = svread_hor_za32_f32_m(T7, pg, 3, 7);
		T8 = svread_hor_za32_f32_m(T8, pg, 3, 8);

		T0 = svmul_m(pg, T0, DT);
		T1 = svmul_m(pg, T1, DT);
		T2 = svmul_m(pg, T2, DT);
		T3 = svmul_m(pg, T3, DT);
		T4 = svmul_m(pg, T4, DT);
		T5 = svmul_m(pg, T5, DT);
		T6 = svmul_m(pg, T6, DT);
		T7 = svmul_m(pg, T7, DT);
		T8 = svmul_m(pg, T8, DT);

		svst1(pg, &(h_W[index+num*0]), T0);						
		svst1(pg, &(h_W[index+num*1]), T1);						
		svst1(pg, &(h_W[index+num*2]), T2);						
		svst1(pg, &(h_W[index+num*3]), T3);						
		svst1(pg, &(h_W[index+num*4]), T4);						
		svst1(pg, &(h_W[index+num*5]), T5);						
		svst1(pg, &(h_W[index+num*6]), T6);						
		svst1(pg, &(h_W[index+num*7]), T7);						
		svst1(pg, &(h_W[index+num*8]), T8);						


	END_CALCULATE3D( )												
}			


void waveDeriv( GRID grid, WAVE wave, FLOAT * CJM, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int FB1, int FB2, int FB3, float DT )	
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	float rDH = grid.rDH;

#ifdef GPU_CUDA
	int nx = _nx_ - 2 * HALO;
	int ny = _ny_ - 2 * HALO;
	int nz = _nz_ - 2 * HALO;


#ifdef XFAST
	dim3 threads( 32, 4, 4);
#endif
#ifdef ZFAST
	dim3 threads( 1, 8, 64);
#endif
	dim3 blocks;
	blocks.x = ( nx + threads.x - 1 ) / threads.x;
	blocks.y = ( ny + threads.y - 1 ) / threads.y;
	blocks.z = ( nz + threads.z - 1 ) / threads.z;
	
	//cout << "X = " << blocks.x << "Y = " << blocks.y << "Z = " << blocks.z << endl;

	wave_deriv <<< blocks, threads >>>
	( wave.h_W, wave.W, CJM, 
#ifdef PML
	pml_beta,
#endif  //PML
	_nx_, _ny_, _nz_, rDH, FB1, FB2, FB3, DT );	

	CHECK( cudaDeviceSynchronize( ) );

#else   //GPU_CUDA

	wave_deriv 
	( wave.h_W, wave.W, CJM, 
#ifdef PML
	pml_beta,
#endif //PML
	_nx_, _ny_, _nz_, rDH, FB1, FB2, FB3, DT );	

#endif //GPU_CUDA
}

