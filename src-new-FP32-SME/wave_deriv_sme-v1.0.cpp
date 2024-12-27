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

		L_SME( Vx_xi, W, 0, FB1, xi, 0 ); 				
		L_SME( Vy_xi, W, 1, FB1, xi, 0 ); 				
		L_SME( Vz_xi, W, 2, FB1, xi, 0 ); 				
		L_SME(Txx_xi, W, 3, FB1, xi, 0 ); 				
		L_SME(Tyy_xi, W, 4, FB1, xi, 0 ); 				
		L_SME(Tzz_xi, W, 5, FB1, xi, 0 ); 				
		L_SME(Txy_xi, W, 6, FB1, xi, 0 ); 				
		L_SME(Txz_xi, W, 7, FB1, xi, 0 ); 				
		L_SME(Tyz_xi, W, 8, FB1, xi, 0 ); 				

		L_SME( Vx_et, W, 0, FB2, et, 1 );				
		L_SME( Vy_et, W, 1, FB2, et, 1 );				
		L_SME( Vz_et, W, 2, FB2, et, 1 );				
		L_SME(Txx_et, W, 3, FB2, et, 1 );				
		L_SME(Tyy_et, W, 4, FB2, et, 1 );				
		L_SME(Tzz_et, W, 5, FB2, et, 1 );				
		L_SME(Txy_et, W, 6, FB2, et, 1 );				
		L_SME(Txz_et, W, 7, FB2, et, 1 );				
		L_SME(Tyz_et, W, 8, FB2, et, 1 );				

   		L_SME( Vx_zt, W, 0, FB3, zt, 2 );				
   		L_SME( Vy_zt, W, 1, FB3, zt, 2 );				
   		L_SME( Vz_zt, W, 2, FB3, zt, 2 );				
  		L_SME(Txx_zt, W, 3, FB3, zt, 2 );				
  		L_SME(Tyy_zt, W, 4, FB3, zt, 2 );				
  		L_SME(Tzz_zt, W, 5, FB3, zt, 2 );				
  		L_SME(Txy_zt, W, 6, FB3, zt, 2 );				
  		L_SME(Txz_zt, W, 7, FB3, zt, 2 );				
  		L_SME(Tyz_zt, W, 8, FB3, zt, 2 );

#ifdef PML
		pml_beta_x = svld1( pg, &(pml_beta.x[i]) );
		pml_beta_y = pml_beta.y[j];
		pml_beta_z = pml_beta.z[k];
#endif


#ifdef PML
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 0); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 0, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 1); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 1, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 2); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 2, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 3); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 3, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 4); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 4, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 5); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 5, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 6); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 6, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 7); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 7, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 0, 8); TMP = svmul_m( pg, TMP, pml_beta_x ); svwrite_hor_za32_m(0, 8, pg,TMP);
                                                                                                                         
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 0); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 0, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 1); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 1, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 2); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 2, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 3); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 3, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 4); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 4, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 5); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 5, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 6); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 6, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 7); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 7, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 1, 8); TMP = svmul_m( pg, TMP, pml_beta_y ); svwrite_hor_za32_m(1, 8, pg,TMP);

		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 0); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 0, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 1); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 1, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 2); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 2, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 3); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 3, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 4); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 4, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 5); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 5, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 6); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 6, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 7); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 7, pg,TMP);
		TMP = svread_hor_za32_f32_m(TMP, pg, 2, 8); TMP = svmul_m( pg, TMP, pml_beta_z ); svwrite_hor_za32_m(2, 8, pg,TMP);
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
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 3)/*Txx_xi*/, xi_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 6)/*Txy_xi*/, xi_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 7)/*Txz_xi*/, xi_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 3)/*Txx_et*/, et_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 6)/*Txy_et*/, et_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 7)/*Txz_et*/, et_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 3)/*Txx_zt*/, zt_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 6)/*Txy_zt*/, zt_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 7)/*Txz_zt*/, zt_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 1); //dVy/dt
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 6)/*Txy_xi*/, xi_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 4)/*Tyy_xi*/, xi_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 8)/*Tyz_xi*/, xi_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 6)/*Txy_et*/, et_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 4)/*Tyy_et*/, et_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 8)/*Tyz_et*/, et_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 6)/*Txy_zt*/, zt_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 4)/*Tyy_zt*/, zt_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 8)/*Tyz_zt*/, zt_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);

		pgIdx = svcmpeq(svptrue_b32(), svindex_s32(0, 1), 2); //dVz/dt  
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 7)/*Txz_xi*/, xi_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 8)/*Tyz_xi*/, xi_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 0, 5)/*Tzz_xi*/, xi_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 7)/*Txz_et*/, et_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 8)/*Tyz_et*/, et_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 1, 5)/*Tzz_et*/, et_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 7)/*Txz_zt*/, zt_x); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 8)/*Tyz_zt*/, zt_y); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		TMP = svmul_m(pg, svread_hor_za32_f32_m(TMP, pg, 2, 5)/*Tzz_zt*/, zt_z); svmopa_za32_m(3, pgIdx, pg, buoyancy, TMP);
		

		//Tile 0:xi 1:et 2:zt 3:d /dt
		lam2mu = svmul_m(pg, mu, 2.0f);
		lam2mu = svadd_m(pg, lam2mu, lambda);
		svfloat32_t Vx_xi = svread_hor_za32_f32_m(Vx_xi, pg, 0, 0);/**/
        svfloat32_t Vx_et = svread_hor_za32_f32_m(Vx_et, pg, 1, 0);/**/
        svfloat32_t Vx_zt = svread_hor_za32_f32_m(Vx_zt, pg, 2, 0);/**/
        svfloat32_t Vy_xi = svread_hor_za32_f32_m(Vy_xi, pg, 0, 1);/**/
        svfloat32_t Vy_et = svread_hor_za32_f32_m(Vy_et, pg, 1, 1);/**/
        svfloat32_t Vy_zt = svread_hor_za32_f32_m(Vy_zt, pg, 2, 1);/**/
        svfloat32_t Vz_xi = svread_hor_za32_f32_m(Vz_xi, pg, 0, 2);/**/
        svfloat32_t Vz_et = svread_hor_za32_f32_m(Vz_et, pg, 1, 2);/**/
        svfloat32_t Vz_zt = svread_hor_za32_f32_m(Vz_zt, pg, 2, 2);/**/

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

		/*
		Z0 = svread_hor_za32_f32_m(TMP, pg, 3, 3); 
		Z1 = svread_hor_za32_f32_m(TMP, pg, 3, 4);
		Z2 = svread_hor_za32_f32_m(TMP, pg, 3, 5);

		TMP = svadd_m(pg, Z0, Z1); 
		TMP = svadd_m(pg, TMP, Z2); 
		TMP = svmul_m(pg, TMP, lambda);

		Z0 = svmul_m(pg, Z0, mu); Z0 = svmul_m(pg, Z0, 2.0f);
		Z1 = svmul_m(pg, Z1, mu); Z1 = svmul_m(pg, Z1, 2.0f);
		Z2 = svmul_m(pg, Z2, mu); Z2 = svmul_m(pg, Z2, 2.0f);

		Z0 = svadd_m(pg, TMP, Z0);	svwrite_hor_za32_m(3, 3, pg, Z0);
		Z1 = svadd_m(pg, TMP, Z1);	svwrite_hor_za32_m(3, 4, pg, Z1);
		Z2 = svadd_m(pg, TMP, Z2);	svwrite_hor_za32_m(3, 5, pg, Z2);
		*/
		

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

