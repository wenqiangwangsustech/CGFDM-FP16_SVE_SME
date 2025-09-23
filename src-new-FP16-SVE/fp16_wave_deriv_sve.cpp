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
#include "fp16_operator.h"

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

svuint16_t setSVE16LDIdx()
{
	svbool_t pg = svptrue_b16();
	int len = svcntb() / sizeof(__fp16);
	uint16_t pos[len];
	for (uint64_t pi = 0; pi < len/2; pi ++ )
	{
		pos[pi * 2] = pi;
		pos[pi * 2+1] = pi + len/2;
	}
	svuint16_t sveLD = svld1(pg, pos);
	return sveLD;
}

svuint16_t setSVE16STIdx()
{
	svbool_t pg = svptrue_b16();
	int len = svcntb() / sizeof(__fp16);
	uint16_t pos[len];
	for (uint64_t pi = 0; pi < len/2; pi ++ )
	{
		pos[pi] = pi * 2;
		pos[pi+len/2] = pi*2+1;
	}
	svuint16_t sveST = svld1(pg, pos);
	return sveST;
}

__GLOBAL__															
void wave_deriv( FLOAT * h_W, FLOAT * W, FLOAT * CJM, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int _nx_, int _ny_, int _nz_, float rDH, int FB1, int FB2, int FB3, float DT )	
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
	svfloat32_t Vx1;	svfloat32_t Vx2;	svfloat32_t Vx3;
	svfloat32_t Vy1;	svfloat32_t Vy2;	svfloat32_t Vy3;
	svfloat32_t Vz1;	svfloat32_t Vz2;	svfloat32_t Vz3;
	svfloat32_t Txx1;	svfloat32_t Txx2;   svfloat32_t Txx3;
	svfloat32_t Tyy1;	svfloat32_t Tyy2;	svfloat32_t Tyy3;
	svfloat32_t Tzz1;   svfloat32_t Tzz2;	svfloat32_t Tzz3;
	svfloat32_t Txy1;   svfloat32_t Txy2;	svfloat32_t Txy3;
	svfloat32_t Txz1;   svfloat32_t Txz2;	svfloat32_t Txz3;
	svfloat32_t Tyz1;   svfloat32_t Tyz2;	svfloat32_t Tyz3;

	svfloat32_t Tx1;
	svfloat32_t Tx2;
	svfloat32_t Tx3;
	svfloat32_t Ty1;
	svfloat32_t Ty2;
	svfloat32_t Ty3;
	svfloat32_t Tz1;
	svfloat32_t Tz2;
	svfloat32_t Tz3;


	svfloat32_t HW0;
	svfloat32_t HW1;
	svfloat32_t HW2;
	svfloat32_t HW3;
	svfloat32_t HW4;
	svfloat32_t HW5;
	svfloat32_t HW6;
	svfloat32_t HW7;
	svfloat32_t HW8;


	int len = svcntw();

	svfloat32_t Z_1, Z0, Z1, Z2, Z3, Z4;

	svbool_t pg;
	svbool_t pgld = svwhilelt_b16(0, len);
	svbool_t pg16 = svptrue_b16();

	svfloat16_t sveFP16;
	svuint16_t sveLD = setSVE16LDIdx();
	svuint16_t sveST = setSVE16STIdx();
		
	long long num = _nx_ * _ny_ * _nz_; 
	
	CALCULATE3D_SVE( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )  
	
	
	
	END_CALCULATE3D( )												

	CALCULATE3D_SVE( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )  
		pg = svwhilelt_b32( i, _nx );
		pgld = svwhilelt_b16(i, _nx);
		if (_nx - i > len)
			pgld = svwhilelt_b16(0, len);
		
		index = INDEX( i, j, k );

		L_SVE_FP16( Vx_xi, W, 0, FB1, xi ); 				
		L_SVE_FP16( Vy_xi, W, 1, FB1, xi ); 				
		L_SVE_FP16( Vz_xi, W, 2, FB1, xi ); 				
		L_SVE_FP16(Txx_xi, W, 3, FB1, xi ); 				
		L_SVE_FP16(Tyy_xi, W, 4, FB1, xi ); 				
		L_SVE_FP16(Tzz_xi, W, 5, FB1, xi ); 				
		L_SVE_FP16(Txy_xi, W, 6, FB1, xi ); 				
		L_SVE_FP16(Txz_xi, W, 7, FB1, xi ); 				
		L_SVE_FP16(Tyz_xi, W, 8, FB1, xi ); 				

		L_SVE_FP16( Vx_et, W, 0, FB2, et );				
		L_SVE_FP16( Vy_et, W, 1, FB2, et );				
		L_SVE_FP16( Vz_et, W, 2, FB2, et );				
		L_SVE_FP16(Txx_et, W, 3, FB2, et );				
		L_SVE_FP16(Tyy_et, W, 4, FB2, et );				
		L_SVE_FP16(Tzz_et, W, 5, FB2, et );				
		L_SVE_FP16(Txy_et, W, 6, FB2, et );				
		L_SVE_FP16(Txz_et, W, 7, FB2, et );				
		L_SVE_FP16(Tyz_et, W, 8, FB2, et );				

   		L_SVE_FP16( Vx_zt, W, 0, FB3, zt );				
   		L_SVE_FP16( Vy_zt, W, 1, FB3, zt );				
   		L_SVE_FP16( Vz_zt, W, 2, FB3, zt );				
  		L_SVE_FP16(Txx_zt, W, 3, FB3, zt );				
  		L_SVE_FP16(Tyy_zt, W, 4, FB3, zt );				
  		L_SVE_FP16(Tzz_zt, W, 5, FB3, zt );				
  		L_SVE_FP16(Txy_zt, W, 6, FB3, zt );				
  		L_SVE_FP16(Txz_zt, W, 7, FB3, zt );				
  		L_SVE_FP16(Tyz_zt, W, 8, FB3, zt );

#ifdef PML
		pml_beta_x = svld1( pg, &(pml_beta.x[i]) );
		pml_beta_y = pml_beta.y[j];
		pml_beta_z = pml_beta.z[k];
#endif


#ifdef PML
		 //Vx_xi = svmul_m( pg,  Vx_xi, 1.0f/*pml_beta_x*/ );
		 //Vy_xi = svmul_m( pg,  Vy_xi, 1.0f/*pml_beta_x*/ );
		 //Vz_xi = svmul_m( pg,  Vz_xi, 1.0f/*pml_beta_x*/ );
		 //printf("I am wangwq");
		 Vx_xi = svmul_m( pg,  Vx_xi, pml_beta_x );
		 Vy_xi = svmul_m( pg,  Vy_xi, pml_beta_x );
		 Vz_xi = svmul_m( pg,  Vz_xi, pml_beta_x );
		 Vx_et = svmul_m( pg,  Vx_et, pml_beta_y );
		 Vy_et = svmul_m( pg,  Vy_et, pml_beta_y );
		 Vz_et = svmul_m( pg,  Vz_et, pml_beta_y );
		 Vx_zt = svmul_m( pg,  Vx_zt, pml_beta_z );
		 Vy_zt = svmul_m( pg,  Vy_zt, pml_beta_z );
		 Vz_zt = svmul_m( pg,  Vz_zt, pml_beta_z );

		//Txx_xi = svmul_m( pg, Txx_xi, 1.0f/*pml_beta_x*/ );
		//Tyy_xi = svmul_m( pg, Tyy_xi, 1.0f/*pml_beta_x*/ );
		//Tzz_xi = svmul_m( pg, Tzz_xi, 1.0f/*pml_beta_x*/ );
		//Txy_xi = svmul_m( pg, Txy_xi, 1.0f/*pml_beta_x*/ );
		//Txz_xi = svmul_m( pg, Txz_xi, 1.0f/*pml_beta_x*/ );
		//Tyz_xi = svmul_m( pg, Tyz_xi, 1.0f/*pml_beta_x*/ );
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

		sveFP16 = svld1(pgld, &(CJM[index+num*0])); sveFP16 = svtbl(sveFP16, sveLD); xi_x = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*1])); sveFP16 = svtbl(sveFP16, sveLD); xi_y = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*2])); sveFP16 = svtbl(sveFP16, sveLD); xi_z = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*3])); sveFP16 = svtbl(sveFP16, sveLD); et_x = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*4])); sveFP16 = svtbl(sveFP16, sveLD); et_y = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*5])); sveFP16 = svtbl(sveFP16, sveLD); et_z = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*6])); sveFP16 = svtbl(sveFP16, sveLD); zt_x = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*7])); sveFP16 = svtbl(sveFP16, sveLD); zt_y = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*8])); sveFP16 = svtbl(sveFP16, sveLD); zt_z = svcvt_f32_x(pg, sveFP16);

		sveFP16 = svld1(pgld, &(CJM[index+num*10])); sveFP16 = svtbl(sveFP16, sveLD); mu 	   = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*11])); sveFP16 = svtbl(sveFP16, sveLD); lambda   = svcvt_f32_x(pg, sveFP16);
		sveFP16 = svld1(pgld, &(CJM[index+num*12])); sveFP16 = svtbl(sveFP16, sveLD); buoyancy = svcvt_f32_x(pg, sveFP16);
		buoyancy = svmul_m(pg, buoyancy, Crho);

		DOT_PRODUCT3D_SVE(Vx1, xi_x, xi_y, xi_z, Txx_xi, Txy_xi, Txz_xi );
		DOT_PRODUCT3D_SVE(Vx2, et_x, et_y, et_z, Txx_et, Txy_et, Txz_et );
		DOT_PRODUCT3D_SVE(Vx3, zt_x, zt_y, zt_z, Txx_zt, Txy_zt, Txz_zt );
		DOT_PRODUCT3D_SVE(Vy1, xi_x, xi_y, xi_z, Txy_xi, Tyy_xi, Tyz_xi );
		DOT_PRODUCT3D_SVE(Vy2, et_x, et_y, et_z, Txy_et, Tyy_et, Tyz_et );
		DOT_PRODUCT3D_SVE(Vy3, zt_x, zt_y, zt_z, Txy_zt, Tyy_zt, Tyz_zt );
		DOT_PRODUCT3D_SVE(Vz1, xi_x, xi_y, xi_z, Txz_xi, Tyz_xi, Tzz_xi );
		DOT_PRODUCT3D_SVE(Vz2, et_x, et_y, et_z, Txz_et, Tyz_et, Tzz_et );
		DOT_PRODUCT3D_SVE(Vz3, zt_x, zt_y, zt_z, Txz_zt, Tyz_zt, Tzz_zt );

		Vx1 = svmul_m(pg, Vx1, buoyancy);
		Vx2 = svmul_m(pg, Vx2, buoyancy);
		Vx3 = svmul_m(pg, Vx3, buoyancy);
		Vy1 = svmul_m(pg, Vy1, buoyancy);
		Vy2 = svmul_m(pg, Vy2, buoyancy);
		Vy3 = svmul_m(pg, Vy3, buoyancy);
		Vz1 = svmul_m(pg, Vz1, buoyancy);
		Vz2 = svmul_m(pg, Vz2, buoyancy);
		Vz3 = svmul_m(pg, Vz3, buoyancy);

		DOT_PRODUCT3D_SVE(Tx1, xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi );
		DOT_PRODUCT3D_SVE(Tx2, et_x, et_y, et_z, Vx_et, Vy_et, Vz_et );
		DOT_PRODUCT3D_SVE(Tx3, zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt );
		DOT_PRODUCT3D_SVE(Ty1, xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi );
		DOT_PRODUCT3D_SVE(Ty2, et_x, et_y, et_z, Vx_et, Vy_et, Vz_et );
		DOT_PRODUCT3D_SVE(Ty3, zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt );
		DOT_PRODUCT3D_SVE(Tz1, xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi );
		DOT_PRODUCT3D_SVE(Tz2, et_x, et_y, et_z, Vx_et, Vy_et, Vz_et );
		DOT_PRODUCT3D_SVE(Tz3, zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt );

		Txx1 = svmul_m(pg, Tx1, lambda);
		Txx2 = svmul_m(pg, Tx2, lambda);
		Txx3 = svmul_m(pg, Tx3, lambda);
		Tyy1 = svmul_m(pg, Ty1, lambda);
		Tyy2 = svmul_m(pg, Ty2, lambda);
		Tyy3 = svmul_m(pg, Ty3, lambda);
		Tzz1 = svmul_m(pg, Tz1, lambda);
		Tzz2 = svmul_m(pg, Tz2, lambda);
		Tzz3 = svmul_m(pg, Tz3, lambda);

		Tx1 = svmul_m(pg, xi_x, Vx_xi);
		Tx2 = svmul_m(pg, et_x, Vx_et);
		Tx3 = svmul_m(pg, zt_x, Vx_zt);
		Ty1 = svmul_m(pg, xi_y, Vy_xi);
		Ty2 = svmul_m(pg, et_y, Vy_et);
		Ty3 = svmul_m(pg, zt_y, Vy_zt);
		Tz1 = svmul_m(pg, xi_z, Vz_xi);
		Tz2 = svmul_m(pg, et_z, Vz_et);
		Tz3 = svmul_m(pg, zt_z, Vz_zt);

		Tx1 = svmul_m(pg, Tx1, mu);	Tx1 = svmul_m(pg, Tx1, 2.0f);
		Tx2 = svmul_m(pg, Tx2, mu);	Tx2 = svmul_m(pg, Tx2, 2.0f);
		Tx3 = svmul_m(pg, Tx3, mu);	Tx3 = svmul_m(pg, Tx3, 2.0f);
		Ty1 = svmul_m(pg, Ty1, mu);	Ty1 = svmul_m(pg, Ty1, 2.0f);
		Ty2 = svmul_m(pg, Ty2, mu);	Ty2 = svmul_m(pg, Ty2, 2.0f);
		Ty3 = svmul_m(pg, Ty3, mu);	Ty3 = svmul_m(pg, Ty3, 2.0f);
		Tz1 = svmul_m(pg, Tz1, mu);	Tz1 = svmul_m(pg, Tz1, 2.0f);
		Tz2 = svmul_m(pg, Tz2, mu);	Tz2 = svmul_m(pg, Tz2, 2.0f);
		Tz3 = svmul_m(pg, Tz3, mu);	Tz3 = svmul_m(pg, Tz3, 2.0f);

		Txx1 = svadd_m(pg, Txx1, Tx1);
		Txx2 = svadd_m(pg, Txx2, Tx2);
		Txx3 = svadd_m(pg, Txx3, Tx3);
		Tyy1 = svadd_m(pg, Tyy1, Ty1);
		Tyy2 = svadd_m(pg, Tyy2, Ty2);
		Tyy3 = svadd_m(pg, Tyy3, Ty3);
		Tzz1 = svadd_m(pg, Tzz1, Tz1);
		Tzz2 = svadd_m(pg, Tzz2, Tz2);
		Tzz3 = svadd_m(pg, Tzz3, Tz3);

		DOT_PRODUCT2D_SVE(Txy1, xi_y, xi_x, Vx_xi, Vy_xi );
		DOT_PRODUCT2D_SVE(Txy2, et_y, et_x, Vx_et, Vy_et );
		DOT_PRODUCT2D_SVE(Txy3, zt_y, zt_x, Vx_zt, Vy_zt );
		DOT_PRODUCT2D_SVE(Txz1, xi_z, xi_x, Vx_xi, Vz_xi );
		DOT_PRODUCT2D_SVE(Txz2, et_z, et_x, Vx_et, Vz_et );
		DOT_PRODUCT2D_SVE(Txz3, zt_z, zt_x, Vx_zt, Vz_zt );
		DOT_PRODUCT2D_SVE(Tyz1, xi_z, xi_y, Vy_xi, Vz_xi );
		DOT_PRODUCT2D_SVE(Tyz2, et_z, et_y, Vy_et, Vz_et );
		DOT_PRODUCT2D_SVE(Tyz3, zt_z, zt_y, Vy_zt, Vz_zt );

		Txy1 = svmul_m(pg, Txy1, mu);
        Txy2 = svmul_m(pg, Txy2, mu);
        Txy3 = svmul_m(pg, Txy3, mu);
        Txz1 = svmul_m(pg, Txz1, mu);
        Txz2 = svmul_m(pg, Txz2, mu);
        Txz3 = svmul_m(pg, Txz3, mu);
        Tyz1 = svmul_m(pg, Tyz1, mu);
        Tyz2 = svmul_m(pg, Tyz2, mu);
        Tyz3 = svmul_m(pg, Tyz3, mu);
		
		SUM9VAR(HW0, Vx1 , Vx2 , Vx3 , DT);	
		SUM9VAR(HW1, Vy1 , Vy2 , Vy3 , DT);	
		SUM9VAR(HW2, Vz1 , Vz2 , Vz3 , DT);	
		SUM9VAR(HW3, Txx1, Txx2, Txx3, DT);	
		SUM9VAR(HW4, Tyy1, Tyy2, Tyy3, DT);	
		SUM9VAR(HW5, Tzz1, Tzz2, Tzz3, DT);	
		SUM9VAR(HW6, Txy1, Txy2, Txy3, DT);	
		SUM9VAR(HW7, Txz1, Txz2, Txz3, DT);	
		SUM9VAR(HW8, Tyz1, Tyz2, Tyz3, DT);	

		sveFP16 = svcvt_f16_x(pg16, HW0); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*0]), sveFP16); 
		sveFP16 = svcvt_f16_x(pg16, HW1); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*1]), sveFP16); 
		sveFP16 = svcvt_f16_x(pg16, HW2); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*2]), sveFP16); 
		sveFP16 = svcvt_f16_x(pg16, HW3); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*3]), sveFP16); 
		sveFP16 = svcvt_f16_x(pg16, HW4); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*4]), sveFP16); 
		sveFP16 = svcvt_f16_x(pg16, HW5); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*5]), sveFP16); 
		sveFP16 = svcvt_f16_x(pg16, HW6); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*6]), sveFP16); 
		sveFP16 = svcvt_f16_x(pg16, HW7); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*7]), sveFP16); 
		sveFP16 = svcvt_f16_x(pg16, HW8); sveFP16 = svtbl(sveFP16, sveST); svst1(pgld, &(h_W[index+num*8]), sveFP16); 
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

