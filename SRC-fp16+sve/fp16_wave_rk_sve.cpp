#include "header.h"
/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:propagate.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-03
*   Discription:
*
================================================================*/
#include <arm_sve.h>
typedef void (*WAVE_RK_FUNC_FLOAT )( FLOAT * h_W, FLOAT * W, FLOAT * t_W, FLOAT * m_W, uint64_t num );
__GLOBAL__
void wave_rk0( FLOAT * h_W, FLOAT * W, FLOAT * t_W, FLOAT * m_W, uint64_t WStride )
{
	//printf( "wave_rk0\n" );
#ifdef GPU_CUDA
	uint64_t i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	uint64_t i = 0;
#endif
	//float h_w, w, t_w, m_w;
	svfloat16_t h_w, w, t_w, m_w;
		
	svfloat16_t svbeta1  = svdup_n_f16(beta1);
	svfloat16_t svalpha2 = svdup_n_f16(alpha2);

	uint64_t len = svcnth();
	svbool_t pg;

	CALCULATE1D_SVE( i, len, 0, WStride )
		pg = svwhilelt_b16( i, WStride );
		h_w = svld1(pg, &(h_W[i]));
		m_w = svld1(pg, &(W[i]));
		
		t_w = svmad_m(pg, h_w, svbeta1,  m_w);
		w   = svmad_m(pg, h_w, svalpha2, m_w);

		//t_w = m_w + beta1  * h_w;
		//w   = m_w + alpha2 * h_w;

		svst1(pg, &(m_W[i]), m_w);
		svst1(pg, &(t_W[i]), t_w);
		svst1(pg, &(W[i]), w);
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_rk1( FLOAT * h_W, FLOAT * W, FLOAT * t_W, FLOAT * m_W, uint64_t WStride)
{
	//printf( "wave_rk1\n" );
#ifdef GPU_CUDA
	uint64_t i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	uint64_t i = 0;
#endif
	svfloat16_t h_w, w, t_w, m_w;
	uint64_t len = svcnth();
	svbool_t pg;
	svfloat16_t svbeta2  = svdup_n_f16(beta2);
	svfloat16_t svalpha3 = svdup_n_f16(alpha3);
	CALCULATE1D_SVE( i, len, 0, WStride )
		pg = svwhilelt_b16( i, WStride );
		h_w = svld1(pg, &(h_W[i]));
		t_w = svld1(pg, &(t_W[i]));
		m_w = svld1(pg, &(m_W[i]));
		
		//t_w += beta2 * h_w;
		//w   = m_w + alpha3 * h_w;
		t_w = svmad_m(pg, h_w, svbeta2,  t_w);	
		w   = svmad_m(pg, h_w, svalpha3, m_w);

		//t_W[i] = t_w;
		//W[i] = w;
		svst1(pg, &(t_W[i]), t_w);
		svst1(pg, &(W[i]), w);
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_rk2( FLOAT * h_W, FLOAT * W, FLOAT * t_W, FLOAT * m_W, uint64_t WStride )
{
	//printf( "wave_rk2\n" );
#ifdef GPU_CUDA
	uint64_t i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	uint64_t i = 0;
#endif
	svfloat16_t h_w, w, t_w, m_w;
	uint64_t len = svcnth();
	svbool_t pg;
	svfloat16_t svbeta3 = svdup_n_f16(beta3);
	CALCULATE1D_SVE( i, len, 0, WStride )
		pg = svwhilelt_b16( i, WStride );
		h_w = svld1(pg, &(h_W[i]));
		t_w = svld1(pg, &(t_W[i]));
		m_w = svld1(pg, &(m_W[i]));

		t_w = svmad_m(pg, h_w, svbeta3,  t_w);	
		w = svadd_m(pg, m_w, h_w);

		svst1(pg, &(t_W[i]), t_w);
		svst1(pg, &(W[i]), w);
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_rk3( FLOAT * h_W, FLOAT * W, FLOAT * t_W, FLOAT * m_W, uint64_t WStride )
{
	//printf( "wave_rk3\n" );
#ifdef GPU_CUDA
	uint64_t i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	uint64_t i = 0;
#endif
	svfloat16_t h_w, w, t_w, m_w;
	uint64_t len = svcnth();
	svbool_t pg;
	svfloat16_t svbeta4  = svdup_n_f16(beta4);
	CALCULATE1D_SVE( i, len, 0, WStride )
		pg = svwhilelt_b16( i, WStride );
		h_w = svld1(pg, &(h_W[i]));
		t_w = svld1(pg, &(t_W[i]));

		w = svmad_m(pg, h_w, svbeta4,  t_w);	

		svst1(pg, &(W[i]), w);
	END_CALCULATE1D( )
}

void waveRk( GRID grid, int irk, WAVE wave )
{
	WAVE_RK_FUNC_FLOAT wave_rk[4] = { wave_rk0, wave_rk1, wave_rk2, wave_rk3 };
	uint64_t num = grid._nx_ * grid._ny_ * grid._nz_ * WSIZE;

#ifdef GPU_CUDA
	dim3 threads( 1024, 1, 1 );
	dim3 blocks;
	blocks.x = ( num + threads.x - 1 ) / threads.x;
	blocks.y = 1;
	blocks.z = 1;
	wave_rk[irk]<<< blocks, threads >>>( wave.h_W, wave.W, wave.t_W, wave.m_W, num );
	//wave_rk<<< blocks, threads >>>( h_W, W, t_W, m_W, num, DT );
	CHECK( cudaDeviceSynchronize( ) );
#else
	wave_rk[irk]( wave.h_W, wave.W, wave.t_W, wave.m_W, num );
#endif

}
