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
	float h_w, w, t_w, m_w;
	CALCULATE1D( i, 0, WStride ) 
		h_w = (float)h_W[i];
		m_w = (float)W[i];

		t_w = m_w + beta1  * h_w;
		w   = m_w + alpha2 * h_w;

		m_W[i] = m_w;
		t_W[i] = t_w;
		W[i] = w;
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
	float h_w, w, t_w, m_w;
	CALCULATE1D( i, 0, WStride ) 
		h_w = (float)h_W[i];
		t_w = (float)t_W[i];
		m_w = (float)m_W[i];

		t_w += beta2 * h_w;
		w   = m_w + alpha3 * h_w;
		
		t_W[i] = t_w;
		W[i] = w;
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
	float h_w, w, t_w, m_w;
	CALCULATE1D( i, 0, WStride ) 
		h_w = (float)h_W[i];
		t_w = (float)t_W[i];
		m_w = (float)m_W[i];

		t_w += beta3 * h_w;
		w    = m_w + h_w;

		t_W[i] = t_w;
		W[i] = w;
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
	float h_w, w, t_w, m_w;
	CALCULATE1D( i, 0, WStride ) 
		h_w = (float)h_W[i];
		t_w = (float)t_W[i];

		w = t_w +  beta4 * h_w;

		W[i] = w;
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
