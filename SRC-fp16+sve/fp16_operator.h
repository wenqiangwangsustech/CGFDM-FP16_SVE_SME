/***********************************************************
  File Name: fp16_operator.h
  Author: Wenqiang Wang
  Mail: 11849528@mail.sustech.edu.cn
  Created Time: Fri 29 Nov 2024 03:08:46 PM CST
 *******************************************************/

#ifndef FP16_SVE
#define FP16_SVE


#define L_SVE_FP16( RET, W, VAR, VARSIZE, FB, SUB )																	\
offset = svindex_u32(0, VARSIZE * sizeof(FLOAT));																	\
sveInt32 = svld1sh_gather_offset_s32(pg, (int16_t*)&(W[INDEX_##SUB( i, j, k, - FB * 1) * VARSIZE + VAR]), offset);	\
sveArr16 = svreinterpret_f16( sveInt32 );																			\
Z_1 = svcvt_f32_m(TMP, pg, sveArr16 );																				\
																													\
sveInt32 = svld1sh_gather_offset_s32(pg, (int16_t*)&(W[index * VARSIZE + VAR]), offset);							\
sveArr16 = svreinterpret_f16( sveInt32 );																			\
Z0  = svcvt_f32_m(TMP, pg, sveArr16 );																				\
																													\
sveInt32 = svld1sh_gather_offset_s32(pg, (int16_t*)&(W[INDEX_##SUB( i, j, k, + FB * 1) * VARSIZE + VAR]), offset);	\
sveArr16 = svreinterpret_f16( sveInt32 );																			\
Z1  = svcvt_f32_m(TMP, pg, sveArr16 );																				\
																													\
sveInt32 = svld1sh_gather_offset_s32(pg, (int16_t*)&(W[INDEX_##SUB( i, j, k, + FB * 2) * VARSIZE + VAR]), offset);	\
sveArr16 = svreinterpret_f16( sveInt32 );																			\
Z2 = svcvt_f32_m(TMP, pg, sveArr16 );																				\
																													\
sveInt32 = svld1sh_gather_offset_s32(pg, (int16_t*)&(W[INDEX_##SUB( i, j, k, + FB * 3) * VARSIZE + VAR]), offset);	\
sveArr16 = svreinterpret_f16( sveInt32 );																			\
Z3 = svcvt_f32_m(TMP, pg, sveArr16 );															\
																								\
Z_1 = svmul_m(pg, Z_1, (rDH * af_1 * FB) );														\
Z0  = svmul_m(pg, Z0 , (rDH * af0  * FB) );														\
Z1  = svmul_m(pg, Z1 , (rDH * af1  * FB) );														\
Z2  = svmul_m(pg, Z2 , (rDH * af2  * FB) );														\
Z3  = svmul_m(pg, Z3 , (rDH * af3  * FB) );														\
																								\
RET = svadd_m( pg, Z_1, Z0 );																	\
RET = svadd_m( pg, RET, Z1 );																	\
RET = svadd_m( pg, RET, Z2 );																	\
RET = svadd_m( pg, RET, Z3 );																	



#endif //FP16_SVE



