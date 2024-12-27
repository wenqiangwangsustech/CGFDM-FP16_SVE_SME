/***********************************************************
  File Name: fp16_operator.h
  Author: Wenqiang Wang
  Mail: 11849528@mail.sustech.edu.cn
  Created Time: Fri 29 Nov 2024 03:08:46 PM CST
 *******************************************************/

#ifndef FP16_SME
#define FP16_SME

#define L_SME_FP16( W, FB, SUB )																					\
offset = svindex_u32(0, WSIZE * sizeof(FLOAT));																	\
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol0 , pg1, svdup_f32(rDH * af_1 * FB), Zn);	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol1 , pg2, svdup_f32(rDH * af_1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol2 , pg1, svdup_f32(rDH * af_1 * FB), Zn);	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol3 , pg2, svdup_f32(rDH * af_1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol4 , pg1, svdup_f32(rDH * af_1 * FB), Zn);	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol5 , pg2, svdup_f32(rDH * af_1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol6 , pg1, svdup_f32(rDH * af_1 * FB), Zn);	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol7 , pg2, svdup_f32(rDH * af_1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol8 , pg1, svdup_f32(rDH * af_1 * FB), Zn);	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol9 , pg2, svdup_f32(rDH * af_1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol10, pg1, svdup_f32(rDH * af_1 * FB), Zn);	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol11, pg2, svdup_f32(rDH * af_1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol12, pg1, svdup_f32(rDH * af_1 * FB), Zn);	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol13, pg2, svdup_f32(rDH * af_1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol14, pg1, svdup_f32(rDH * af_1 * FB), Zn);	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol15, pg2, svdup_f32(rDH * af_1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, - FB * 1) * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); Vz_##SUB##1 = svmad_m(pg1, svdup_f32(rDH * af_1 * FB), Zn, svdup_f32(0.0f));	\
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, - FB * 1) * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); Vz_##SUB##2 = svmad_m(pg2, svdup_f32(rDH * af_1 * FB), Zm, svdup_f32(0.0f)); \
\
\
\
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol0 , pg1, svdup_f32(rDH * af0 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol1 , pg2, svdup_f32(rDH * af0 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol2 , pg1, svdup_f32(rDH * af0 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol3 , pg2, svdup_f32(rDH * af0 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol4 , pg1, svdup_f32(rDH * af0 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol5 , pg2, svdup_f32(rDH * af0 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol6 , pg1, svdup_f32(rDH * af0 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol7 , pg2, svdup_f32(rDH * af0 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol8 , pg1, svdup_f32(rDH * af0 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol9 , pg2, svdup_f32(rDH * af0 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol10, pg1, svdup_f32(rDH * af0 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol11, pg2, svdup_f32(rDH * af0 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol12, pg1, svdup_f32(rDH * af0 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol13, pg2, svdup_f32(rDH * af0 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol14, pg1, svdup_f32(rDH * af0 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol15, pg2, svdup_f32(rDH * af0 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[index * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); Vz_##SUB##1 = svmad_m(pg1, svdup_f32(rDH * af0 * FB), Zn, Vz_##SUB##1);  \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[index * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); Vz_##SUB##2 = svmad_m(pg2, svdup_f32(rDH * af0 * FB), Zm, Vz_##SUB##2);  \
\
\
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol0 , pg1, svdup_f32(rDH * af1 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol1 , pg2, svdup_f32(rDH * af1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol2 , pg1, svdup_f32(rDH * af1 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol3 , pg2, svdup_f32(rDH * af1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol4 , pg1, svdup_f32(rDH * af1 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol5 , pg2, svdup_f32(rDH * af1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol6 , pg1, svdup_f32(rDH * af1 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol7 , pg2, svdup_f32(rDH * af1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol8 , pg1, svdup_f32(rDH * af1 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol9 , pg2, svdup_f32(rDH * af1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol10, pg1, svdup_f32(rDH * af1 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol11, pg2, svdup_f32(rDH * af1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol12, pg1, svdup_f32(rDH * af1 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol13, pg2, svdup_f32(rDH * af1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol14, pg1, svdup_f32(rDH * af1 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol15, pg2, svdup_f32(rDH * af1 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB) * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); Vz_##SUB##1 = svmad_m(pg1, svdup_f32(rDH * af1 * FB), Zn, Vz_##SUB##1);  \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB) * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); Vz_##SUB##2 = svmad_m(pg2, svdup_f32(rDH * af1 * FB), Zm, Vz_##SUB##2);  \
\
\
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol0 , pg1, svdup_f32(rDH * af2 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol1 , pg2, svdup_f32(rDH * af2 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol2 , pg1, svdup_f32(rDH * af2 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol3 , pg2, svdup_f32(rDH * af2 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol4 , pg1, svdup_f32(rDH * af2 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol5 , pg2, svdup_f32(rDH * af2 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol6 , pg1, svdup_f32(rDH * af2 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol7 , pg2, svdup_f32(rDH * af2 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol8 , pg1, svdup_f32(rDH * af2 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol9 , pg2, svdup_f32(rDH * af2 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol10, pg1, svdup_f32(rDH * af2 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol11, pg2, svdup_f32(rDH * af2 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol12, pg1, svdup_f32(rDH * af2 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol13, pg2, svdup_f32(rDH * af2 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol14, pg1, svdup_f32(rDH * af2 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol15, pg2, svdup_f32(rDH * af2 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 2) * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); Vz_##SUB##1 = svmad_m(pg1, svdup_f32(rDH * af2 * FB), Zn, Vz_##SUB##1);  \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 2) * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); Vz_##SUB##2 = svmad_m(pg2, svdup_f32(rDH * af2 * FB), Zm, Vz_##SUB##2);  \
\
\
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol0 , pg1, svdup_f32(rDH * af3 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 0]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol1 , pg2, svdup_f32(rDH * af3 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol2 , pg1, svdup_f32(rDH * af3 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 1]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol3 , pg2, svdup_f32(rDH * af3 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol4 , pg1, svdup_f32(rDH * af3 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 2]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol5 , pg2, svdup_f32(rDH * af3 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol6 , pg1, svdup_f32(rDH * af3 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 3]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol7 , pg2, svdup_f32(rDH * af3 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol8 , pg1, svdup_f32(rDH * af3 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 4]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol9 , pg2, svdup_f32(rDH * af3 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol10, pg1, svdup_f32(rDH * af3 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 5]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol11, pg2, svdup_f32(rDH * af3 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol12, pg1, svdup_f32(rDH * af3 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 6]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol13, pg2, svdup_f32(rDH * af3 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); svmopa_za32_m(0, pgCol14, pg1, svdup_f32(rDH * af3 * FB), Zn); \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 7]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); svmopa_za32_m(0, pgCol15, pg2, svdup_f32(rDH * af3 * FB), Zm); \
sveInt32 = svld1sh_gather_offset_s32(pg1, (int16_t*)&(W[INDEX_##SUB(i, 		   j, k, + FB * 3) * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zn = svcvt_f32_m(TMP, pg1, sveArr16); Vz_##SUB##1 = svmad_m(pg1, svdup_f32(rDH * af3 * FB), Zn, Vz_##SUB##1);  \
sveInt32 = svld1sh_gather_offset_s32(pg2, (int16_t*)&(W[INDEX_##SUB(i + len32, j, k, + FB * 3) * WSIZE + 8]), offset); sveArr16 = svreinterpret_f16(sveInt32); Zm = svcvt_f32_m(TMP, pg2, sveArr16); Vz_##SUB##2 = svmad_m(pg2, svdup_f32(rDH * af3 * FB), Zm, Vz_##SUB##2);  \



#endif //FP16_SME



