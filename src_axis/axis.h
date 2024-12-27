/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:macro.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2022-10-21
*   Discription:
*
================================================================*/


#define FOR_LOOP1D( i, startI, endI )  \
for ( i = startI; i < endI; i ++ ) { 

#define END_LOOP1D( ) }
#define END_LOOP2D( ) }}
#define END_LOOP3D( ) }}}



#ifdef XFAST

#define Index2D( i, j, nx, ny ) ( ( i ) + ( j ) * ( nx ) )

#define FOR_LOOP2D( i, j, startI, endI, startJ, endJ ) \
for ( j = startJ; j < endJ; j ++ ) {    \
for ( i = startI; i < endI; i ++ ) {    


#define FOR_LOOP3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( k = startK; k < endK; k ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( i = startI; i < endI; i ++ ) {     

#endif


//#define Index2D( i, j, nx, ny ) ( ( i ) + ( j ) * ( nx ) )




#ifdef ZFAST

#define Index2D( i, j, nx, ny ) ( ( j ) + ( i ) * ( ny ) )

#define FOR_LOOP2D( i, j, startI, endI, startJ, endJ ) \
for ( i = startI; i < endI; i ++ ) {	\
for ( j = startJ; j < endJ; j ++ ) {


#define FOR_LOOP3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( i = startI; i < endI; i ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( k = startK; k < endK; k ++ ) {     

#endif





#ifdef XFAST
#define Index3D( i, j, k, nx, ny, nz ) ( ( i ) + ( j ) * ( nx ) + ( k ) * ( ( nx ) * ( ny ) ) )

#define INDEX( i, j, k ) ( ( i ) + ( j ) * _nx_ + ( k ) * ( _nx_ * _ny_ ) )

#define INDEX_xi( i, j, k, offset ) ( ( i  offset ) + ( j ) * _nx_ + ( k ) * ( _nx_ * _ny_ ) )
#define INDEX_et( i, j, k, offset ) ( ( i ) + ( j  offset ) * _nx_ + ( k ) * ( _nx_ * _ny_ ) )
#define INDEX_zt( i, j, k, offset ) ( ( i ) + ( j ) * _nx_ + ( k  offset ) * ( _nx_ * _ny_ ) )

#endif

#ifdef ZFAST
#define Index3D( i, j, k, nx, ny, nz ) ( ( k ) + ( j ) * ( nz ) + ( i ) * ( ( nz ) * ( ny ) ) )

#define INDEX( i, j, k ) ( ( k ) + ( j ) * _nz_ + ( i ) * ( _nz_ * _ny_ ) )

#define INDEX_xi( i, j, k, offset ) ( ( k ) + ( j ) * _nz_ + ( i offset ) * ( _nz_ * _ny_ ) )
#define INDEX_et( i, j, k, offset ) ( ( k ) + ( j offset ) * _nz_ + ( i ) * ( _nz_ * _ny_ ) )
#define INDEX_zt( i, j, k, offset ) ( ( k offset ) + ( j ) * _nz_ + ( i ) * ( _nz_ * _ny_ ) )

#endif

/*
#define L( W, FB, SUB ) (  FB * ( af_1 * W[INDEX_##SUB( i, j, k, - FB * 1 )] + 				\
						 		   af0 * W[index] + 										\
						 		   af1 * W[INDEX_##SUB( i, j, k, + FB * 1 )] + 				\
						 		   af2 * W[INDEX_##SUB( i, j, k, + FB * 2 )] + 				\
						 		   af3 * W[INDEX_##SUB( i, j, k, + FB * 3 )] ) * rDH )
*/

#define L_J_T( J_T, FB ) ( FB * ( af_1 * J_T[3 - FB * 1] + \
					 			  af0  * J_T[3] + 		   \
					 			  af1  * J_T[3 + FB * 1] + \
					 			  af2  * J_T[3 + FB * 2] + \
					 			  af3  * J_T[3 + FB * 3] ) * rDH )

#define L2( W, VAR, VARSIZE, FB, SUB ) ( FB * ( W[INDEX_##SUB( i, j, k, + FB * 1 )* VARSIZE+VAR] - W[index* VARSIZE+VAR] ) * rDH )
#define L3( W, VAR, VARSIZE, FB, SUB ) ( FB * ( Cf1 * W[index* VARSIZE+VAR] + Cf2 * W[INDEX_##SUB( i, j, k, + FB * 1 )* VARSIZE+VAR] + Cf3 * W[INDEX_##SUB( i, j, k, + FB * 2 )* VARSIZE+VAR] ) * rDH ) 

#define DOT_PRODUCT3D( A1, A2, A3, B1, B2, B3 ) ( ( A1 ) * ( B1 ) + ( A2 ) * ( B2 ) + ( A3 ) * ( B3 ) )
#define DOT_PRODUCT2D( A1, A2, B1, B2 ) ( ( A1 ) * ( B1 ) + ( A2 ) * ( B2 ) )


#define DOT_PRODUCT3D_SVE( RET, A1, A2, A3, B1, B2, B3 ) \
Z1 	= svmul_m(pg,  A1, B1); \
Z2 	= svmul_m(pg,  A2, B2); \
Z3 	= svmul_m(pg,  A3, B3); \
RET = svadd_m(pg,  Z1, Z2); \
RET = svadd_m(pg, RET, Z3); 


#define DOT_PRODUCT2D_SVE( RET, A1, A2, B1, B2 ) \
Z1 	= svmul_m(pg, A1, B1); \
Z2 	= svmul_m(pg, A2, B2); \
RET = svadd_m(pg, Z1, Z2);

#define STORE_HW( HW, index, VARSIZE, VAR, TMP, A1, A2, A3, DT ) \
offset = svindex_s32(0, VARSIZE * sizeof(FLOAT));		  \
TMP = svadd_m(pg, A1,  A2);								  \
TMP = svadd_m(pg, TMP, A3);								  \
TMP = svmul_m(pg, TMP, DT);								  \
svst1_scatter_offset(pg, &(HW[index * VARSIZE + VAR]), offset, TMP);						





#define LC( C, VAR, VARSIZE, FB, SUB ) (  FB * ( af_1 * C[INDEX_##SUB( i, j, k, - FB * 1 ) * VARSIZE+VAR] + 		\
						 		   				  af0 * C[index * VARSIZE+VAR] +									\
						 		   				  af1 * C[INDEX_##SUB( i, j, k, + FB * 1 ) * VARSIZE+VAR] + 		\
						 		   				  af2 * C[INDEX_##SUB( i, j, k, + FB * 2 ) * VARSIZE+VAR] + 		\
						 		   				  af3 * C[INDEX_##SUB( i, j, k, + FB * 3 ) * VARSIZE+VAR] ) * rDH )


/*
#define L( W, VAR, VARSIZE, FB, SUB ) (  FB * ( af_1 * W[INDEX_##SUB( i, j, k, - FB * 1 ) * VARSIZE+VAR] + 		\
						 		   				 af0 * W[index * VARSIZE+VAR] +									\
						 		   				 af1 * W[INDEX_##SUB( i, j, k, + FB * 1 ) * VARSIZE+VAR] + 		\
						 		   				 af2 * W[INDEX_##SUB( i, j, k, + FB * 2 ) * VARSIZE+VAR] + 		\
						 		   				 af3 * W[INDEX_##SUB( i, j, k, + FB * 3 ) * VARSIZE+VAR] ) * rDH )

#define L( W, VAR, VARSIZE, FB, SUB ) (  FB * ( af_1   + 		\
						 		   				 af0   +		\
						 		   				 af1   + 		\
						 		   				 af2   + 		\
						 		   				 af3   ) * rDH )
*/
#define L( W, VAR, VARSIZE, FB, SUB ) (  FB * ( af_1 * W[INDEX_##SUB( i, j, k, - FB * 1 ) * VARSIZE+VAR] + 		\
						 		   				 af0 * W[index * VARSIZE+VAR] +									\
						 		   				 af1 * W[INDEX_##SUB( i, j, k, + FB * 1 ) * VARSIZE+VAR] + 		\
						 		   				 af2 * W[INDEX_##SUB( i, j, k, + FB * 2 ) * VARSIZE+VAR] + 		\
						 		   				 af3 * W[INDEX_##SUB( i, j, k, + FB * 3 ) * VARSIZE+VAR] ) * rDH )


#define L_SVE( RET, W, VAR, VARSIZE, FB, SUB )			\
offset = svindex_s32(0, VARSIZE * sizeof(FLOAT));		\
Z_1 = svld1_gather_offset( pg, &(W[INDEX_##SUB( i, j, k, - FB * 1) * VARSIZE + VAR]), offset );	\
Z0  = svld1_gather_offset( pg, &(W[index * VARSIZE + VAR]), offset );			   				\
Z1  = svld1_gather_offset( pg, &(W[INDEX_##SUB( i, j, k, + FB * 1) * VARSIZE + VAR]), offset );	\
Z2  = svld1_gather_offset( pg, &(W[INDEX_##SUB( i, j, k, + FB * 2) * VARSIZE + VAR]), offset );	\
Z3  = svld1_gather_offset( pg, &(W[INDEX_##SUB( i, j, k, + FB * 3) * VARSIZE + VAR]), offset );	\
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





#define LW( W, VAR, FB ) (  FB * ( af_1 * W[index_1+ VAR]  + \
						            af0 * W[index0 + VAR]  + \
						   		    af1 * W[index1 + VAR]  + \
						   		    af2 * W[index2 + VAR]  + \
						   		    af3 * W[index3 + VAR] ) )




#ifdef GPU_CUDA
#define END_CALCULATE1D( ) } 
#define END_CALCULATE2D( ) }
#define END_CALCULATE3D( ) }

#define CALCULATE1D( i, startI, endI )         \
if ( i >= startI && i < endI ) {



#define CALCULATE2D( i, j, startI, endI, startJ, endJ )         \
if ( i >= ( startI ) && i < ( endI ) && j >= ( startJ ) && j < ( endJ ) ) { 

#define CALCULATE3D( i, j, k, startI, endI, startJ, endJ, startK, endK )  \
if ( i >= ( startI ) && i < ( endI ) && j >= ( startJ ) && j < ( endJ ) && k >= ( startK ) && k < ( endK ) ) {

#else //GPU_CUDA

#define END_CALCULATE1D( ) } 
#define END_CALCULATE2D( ) }}
#define END_CALCULATE3D( ) }}}



#define CALCULATE1D( i, startI, endI )         \
    for ( i = ( startI ); i < endI; i ++ ) {


#define CALCULATE1D_SVE( i, len, startI, endI )         \
    for ( i = ( startI ); i < endI; i += len ) {


#define CALCULATE2D( i, j, startI, endI, startJ, endJ ) \
for ( j = startJ; j < endJ; j ++ ) {    \
for ( i = startI; i < endI; i ++ ) {    





#ifdef XFAST


#define CALCULATE3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( k = startK; k < endK; k ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( i = startI; i < endI; i ++ ) {     


#define CALCULATE3D_SVE( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( k = startK; k < endK; k ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( i = startI; i < endI; i += len ) {     



#endif //XFAST



#ifdef YFAST

#define CALCULATE3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( k = startK; k < endK; k ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( i = startI; i < endI; i ++ ) {    

#endif //YFAST


#ifdef ZFAST

#define CALCULATE3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( i = startI; i < endI; i ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( k = startK; k < endK; k ++ ) {     
#endif //ZFAST



#endif //GPU_CUDA




