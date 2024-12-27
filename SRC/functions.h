/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:functions.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-05
*   Discription:
*
================================================================*/
#ifndef __FUNCTIONS__
#define __FUNCTIONS__

/*extern*/ void getParams( PARAMS * params);

void projTrans( double lon_0, double lat_0, GRID grid, float * coord, LONLAT LonLat, float lonMin, float lonMax, float latMin, float latMax );
double interp2d(double x[2], double y[2], double z[4], double x_, double y_ );
/*
double bilinear( double x, double y, double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22 ); 
double bilinearInterp( double x, double y, double x1, double y1, double x2, double y2, double A, double B, double C, double D );
void cart2LonLat( GRID grid, PJ * P,  PJ_DIRECTION PD, COORD coord, LONLAT LonLat );
void preprocessTerrain( PARAMS params, MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, COORD coord );
*/
void init_grid( PARAMS params, GRID * grid, MPI_COORD thisMPICoord );
void printInfo( GRID grid, PARAMS params );
void createDir( PARAMS params );
void modelChecking( PARAMS params  );



void allocCoord( GRID grid, float ** coord, float ** cpu_coord );
void freeCoord( float * coord, float * cpu_coord );


void allocTerrain( GRID grid, float ** terrain, float ** cpu_terrain );
void freeTerrain( float * terrain, float * cpu_terrain );

void constructCoord(MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, PARAMS params, float * coord, float * cpu_coord, float * terrain, float * cpu_terrain );
void preprocessTerrain( PARAMS params, MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, float * cpu_coord, float * cpu_terrain );



void allocMedium( GRID grid, float ** medium, float ** cpu_medium );
void freeMedium( float * medium, float * cpu_medium );


void allocCJM( GRID grid, FLOAT ** CJM );
void freeCJM( FLOAT * CJM);
void  constructCJM( MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, GRID grid, FLOAT *CJM, float *cpu_coord, float *medium, Mat_rDZ mat_rDZ );


void allocMat_rDZ( GRID grid, Mat_rDZ * mat_rDZ );
void freeMat_rDZ( Mat_rDZ mat_rDZ );


void constructMedium( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * cpu_coord, float * cpu_terrain, float * medium, float * cpu_medium );
void readWeisenShenModel( PARAMS params, GRID grid, MPI_COORD thisMPICoord, float * cpu_coord, float * cpu_terrain, float * structure );
void readCrustal_1( PARAMS params, GRID grid, MPI_COORD thisMPICoord, float * coord, float * structure );

void calc_CFL( GRID grid, float * cpu_coord, float * cpu_medium, PARAMS params );

void init_MultiSource( PARAMS params, GRID grid, MPI_COORD thisMPICoord, float * coord, float * terrain, long long ** srcIndex, float ** momentRate, float ** momentRateSlice, SOURCE_FILE_INPUT * ret_src_in );


void addMomenteRate( GRID grid, SOURCE_FILE_INPUT src_in, 
		FLOAT * hW, long long * srcIndex, 
		float * momentRate, float * momentRateSlice, 
		int it, int irk, float DT, float DH, 
		float * gaussFactor, int nGauss, int flagSurf, int * useMultiSource );


void finish_MultiSource( long long * srcIndex, float * momentRate, float * momentRateSlice, long long npts );

void calculateMomentRate( SOURCE_FILE_INPUT src_in, FLOAT * CJM, float * momentRate, long long * srcIndex, float DH );

void allocGaussFactor( float ** gaussFactor, int nGauss );
void gaussSmooth( float * gaussFactor, int nGauss );
void freeGaussFactor( float * gaussFactor  );

/*
#ifdef FREE_SURFACE
void solveContravariantJac( MPI_Comm comm_cart, GRID grid, float * con, float * coord, float * Jac, 
float * medium, Mat_rDZ mat_rDZ );
#else
void solveContravariantJac( MPI_Comm comm_cart, GRID grid, float * con, float * coord, float * Jac );
#endif
*/
  

void allocWave( GRID grid, WAVE * wave );
void freeWave( WAVE wave );

void waveDeriv( GRID grid, WAVE wave, FLOAT * CJM, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int FB1, int FB2, int FB3, float DT );	

void freeSurfaceDeriv( 
	GRID grid, WAVE wave, FLOAT * CJM, Mat_rDZ mat_rDZ, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int FB1, int FB2, int FB3, float DT );							



void waveRk( GRID grid, int irk, WAVE wave );

void allocPMLParameter( GRID grid, PML_ALPHA * pml_alpha, PML_BETA *pml_beta, PML_D * pml_d  );

void freePMLParamter( PML_ALPHA pml_alpha, PML_BETA pml_beta, PML_D pml_d );

void init_pml_parameter( PARAMS params, GRID grid, MPI_BORDER border, PML_ALPHA pml_alpha, PML_BETA pml_beta, PML_D pml_d );


void allocPML( GRID grid, AUX6SURF * Aux6, MPI_BORDER border );
void freePML( MPI_BORDER border, AUX6SURF Aux6  );

void pmlDeriv( GRID grid, WAVE wave, FLOAT * CJM, AUX6 Aux6, PML_ALPHA pml_alpha, 
PML_BETA pml_beta, PML_D pml_d, MPI_BORDER border, int FB1, int FB2, int FB3, float DT );
void pmlRk( GRID grid, MPI_BORDER border, int irk, AUX6 Aux6 );

void pmlFreeSurfaceDeriv( GRID grid, WAVE wave, 
						  FLOAT * CJM, AUX6 Aux6, Mat_rDZ mat_rDZ,
						  PML_D pml_d, MPI_BORDER border, int FB1, int FB2, float DT );


void finalize_MPI( MPI_Comm * comm_cart );
void init_MPI( int *argc, char *** argv, PARAMS params, MPI_Comm * comm_cart, MPI_COORD * thisMPICoord, MPI_NEIGHBOR * mpiNeigbor );

void init_gpu( int PX, int PY, int PZ  );
void run( MPI_Comm comm_cart, MPI_COORD thisMPICoord, MPI_NEIGHBOR mpiNeighbor, GRID grid, PARAMS params );
void propagate( 
MPI_Comm comm_cart, MPI_COORD thisMPICoord, MPI_NEIGHBOR mpiNeighbor,
GRID grid, PARAMS params, 
WAVE wave, Mat_rDZ mat_rDZ, FLOAT * CJM,
SOURCE_FILE_INPUT src_in, long long * srcIndex, float * momentRate, float * momentRateSlice, 
SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu );



void allocStation( STATION * station, STATION * cpu_station, int stationNum, int NT );
void freeStation( STATION station, STATION cpu_station );

int readStationIndex( GRID grid ) ;
void initStationIndex( GRID grid, STATION station );
void stationCPU2GPU( STATION station, STATION station_cpu, int stationNum );
void storageStation( GRID grid, int NT, int stationNum, STATION station, FLOAT * W, int it );
void stationGPU2CPU( STATION station, STATION station_cpu, int stationNum, int NT );
void stationWrite( PARAMS params, GRID grid, MPI_COORD thisMPICoord, STATION station, int NT, int stationNum );


void loadPointSource( GRID grid, SOURCE S, FLOAT * h_W, FLOAT * CJM, int it, int irk, float DT, float DH, float rickerfc );

void GaussField( GRID grid, FLOAT * W );



void locateSlice( PARAMS params, GRID grid, SLICE * slice );

void locateSource( PARAMS params, GRID grid, SOURCE * source );
void locateFreeSurfSlice( GRID grid, SLICE * slice );

void allocSliceData( GRID grid, SLICE slice, SLICE_DATA * sliceData, SLICE_DATA * sliceDataCpu );
void freeSliceData( GRID grid, SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu  );
void data2D_output_bin( GRID grid, SLICE slice, 
						MPI_COORD thisMPICoord, 
						float * datain, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu,	
						const char * name );
void data2D_XYZ_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, FLOAT * wave, SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu, char var, int it );
void data2D_Model_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, float * coord, float * medium, SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu );





void allocSendRecv( GRID grid, MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA * sr, int VARSIZE );
void freeSendRecv( MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA sr );
void mpiSendRecv( MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, GRID grid, float * con, SEND_RECV_DATA sr, int VARSIZE );

void FLOAT_allocSendRecv( GRID grid, MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA_FLOAT * sr, int VARSIZE );
void FLOAT_freeSendRecv( MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA_FLOAT sr );
void FLOAT_mpiSendRecv( MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, GRID grid, FLOAT * wave, SEND_RECV_DATA_FLOAT sr, int VARSIZE );



void allocatePGV( GRID grid, float ** pgv, float ** cpu_pgv );
void freePGV( float * pgv, float * cpu_pgv );
void outputPGV( PARAMS params, GRID grid, MPI_COORDINATE thisMPICoord, float * pgv, float * cpuPgv );
void comparePGV( GRID grid, MPI_COORDINATE thisMPICoord, FLOAT * W, float * pgv );

#endif //__FUNCTIONS__
