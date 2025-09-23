/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:init_grid.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-06
*   Discription:
*
================================================================*/
#include "header.h"

void init_heter_grid( PARAMS params, GRID * grid, MPI_COORD thisMPICoord )
{
	int resX = 0;	
	int resY = 0;
	int resZ = 0;	
	
	grid->PX = params.PX;
	grid->PY = params.PY;
	grid->PZ = params.PZ;

	grid->_NX_ = params.NX + 2 * HALO;
	grid->_NY_ = params.NY + 2 * HALO;
	grid->_NZ_ = params.NZ + 2 * HALO;
	
	grid->_NX = params.NX + HALO;
	grid->_NY = params.NY + HALO;
	grid->_NZ = params.NZ + HALO;

	grid->NX = params.NX;
	grid->NY = params.NY;
	grid->NZ = params.NZ;
	
	if ( grid->PX >= 3)
	{
		if ( (thisMPICoord.X == 0) || (thisMPICoord.X == grid->PX - 1) )
		{
			if (thisMPICoord.X == 0)
			{
				grid->nx = params.nxl;
				grid->frontNX = 0;
			}
			if (thisMPICoord.X == grid->PX - 1)
			{
				grid->nx = params.nxr;
				grid->frontNX = grid->NX - grid->nx;
			}
		}	
		else
		{	
			int ResNX = params.NX - params.nxl - params.nxr;
			grid->nx = ResNX / (params.PX-2);
			resX = ResNX % (params.PX-2);
			if ( thisMPICoord.X - 1 < resX )
			{
				grid->nx ++;
				grid->frontNX = (thisMPICoord.X-1) * grid->nx + params.nxl;
			}	
			else
			{	
				grid->frontNX = resX * ( grid->nx + 1 ) + ( thisMPICoord.X - 1 - resX ) * grid->nx + params.nxl;
			}
		}
	}
	
	if ( grid->PY >= 3)
	{
		if ( (thisMPICoord.Y == 0) || (thisMPICoord.Y == grid->PY - 1) )
		{
			if (thisMPICoord.Y == 0)
			{
				grid->ny = params.nyl;
				grid->frontNY = 0;
			}
			if (thisMPICoord.Y == grid->PY - 1)
			{
				grid->ny = params.nyr;
				grid->frontNY = grid->NY - grid->ny;
			}
		}	
		else
		{
			int ResNY = params.NY - params.nyl - params.nyr;
			grid->ny = ResNY / (params.PY - 2);
			resY = ResNY % (params.PY - 2);
			if ( thisMPICoord.Y - 1 < resY )
			{
				grid->ny ++;
				grid->frontNY = (thisMPICoord.Y-1) * grid->ny + params.nyl;
			}	
			else
			{	
				grid->frontNY = resY * ( grid->ny + 1 ) + ( thisMPICoord.Y - 1 - resY ) * grid->ny + params.nyl;
			}
		}
	}
	
	if ( grid->PZ >= 3)
	{
		if ( (thisMPICoord.Z == 0) || (thisMPICoord.Z == grid->PZ - 1) )
		{
			if (thisMPICoord.Z == 0)
			{
				grid->nz = params.nzl;
				grid->frontNZ = 0;
			}
			if (thisMPICoord.Z == grid->PZ - 1)
			{
				grid->nz = params.nzr;
				grid->frontNZ = grid->NZ - grid->nz;
			}
		}
		else
		{
			int ResNZ = params.NZ - params.nzl - params.nzr;
			grid->nz = ResNZ / (params.PZ - 2);
			resZ = ResNZ % (params.PZ - 2);
			if ( thisMPICoord.Z - 1 < resZ )
			{
				grid->nz ++;
				grid->frontNZ = (thisMPICoord.Z - 1) * grid->nz + params.nzl;
			}	
			else
			{	
				grid->frontNZ = resZ * ( grid->nz + 1 ) + ( thisMPICoord.Z - 1 - resZ ) * grid->nz + params.nzl;
			}
		}
	}

	MPI_Barrier( MPI_COMM_WORLD );

	grid->_frontNX = grid->frontNX + HALO;
	grid->_frontNY = grid->frontNY + HALO;
	grid->_frontNZ = grid->frontNZ + HALO;


	grid->_nx = grid->nx + HALO;
	grid->_ny = grid->ny + HALO;
	grid->_nz = grid->nz + HALO;

	grid->_nx_ = grid->nx + 2 * HALO;
	grid->_ny_ = grid->ny + 2 * HALO;
	grid->_nz_ = grid->nz + 2 * HALO;
	
	grid->originalX = params.centerX;
	grid->originalY = params.centerY;

	grid->_originalX = grid->originalX + HALO;
	grid->_originalY = grid->originalY + HALO;

	grid->halo = HALO;

	grid->DH = params.DH;
	grid->rDH = 1.0 / grid->DH;
	//printf( "dh = %e\n", grid->DH );
	

	grid->nPML = params.nPML;
}
