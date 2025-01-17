/**************************************************
**** All rights reserved.
**** Author: Wenqiang Wang
**** Department of Earth and Space Sciences, Sustech
**** E-mail: 11849528@mail.sustech.edu.cn
**** Date:		2019-03-27
**** Update:	2019-03-31
**** update:		2019-04-11
**************************************************/
#ifndef HEADER_H
#define HEADER_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


#ifdef GPU_CUDA
#include <cublas_v2.h>
#include <cuda_fp16.h>
//#include <cuda_bf16.h>
#endif

#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
#include <list>
#include <map>
#include <vector>


#if __GNUC__
#include <sys/stat.h>
#include <sys/types.h>
#elif _MSC_VER
#include <windows.h>
#include <direct.h>
#endif


#include <proj.h>
#include "macro.h"
#include "struct.h"
#include "cJSON.h"


#include "functions.h"


#endif // !HEADER_H

