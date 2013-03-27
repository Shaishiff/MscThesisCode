
#ifndef __DEFS_H_
#define __DEFS_H_

#define FILE_NAME_BUFFER_SIZE	1024

//////////////////////////////////////////////////////////////////////////////////////////
// Includes
//////////////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "math.h"
#include "time.h"
using namespace std;

#define EPSILON	double(0.000001)

//////////////////////////////////////////////////////////////////////////////////////////
// MPI defines
//////////////////////////////////////////////////////////////////////////////////////////

const int MPI_MASTER = 0;

// Msg tags:
const int MPI_DIE_TAG = 1000;
const int MPI_JOB_1_TAG = 2000;
const int MPI_JOB_2_TAG = 3000;
const int MPI_RESULT_TAG = 4000;

//////////////////////////////////////////////////////////////////////////////////////////
// GA defines and consts
//////////////////////////////////////////////////////////////////////////////////////////

#define Npop 20
#define MaxIterations 5
#define SelectionRate 0.4
#define NsurvivingPopulation (int)floor(Npop*SelectionRate + 0.5)
#define Nmates (Npop - NsurvivingPopulation)
#define SAMPLING_INTERVALS	1
#define MAX_REPEATING_COSTS_FOR_DEAD_END	20

//////////////////////////////////////////////////////////////////////////////////////////
// Simulation consts
//////////////////////////////////////////////////////////////////////////////////////////

// Space parameters
#define dW 0.01 // cm/node
#define dH 0.01 // cm/node
#define W 1.4 // Width, x[cm]
#define H 1.4 // Height, y[cm]
const int Nw = (int)ceil(W/dW); // In indexes
const int Nh = (int)ceil(H/dH); // In indexes
#define Nw_with_border (Nw+2)
#define Nh_with_border (Nh+2)
#define MeasurementMarginIndexes 20
#define Min_w_Fibroblast (MeasurementMarginIndexes + 1)
#define Min_h_Fibroblast (MeasurementMarginIndexes + 1)
#define Max_w_Fibroblast (Nw - (MeasurementMarginIndexes + 1))
#define Max_h_Fibroblast (Nh - (MeasurementMarginIndexes + 1))
#define nHPartitionSize (ceil((Nh - MeasurementMarginIndexes)/4) + 5)
#define nWPartitionSize (ceil((Nw - MeasurementMarginIndexes)/4) + 5)

// Time parameters:
const double dt = 0.005; // millisec
const double MaxSimulationTime = 50.0; // Total time of simulation
const int Nt = (int)(ceil(MaxSimulationTime/dt));

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __DEFS_H_
