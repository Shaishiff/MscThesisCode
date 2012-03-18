
#ifndef __DEFS_H_
#define __DEFS_H_

#define FILE_NAME_BUFFER_SIZE	1024

//////////////////////////////////////////////////////////////////////////////////////////
// Includes
//////////////////////////////////////////////////////////////////////////////////////////

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include "time.h"
#include <algorithm>
#include <vector>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
// MPI defines
//////////////////////////////////////////////////////////////////////////////////////////

const int MPI_MASTER = 0;

// Msg tags:
const int MPI_FLAG_MSG_TAG = 1000;
const int MPI_JOB_MSG_TAG = 2000;
const int MPI_RESULT_MSG_TAG = 3000;	

// Flags:
const int MPI_FLAG_QUIT = 1000;
const int MPI_FLAG_START_JOB = 2000;

//////////////////////////////////////////////////////////////////////////////////////////
// GA defines and consts
//////////////////////////////////////////////////////////////////////////////////////////

const int TargetCenterH = 5;
const int TargetCenterW = 5;
const int TargetHeight = 25;
const int TargetWidth = 25;

#define MAX_NUMBER_OF_THREADS 100
#define Npop 4//10
#define MaxIterations 1
#define SelectionRate 0.5//0.4
#define NsurvivingPopulation (int)floor(Npop*SelectionRate +0.5)
#define Nmates (Npop - NsurvivingPopulation)

//////////////////////////////////////////////////////////////////////////////////////////
// Physiological consts
//////////////////////////////////////////////////////////////////////////////////////////

const double Cm = 1.0; // membrance capacitance - microF/cm^2
const double Am = 3000.0; // 1/cm
const double sigma = 0.75; // 0.3; 3.0*10^-4; mS/cm

//////////////////////////////////////////////////////////////////////////////////////////
// APD & CV Restitution
//////////////////////////////////////////////////////////////////////////////////////////

const double t_fast  = 0.249; // tau_fast
const double tf_opn1 =  40; // tau_f_open_1
const double tf_opn2 =  82.5; // tau_f_open_2
const double tf_cls =  5.75; // tau_f_close
const double t_slow =226.9; // tau_slow
const double ts_opn =  100; // tau_s_open
const double ts_cls =  300; // tau_s_close
const double t_ug1 = 64.7; // tau_ung_1
const double t_ug2 = 222.9; // tau_ung_2
const double Vrest = -81.2; // in  mV
const double Vpeak = 3.6; // in  mV
const double V_f = 0.16;  
const double V_s = V_f; 
const double V_u = V_f;
const double Vf_open = 0.04;
const double V_slow  =  0.85;
const double kwm = 10;

//////////////////////////////////////////////////////////////////////////////////////////
// Simulation consts
//////////////////////////////////////////////////////////////////////////////////////////

// Time parameters:
const double dt = 0.02; // millisec
const double TotalSimulationTime = 200; // Total time of simulation
const int Nt = (int)(ceil(TotalSimulationTime/dt));

// Space parameters
#define dW 0.05 // cm/node
#define dH 0.05 // cm/node
#define W 2 // Width, x[cm]
#define H 2 // Height, y[cm]
const int Nw = (int)ceil(W/dW);
const int Nh = (int)ceil(H/dH);

struct ProtocolParams
{
	double m_Amp;
	double m_TotalTime;
	double m_BeginTime;
	double m_hStart;
	double m_hEnd;
	double m_wStart;
	double m_wEnd;
};

// Protocol parameters:
const double S1Amp = 500.0; // microA/cm^3
const double S1TotalTime = 5; // milliseconds
const double S1BeginTime = 5.0; // milliseconds
const double S1hStart = 0.0;
const double S1hEnd = dH;
const double S1wStart = 0.0;
const double S1wEnd = W;

struct S1Protocol : public ProtocolParams
{
	S1Protocol()
	{
		m_Amp = S1Amp;
		m_TotalTime = S1TotalTime;
		m_BeginTime = S1BeginTime;
		m_hStart = S1hStart;
		m_hEnd = S1hEnd;
		m_wStart = S1wStart;
		m_wEnd = S1wEnd;
	}
};

const double S2Amp = 500.0; // microA/cm^3
const double S2TotalTime = 5; // milliseconds
const double S2BeginTime = 5.0; // milliseconds
const double S2hStart = 0.0;
const double S2hEnd = H;
const double S2wStart = 0.0;
const double S2wEnd = dW;

struct S2Protocol : public ProtocolParams
{
	S2Protocol()
	{
		m_Amp = S2Amp;
		m_TotalTime = S2TotalTime;
		m_BeginTime = S2BeginTime;
		m_hStart = S2hStart;
		m_hEnd = S2hEnd;
		m_wStart = S2wStart;
		m_wEnd = S2wEnd;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __DEFS_H_

