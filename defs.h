
#ifndef __DEFS_H_
#define __DEFS_H_

#define FILE_NAME_BUFFER_SIZE	1024

//////////////////////////////////////////////////////////////////////////////////////////
// Includes
//////////////////////////////////////////////////////////////////////////////////////////

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
const int MPI_DIE_TAG = 1000;
const int MPI_JOB_1_TAG = 2000;
const int MPI_JOB_2_TAG = 3000;
const int MPI_RESULT_TAG = 4000;

//////////////////////////////////////////////////////////////////////////////////////////
// GA defines and consts
//////////////////////////////////////////////////////////////////////////////////////////

#define Npop 15
#define MaxIterations 100
#define SelectionRate 0.4
#define NsurvivingPopulation (int)floor(Npop*SelectionRate +0.5)
#define Nmates (Npop - NsurvivingPopulation)

//////////////////////////////////////////////////////////////////////////////////////////
// Physiological consts
//////////////////////////////////////////////////////////////////////////////////////////

const double Cm = 1.0; // membrance capacitance - microF/cm^2
const double Am = 3000.0; // 1/cm
const double sigma = 13.865; // 0.3; 3.0*10^-4; mS/cm

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

// Space parameters
#define dW 0.01 // cm/node
#define dH 0.01 // cm/node
#define W 1 // Width, x[cm]
#define H 1 // Height, y[cm]
const int Nw = (int)ceil(W/dW);
const int Nh = (int)ceil(H/dH);
#define Nw_with_border (Nw+2)
#define Nh_with_border (Nh+2)

// The values here are in cell indexes:
const int TargetCenterH = 40;
const int TargetCenterW = 40;
const int TargetHeight = 20;
const int TargetWidth = 20;

// Time parameters:
const double dt = 0.005; // millisec
const double TotalSimulationTime = max(W,H)*400; // Total time of simulation
const int Nt = (int)(ceil(TotalSimulationTime/dt));

double** CreateMat();
void DestroyMat(double** in);
void PrintMat(double** mat);
bool SaveMatToFile(double** mat, char* fileName);

//////////////////////////////////////////////////////////////////////////////////////////
// Protocols
//////////////////////////////////////////////////////////////////////////////////////////

struct ProtocolParams
{
	double m_Amp;
	double m_TotalTime;
	double m_BeginTime;
	double m_hStart;
	double m_hEnd;
	double m_wStart;
	double m_wEnd;
	double m_hMeasureStart;
	double m_hMeasureEnd;
	double m_wMeasureStart;
	double m_wMeasureEnd;
};

// Protocol parameters:
const double S1Amp = 5000.0;//500.0; // microA/cm^3
const double S1TotalTime = 5; // milliseconds
const double S1BeginTime = 2.0; // milliseconds
const double S1hStart = 0.0;
const double S1hEnd = dH*2;
const double S1wStart = 0.0;
const double S1wEnd = W;
const double S1hMeasureStart = H;
const double S1hMeasureEnd = H;
const double S1wMeasureStart = 0.0;
const double S1wMeasureEnd = W;

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
		m_hMeasureStart = S1hMeasureStart;
		m_hMeasureEnd = S1hMeasureEnd;
		m_wMeasureStart = S1wMeasureStart;
		m_wMeasureEnd = S1wMeasureEnd;
	}
};

const double S2Amp = 5000.0;//500.0; // microA/cm^3
const double S2TotalTime = 5; // milliseconds
const double S2BeginTime = 2.0; // milliseconds
const double S2hStart = 0.0;
const double S2hEnd = H;
const double S2wStart = 0.0;
const double S2wEnd = dW*2;
const double S2hMeasureStart = 0.0;
const double S2hMeasureEnd = H;
const double S2wMeasureStart = W;
const double S2wMeasureEnd = W;

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
		m_hMeasureStart = S2hMeasureStart;
		m_hMeasureEnd = S2hMeasureEnd;
		m_wMeasureStart = S2wMeasureStart;
		m_wMeasureEnd = S2wMeasureEnd;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __DEFS_H_

