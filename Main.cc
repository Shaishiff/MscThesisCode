
// Desktop Version
// Symmetric propagation with central stimulation

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Nhumatr.h" 

const double PI = 3.1415926;

// Saving parameters
//#define out_file_name "voltage.%d"
//#define out_file_name_current "current.%d"
//const int Save_W = 100; //210; // 42; // width of the saved frame
//const int Save_H = 100; //560; // 112; // hight of the saved frame
//const int Save_dW = 1; // 5; // width resolution of saved elements (Save_W*Save_dW = W MUST !!) 
//const int Save_dH = 1; // 5; // height resolution of saved elements (Save_H*Save_dH = H MUST !!)

// Time parameters
const double MaxSimulationTime = 400.0; // Total time of simulation
const int Nt = (int)(ceil(MaxSimulationTime/dt));

//const double save_period = 1.0; // milliseconds
//const double send_period = 1.0; // milliseconds - For collecting the data from Workers to MASTER for saving - Make equal to save_period !!!
//const double time_step = (double) dt; // milliseconds - Defind in const.h that is included in Nhumatr.h
//const double max_time = 500.0; // milliseconds - Total time of simulation
//double time_to_save = 0.0; // Initialization of saving
//double time_to_send = 0.0; // Initialization of sending
//double sim_time = 0.0; // Used locally/independently in MASTER and Workers

// Space parameters
//const double Width = 10.0; //10.5; // mm
//const double Height = 10.0; //28.0; // mm
//const double dW = 0.1; // mm/node
//const double dH = 0.1; // mm/node
//const int W = 142; //210; // = Width/dW; // Width, x
//const int H = 142; //560; // = Height/dH; // Height, y

#define dW 0.01 // cm/node
#define dH 0.01 // cm/node
const double dH2 = dH*dH;
const double dW2 = dW*dW;
#define W 1.4 // Width, x[cm]
#define H 1.4 // Height, y[cm]
const int Nw = (int)ceil(W/dW); // In indexes
const int Nh = (int)ceil(H/dH); // In indexes
#define Nw_with_border (Nw+2)
#define Nh_with_border (Nh+2)
#define MeasurementMarginIndexes 20
#include "Protocols.h"

#define INVALID_RISE_TIME	1000.0
double Diff = 0.04*1.4; // was 0.013

/*
// Reaction - Diffusion details
// const double Cm = (double) 100.0e-12;
const double Cm_factor = (double) 1.0/(Cm); 

// Difusion  parameters
struct D_tensor_type
{
	D_tensor_type()
	{
		Dii = 0.0;
		Djj = 0.0;
		Dij = 0.0;
	}
	double Dii, Djj, Dij;
};
*/

// Fibrosis parameters

/*
const double Rm_fibro_low = (double) 5.0e+09; // Ohm Kamkin A. Experimental Physiology 1999
const double Rm_fibro_high = (double) 25.0e+09; // Ohm Kamkin A. Experimental Physiology 1999
const double Cm_fibro_factor_low = (double) 1.0/((double) Rm_fibro_low*Cm);
const double Cm_fibro_factor_high = (double) 1.0/((double) Rm_fibro_high*Cm);
const double Vrm_fibro = (double) -15.9; // mV Kamkin A. Experimental Physiology 1999
*/

// Stimulation parameters
//const double S1_amp = -142.4;//3.0752 -100.0;//-100*Cm; // pA
/*
const double S1_amp = -100.0;//3.0752 -100.0;//-100*Cm; // pA
const double S1_total_time = 5.0; // 50.0 milliseconds
const double S1_begin = 2.0; // milliseconds
*/

/*
const double S2_amp = 0.0;//-100.0;//-3000.0;//-3000.0; //-80.0*Cm; // pA
const double S2_total_time = 4.0; // milliseconds
const double S2_begin = 200.0; // milliseconds
*/

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double** inFibroblastMat = NULL;
double** outRiseTimeMat = NULL;
double dVdh = 0.0;
double dVdw = 0.0;

const int Protocol = 1; // Can be either 1 or 2
const int m_nHStart = 50;
const int m_nWStart = 40;
const int m_nHEnd = 64;
const int m_nWEnd = 64;

state_variables** Node = NULL;
state_variables** pNewNode = NULL;
state_variables** pTempNode = NULL;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void ShowOutput()
{	
	for (int iH = 1; iH < Nh+1; iH += 4)
	{
		for (int iW = 1; iW < Nw+1; iW += 4)
		{
			if (inFibroblastMat[iH][iW] != 0.0)	printf("  ");
			else if (Node[iH][iW].V < -80.0) printf("88");
			else if (Node[iH][iW].V < -70.0) printf("77");
			else if (Node[iH][iW].V < -60.0) printf("66");
			else if (Node[iH][iW].V < -50.0) printf("55");
			else if (Node[iH][iW].V < -40.0) printf("44");
			else if (Node[iH][iW].V < -30.0) printf("33");
			else if (Node[iH][iW].V < -20.0) printf("22");
			else if (Node[iH][iW].V < -10.0) printf("11");
			else if (Node[iH][iW].V < 0.0) printf("00");
			else if (Node[iH][iW].V < 5.0) printf("AA");
			else if (Node[iH][iW].V < 10.0) printf("BB");
			else if (Node[iH][iW].V < 15.0) printf("CC");
			else if (Node[iH][iW].V < 20.0) printf("DD");
			else if (Node[iH][iW].V < 25.0) printf("++");
			else if (Node[iH][iW].V < 25.0) printf("**");
			else printf("^^");
			
		}
		printf("\n");
	}
	printf("\n"); // Separator
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void SaveToFile(double sim_time)
{
	double sim_time_to_save = floor(sim_time + 0.5);
	char fileName[1024] = {0};
	sprintf(fileName, "Output\\ModelOutput_at_%.2f.txt", sim_time_to_save);
	FILE* pFile = fopen(fileName, "w");
	if(pFile == NULL)
	{
		return;
	}

	for (int iH = 1; iH < Nh+1; ++iH)
	{
		for (int iW = 1; iW < Nw+1; ++iW)
		{
			fprintf(pFile, "%4.3f ", Node[iH][iW].V);
		}
		fprintf(pFile, "\n");
	}	
	fclose(pFile);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double** CreateMat()
{
	double* data = (double *)malloc(Nh_with_border*Nw_with_border*sizeof(double));
	double** mat = (double **)malloc(Nh_with_border*sizeof(double*));
	for (int i = 0; i < Nh_with_border; i++)
	{
		mat[i] = &(data[Nw_with_border*i]);
	}

	for (int iH = 0; iH < Nh_with_border; ++iH)
	{
		for (int iW = 0; iW < Nw_with_border; ++iW)
		{
			mat[iH][iW] = 0.0;
		}
	}

	return mat;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

state_variables** CreateNode()
{
	state_variables* data = (state_variables *)malloc(Nh_with_border*Nw_with_border*sizeof(state_variables));
	state_variables** mat = (state_variables **)malloc(Nh_with_border*sizeof(state_variables*));
	for (int i = 0; i < Nh_with_border; i++)
	{
		mat[i] = &(data[Nw_with_border*i]);
	}
	
	return mat;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CreateFibroblastBorders()
{
	for (int iH = 0; iH < Nh_with_border; ++iH)
	{
		for (int iW = 0; iW < Nw_with_border; ++iW)
		{								
			if((iH == 0) || (iH == Nh+1) || (iW == 0) || (iW == Nw+1))
			{
				inFibroblastMat[iH][iW] = 1.0;
			}
			else
			{
				inFibroblastMat[iH][iW] = 0.0;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CalculateDer(int iH, int iW, state_variables** Node)
{
	if(inFibroblastMat[(iH-1)][iW] == 1.0)
	{
		dVdh = (Node[(iH+1)][iW].V - Node[iH][iW].V)/dH2;
	}
	else if (inFibroblastMat[(iH + 1)][iW] == 1.0)
	{
		dVdh = (Node[(iH-1)][iW].V - Node[iH][iW].V)/dH2;
	}
	else
	{
		dVdh = (Node[(iH-1)][iW].V + Node[(iH+1)][iW].V - 2.0*Node[iH][iW].V)/dH2;
	}

	if(inFibroblastMat[iH][iW-1] == 1)
	{
		dVdw = (Node[iH][iW+1].V - Node[iH][iW].V)/dW2;
	}
	else if (inFibroblastMat[iH][iW+1] == 1)
	{
		dVdw  = (Node[iH][iW-1].V - Node[iH][iW].V)/dW2;
	}
	else
	{
		dVdw  = (Node[iH][iW-1].V + Node[iH][iW+1].V - 2.0*Node[iH][iW].V)/dW2;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CreateFibroblastPatch()
{
	for (int iH = m_nHStart; iH <= m_nHEnd; ++iH)
	{	
		for (int iW = m_nWStart; iW <= m_nWEnd; ++iW)
		{
			inFibroblastMat[iH][iW] = 1.0;
		}
	}	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void ClearRiseTimeMat()
{
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		for (int iW = 1; iW < Nw+1; ++iW)
		{            
			if(inFibroblastMat[iH][iW] != 1)
			{
				outRiseTimeMat[iH][iW] = INVALID_RISE_TIME;
			}
			else
			{
				outRiseTimeMat[iH][iW] = -1.0;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

ProtocolParams* ChooseProtocol()
{
	if(Protocol == 1)
	{
		return new S1Protocol();
	}
	else if(Protocol == 2)
	{
		return new S2Protocol();
	}
	return NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void InitNodes()
{
	Get_state_variables_initial_condition(); // Put initial conditions in the vector state[...]
	for (int iH = 0; iH < Nh_with_border; ++iH)
	{
		for (int iW = 0; iW < Nw_with_border; ++iW)
		{
			Assign_node_initial_condition(Node[iH][iW]);
			Assign_node_initial_condition(pNewNode[iH][iW]);		
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void WriteFibroblastMatFiles()
{
	FILE* pFile = fopen("Output\\TargetFibroblastMatReadable.txt", "w");
	if(pFile != NULL)
	{		
		for (int iH = 1; iH < Nh+1; ++iH)
		{
			for (int iW = 1; iW < Nw+1; ++iW)
			{
				fprintf(pFile, "%d", (int)ceil(inFibroblastMat[iH][iW]));
			}
			fprintf(pFile, "\n");
		}	
		fclose(pFile);
	}

	pFile = fopen("Output\\TargetFibroblastMat.txt", "w");
	if(pFile != NULL)
	{		
		for (int iH = 1; iH < Nh+1; ++iH)
		{
			for (int iW = 1; iW < Nw+1; ++iW)
			{
				fprintf(pFile, "%.3f ", inFibroblastMat[iH][iW]);
			}
			fprintf(pFile, "\n");
		}	
		fclose(pFile);
		pFile = NULL;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void WriteRiseTimeFile()
{
	char targetFilename[1024] = {0};
	sprintf(targetFilename, "Output\\TargetFibroblastMatResults%d.txt", Protocol);
	FILE* pFibroBlastFile = fopen(targetFilename, "w");
	if(pFibroBlastFile != NULL)
	{		
		for (int iH = 1; iH < Nh+1; ++iH)
		{
			for (int iW = 1; iW < Nw+1; ++iW)
			{
				fprintf(pFibroBlastFile, "%.3f ", outRiseTimeMat[iH][iW]);
			}
			fprintf(pFibroBlastFile, "\n");
		}	
		fclose(pFibroBlastFile);
		pFibroBlastFile = NULL;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	ProtocolParams* pProt = ChooseProtocol();
	
	inFibroblastMat = CreateMat();	
	CreateFibroblastBorders();
	CreateFibroblastPatch();
	WriteFibroblastMatFiles();
		
	outRiseTimeMat = CreateMat();
	ClearRiseTimeMat();
	double dFirstRiseTime = 0.0;
	
	Node = CreateNode();
	pNewNode = CreateNode();
	InitNodes();

	// Start the temporal loop.
	for (int iT = 0; iT < Nt; ++iT)
	{
		double curTime = iT*dt;
		bool bS1 = ((curTime >= pProt->m_BeginTime) && (curTime <= (pProt->m_BeginTime+pProt->m_TotalTime)));

		if(iT%20 == 0 && iT != 0)
		{
			printf("%.2f ... ", curTime);
		}
		if(iT%200 == 0)
		{
			printf("Mat at time: %.2f\n", curTime);
			SaveToFile(curTime);
			ShowOutput();
			printf("Time: ", curTime);
		}

		for (int iW = 1; iW < Nw+1; ++iW)		
		{
			for (int iH = 1; iH < Nh+1; ++iH)
			{            
				if(inFibroblastMat[iH][iW] != 1.0)
				{	
					// Check if we need to add a stimulation current
					// according to the protocol.
					double Jstim = 0.0;
					if(bS1)
					{
						if(iH >= pProt->m_hStart &&
						   iH <= pProt->m_hEnd &&
						   iW >= pProt->m_wStart &&
						   iW <= pProt->m_wEnd)
						{
							Jstim = S1Amp*(-100.0);
						}
					}
				
					if((outRiseTimeMat[iH][iW] == INVALID_RISE_TIME) && (Node[iH][iW].V > 0.0))
					{
						if(dFirstRiseTime == 0.0)
						{
							dFirstRiseTime = iT*dt;
						}
						outRiseTimeMat[iH][iW] = iT*dt - dFirstRiseTime;
					}

					double I = Total_transmembrane_currents(Node[iH][iW], Jstim);
					CalculateDer(iH, iW, Node);
					pNewNode[iH][iW].V = Node[iH][iW].V + dt*(Diff*(dVdh + dVdw) - I);					
				}			
			}			
		}
		
		pTempNode = Node;
		Node = pNewNode;
		pNewNode = Node;

		// Check if we can stop the simulation.
		bool bStopSim = true;
		for (double curH = pProt->m_hMeasureStart; curH <= pProt->m_hMeasureEnd; curH += dH)
		{
			for (double curW = pProt->m_wMeasureStart; curW <= pProt->m_wMeasureEnd; curW += dW)
			{            									
				int iH = (int)ceil(curH/dH);
				int iW = (int)ceil(curW/dW);
				if(outRiseTimeMat[iH][iW] == INVALID_RISE_TIME)
				{
					bStopSim = false;
					break;
				}
			}
		}

		if(bStopSim)
		{
			return 0;
		}
	} // sim_time loop

	WriteRiseTimeFile();
	
	return 0;
} // of main

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
