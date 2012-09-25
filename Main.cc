
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "Nhumatr.h" 

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// Time parameters
#define SAVE_OUTPUT_PERIOD 1.0 // milliseconds - For collecting the data from Workers to MASTER for saving - Make equal to save_period !!!
#define MAX_SIMULATION_TIME 450.0 // milliseconds - Total time of simulation
#define INVALID_RISE_TIME 1000.0
#define FIBROBLAST_RISE_TIME (-1.0)

// Fibroblasts
#define FIBROBLAST_H_START 0
#define FIBROBLAST_H_END 0
#define FIBROBLAST_W_START 0
#define FIBROBLAST_W_END 0

// Stimulation parameters
#define STIMULATION_AMP -142.4 //3.0752 -100.0;//-100*Cm; // pA
#define STIMULATION_TOTAL_TIME 2.0 // 50.0 milliseconds
#define STIMULATION_BEGIN 2.0 // milliseconds

#define PROTOCOL 1

#if PROTOCOL == 0
#define STIMULATION_H_START 0
#define STIMULATION_H_END 0
#define STIMULATION_W_START 0
#define STIMULATION_W_END 0
#else
#if PROTOCOL == 2
#define STIMULATION_H_START 1
#define STIMULATION_H_END (H-1)
#define STIMULATION_W_START 1
#define STIMULATION_W_END 2
#else 
#if PROTOCOL == 1
#define STIMULATION_H_START 1
#define STIMULATION_H_END 2	
#define STIMULATION_W_START 1
#define STIMULATION_W_END (W-1)
#else
#error PROTOCOL must defined as 0, 1 or 2
#endif
#endif
#endif

// Space parameters
#define dW 0.1 // mm/node
#define dH 0.1 // mm/node
#define W 142
#define H 142

// Diffusion  parameters
#define DIFFUSION_COEF 0.056
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

// Global variables
double* g_pFibroblastMat = (double*)malloc(W*H*sizeof(double));
D_tensor_type* g_pTensorMat = (D_tensor_type*)malloc(W*H*sizeof(D_tensor_type));
state_variables* g_pStateVars = (state_variables*)malloc(W*H*sizeof(state_variables));
double* g_pStimulation = (double*)malloc(W*H*sizeof(double));
double* g_pRiseTimeMat = (double*)malloc(W*H*sizeof(double));

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// Helping macros
#define START_LOOP \
	int CUR_INDEX = 0; \
	for(int i = 0; i < W; i++) \
{ \
	for(int j = 0; j < H; j++) \
{ \
	CUR_INDEX = i*H + j;

#define END_LOOP } }

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void ShowVmMat()
{	
	for(int i = 0; i < W; i += 4)	
	{	
		for (int j = 0; j < H; j += 4)
		{
			int nCurIndex = i*H + j;
			if (g_pFibroblastMat[nCurIndex] != 0.0)	printf(" ");
			else
			{
				double dVm = g_pStateVars[nCurIndex].V;
				if (dVm < -80.0) printf("88");
				else if (dVm < -70.0) printf("77");
				else if (dVm < -60.0) printf("66");
				else if (dVm < -50.0) printf("55");
				else if (dVm < -40.0) printf("44");
				else if (dVm < -30.0) printf("33");
				else if (dVm < -20.0) printf("22");
				else if (dVm < -10.0) printf("11");
				else if (dVm < 0.0) printf("00");
				else if (dVm < 5.0) printf("AA");
				else if (dVm < 10.0) printf("BB");
				else if (dVm < 15.0) printf("CC");
				else if (dVm < 20.0) printf("DD");
				else if (dVm < 25.0) printf("EE");
				else if (dVm < 30.0) printf("**");
				else printf("^^");
			}
		}
		printf("\n");
	}
	printf("\n");
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void SaveModelOutputToFile(double sim_time)
{
	double sim_time_to_save = floor(sim_time + 0.5);
	char fileName[1024] = {0};
	sprintf(fileName, "Output\\ModelOutput_at_%.2f.txt", sim_time_to_save);
	FILE* pFile = fopen(fileName, "w");
	if(pFile == NULL)
	{
		return;
	}

	for(int i = 1; i < H-1; i++)
	{	
		for(int j = 1; j < W-1; j++)
		{
			fprintf(pFile, "%f ", g_pStateVars[i*H + j].V);
		}
		fprintf(pFile, "\n");
	}	
	fclose(pFile);
	pFile = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void SaveRiseTimeToFile()
{
	char fileName[1024] = {0};
	sprintf(fileName, "Output\\TargetFibroblastMatResults%d.txt", PROTOCOL);
	FILE* pFile = fopen(fileName, "w");
	if(pFile == NULL)
	{
		return;
	}

	for (int i = 1; i < W-1; i++)
	{	
		for (int j = 1; j < H-1; j++)
		{
			fprintf(pFile, "%.3f ", g_pRiseTimeMat[i*H + j]);
		}
		fprintf(pFile, "\n");
	}	
	fclose(pFile);
	pFile = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void InitStateVars()
{
	ReadStateVariablesInitialConditionFromFile(); // Put initial conditions in the vector state[...]
	START_LOOP	
		AssignInitialCondition(g_pStateVars[CUR_INDEX]);
	END_LOOP	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void InitFibroblastMat()
{
	START_LOOP
		g_pFibroblastMat[CUR_INDEX] = 0.0;
	END_LOOP
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CreateFibroblastBorders()
{
	for(int iBorder = 0; iBorder < W; ++iBorder)
	{
		g_pFibroblastMat[iBorder] = 1.0;
		g_pFibroblastMat[(H-1)*W + iBorder] = 1.0;
	}
	for(int iBorder = 0; iBorder < H; ++iBorder)
	{
		g_pFibroblastMat[iBorder*W] = 1.0;
		g_pFibroblastMat[(iBorder+1)*W - 1] = 1.0;
	}		
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CreateFibroblastPatch()
{
	for (int i = FIBROBLAST_H_START; i <= FIBROBLAST_H_END; ++i)
	{	
		for (int j = FIBROBLAST_W_START; j <= FIBROBLAST_W_END; ++j)
		{
			g_pFibroblastMat[i*H + j] = 1.0;
		}
	}	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void InitTensorMat()
{
	START_LOOP			
		if(g_pFibroblastMat[CUR_INDEX] == 0.0)
		{
			g_pTensorMat[CUR_INDEX].Dii = DIFFUSION_COEF;
			g_pTensorMat[CUR_INDEX].Djj = DIFFUSION_COEF;
			g_pTensorMat[CUR_INDEX].Dij = 0.0;
		}
		else if(g_pFibroblastMat[CUR_INDEX] == 1.0)
		{
			g_pTensorMat[CUR_INDEX].Dii = 0.0;
			g_pTensorMat[CUR_INDEX].Djj = 0.0;
			g_pTensorMat[CUR_INDEX].Dij = 0.0;
		}
		else
		{
			throw;
		}
		END_LOOP
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void InitProtcolMat()
{
	START_LOOP
		if(i >= STIMULATION_H_START && i <= STIMULATION_H_END &&
			j >= STIMULATION_W_START && j <= STIMULATION_W_END)
		{	
			g_pStimulation[CUR_INDEX] = STIMULATION_AMP;
		}
		else
		{
			g_pStimulation[CUR_INDEX] = 0.0;
		}
		END_LOOP
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void InitRiseTimeMat()
{
	START_LOOP	
		if(g_pFibroblastMat[CUR_INDEX] == 0.0)
		{
			g_pRiseTimeMat[CUR_INDEX] = INVALID_RISE_TIME;
		}
		else
		{
			g_pRiseTimeMat[CUR_INDEX] = FIBROBLAST_RISE_TIME;
		}		
		END_LOOP
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{		
	InitFibroblastMat();
	CreateFibroblastBorders();
	CreateFibroblastPatch();
	InitStateVars();
	InitTensorMat();
	InitProtcolMat();
	InitRiseTimeMat();

	double dFirstRiseTime = 0.0;
	double dStimulation = 0.0;
	bool bStimulation = false;
	double dNextTimeToSaveOutput = 0.0;
	double dVm[W*H] = {0};
	double dVm_plus_x = 0.0;
	double dVm_minus_x = 0.0;
	double dVm_plus_y = 0.0;
	double dVm_minus_y = 0.0;
	double dIion = 0.0;
	clock_t startingTime = clock();
	for (double sim_time = 0.0; sim_time < (MAX_SIMULATION_TIME + dt); sim_time += dt) 
	{	
		bStimulation = ((sim_time >= STIMULATION_BEGIN) && (sim_time < STIMULATION_BEGIN + STIMULATION_TOTAL_TIME));

		for (int x = 0; x < W; x++)
		{
			for (int y = 0; y < H; y++)
			{
				if (g_pFibroblastMat[x*H + y] == 0.0)
				{
					dVm_plus_x = (2.0*g_pTensorMat[x*H+y].Dii*g_pTensorMat[(x+1)*H+y].Dii/(g_pTensorMat[x*H+y].Dii+g_pTensorMat[(x+1)*H+y].Dii))*(g_pStateVars[(x+1)*H+y].V-g_pStateVars[x*H+y].V)/dW;
					dVm_plus_x += ((g_pTensorMat[x*H+y].Dij*g_pTensorMat[(x+1)*H+y].Dii+g_pTensorMat[x*H+y].Dii*g_pTensorMat[(x+1)*H+y].Dij)/(g_pTensorMat[x*H+y].Dii+g_pTensorMat[(x+1)*H+y].Dii))*(g_pStateVars[x*H+y+1].V+g_pStateVars[(x+1)*H+y+1].V-g_pStateVars[x*H+y-1].V-g_pStateVars[(x+1)*H+y-1].V)/(4.0*dH);
					dVm_minus_x = (2.0*g_pTensorMat[(x-1)*H+y].Dii*g_pTensorMat[x*H+y].Dii/(g_pTensorMat[(x-1)*H+y].Dii+g_pTensorMat[x*H+y].Dii))*(g_pStateVars[x*H+y].V-g_pStateVars[(x-1)*H+y].V)/dW;
					dVm_minus_x += ((g_pTensorMat[(x-1)*H+y].Dij*g_pTensorMat[x*H+y].Dii+g_pTensorMat[(x-1)*H+y].Dii*g_pTensorMat[x*H+y].Dij)/(g_pTensorMat[(x-1)*H+y].Dii+g_pTensorMat[x*H+y].Dii))*(g_pStateVars[x*H+y+1].V+g_pStateVars[(x-1)*H+y+1].V-g_pStateVars[x*H+y-1].V-g_pStateVars[(x-1)*H+y-1].V)/(4.0*dH);
					dVm_plus_y = (2.0*g_pTensorMat[x*H+y].Djj*g_pTensorMat[x*H+y+1].Djj/(g_pTensorMat[x*H+y].Djj+g_pTensorMat[x*H+y+1].Djj))*(g_pStateVars[x*H+y+1].V-g_pStateVars[x*H+y].V)/dH;
					dVm_plus_y += ((g_pTensorMat[x*H+y].Dij*g_pTensorMat[x*H+y+1].Djj+g_pTensorMat[x*H+y].Djj*g_pTensorMat[x*H+y+1].Dij)/(g_pTensorMat[x*H+y].Djj+g_pTensorMat[x*H+y+1].Djj))*(g_pStateVars[(x+1)*H+y].V+g_pStateVars[(x+1)*H+y+1].V-g_pStateVars[(x-1)*H+y].V-g_pStateVars[(x-1)*H+y+1].V)/(4.0*dW);
					dVm_minus_y = (2.0*g_pTensorMat[x*H+y-1].Djj*g_pTensorMat[x*H+y].Djj/(g_pTensorMat[x*H+y-1].Djj+g_pTensorMat[x*H+y].Djj))*(g_pStateVars[x*H+y].V-g_pStateVars[x*H+y-1].V)/dH;
					dVm_minus_y += ((g_pTensorMat[x*H+y-1].Dij*g_pTensorMat[x*H+y].Djj+g_pTensorMat[x*H+y-1].Djj*g_pTensorMat[x*H+y].Dij)/(g_pTensorMat[x*H+y-1].Djj+g_pTensorMat[x*H+y].Djj))*(g_pStateVars[(x+1)*H+y].V+g_pStateVars[(x+1)*H+y-1].V-g_pStateVars[(x-1)*H+y].V-g_pStateVars[(x-1)*H+y-1].V)/(4.0*dW);
					dVm[x*H+y] = ((dVm_plus_x-dVm_minus_x)*dH+(dVm_plus_y-dVm_minus_y)*dW)/(dW*dH);
				}
			}
		}

		START_LOOP		
			bStimulation ? dStimulation = g_pStimulation[CUR_INDEX] : dStimulation = 0.0;
		if (g_pFibroblastMat[CUR_INDEX] == 0.0)
		{
			dIion = CalcTotalTransmembraneCurrents(g_pStateVars[CUR_INDEX], dStimulation);
			g_pStateVars[CUR_INDEX].V += dt*(dVm[CUR_INDEX] - dIion);
			
			if((g_pRiseTimeMat[CUR_INDEX] == INVALID_RISE_TIME) && (g_pStateVars[CUR_INDEX].V > 0.0))
			{				
				if(dFirstRiseTime == 0.0)
				{
					dFirstRiseTime = sim_time;
				}
				g_pRiseTimeMat[CUR_INDEX] = sim_time - dFirstRiseTime;
			}
		}
		END_LOOP

			if (sim_time >= dNextTimeToSaveOutput)
			{
				clock_t endingTime = clock();
				double runningTime = (endingTime - startingTime)/double(CLOCKS_PER_SEC);
				startingTime = clock();
				dNextTimeToSaveOutput += SAVE_OUTPUT_PERIOD;
				printf("Vm mat at time: %.2f (duration: %.3f)\n", sim_time, runningTime);
				ShowVmMat();
				SaveModelOutputToFile(sim_time);

				bool bStopSim = true;
				START_LOOP					
					if (g_pFibroblastMat[CUR_INDEX] == 0.0)
					{
						if(g_pRiseTimeMat[CUR_INDEX] == INVALID_RISE_TIME)
						{
							bStopSim = false;
							break;
						}
					}
					END_LOOP
						if(bStopSim)
						{
							break;
						}			
			}
	} // sim_time loop

	SaveRiseTimeToFile();
	return 0;
} // of main

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
