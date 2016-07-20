
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
#define RISE_TIME_VM_THRESHOLD (-20.0)

// Fibroblasts patches.
#define FIBROBLAST_PATCHES 2

// (43_43) - (48_48)_ + _(93_93) - (98_98): Done
// (41_41) - (50_50)_ + _(91_91) - (100_100): Done
// (36_36) - (55_55)_ + _(86_86) - (105_105): Done
// (31_31) - (60_60)_ + _(81_81) - (110_110): Done
// (26_26) - (65_65)_ + _(76_76) - (115_115): Done
// (22_22) - (69_69)_ + _(72_72) - (119_119): Done

#define FIBROBLAST_H_START_PATCH_0 36
#define FIBROBLAST_W_START_PATCH_0 FIBROBLAST_H_START_PATCH_0

#define FIBROBLAST_H_END_PATCH_0 55
#define FIBROBLAST_W_END_PATCH_0 FIBROBLAST_H_END_PATCH_0

#define FIBROBLAST_H_START_PATCH_1 86
#define FIBROBLAST_W_START_PATCH_1 FIBROBLAST_H_START_PATCH_1

#define FIBROBLAST_H_END_PATCH_1 105
#define FIBROBLAST_W_END_PATCH_1 FIBROBLAST_H_END_PATCH_1 

#define FIBROBLAST_H_START_PATCH_2 0
#define FIBROBLAST_H_END_PATCH_2 0
#define FIBROBLAST_W_START_PATCH_2 0
#define FIBROBLAST_W_END_PATCH_2 0

#define FIBROBLAST_H_START_PATCH_3 0
#define FIBROBLAST_H_END_PATCH_3 0
#define FIBROBLAST_W_START_PATCH_3 0
#define FIBROBLAST_W_END_PATCH_3 0

// Disc fibroblasts.
#define FIBROBLAST_H_CENTER 60
#define FIBROBLAST_W_CENTER 60
#define FIBROBLAST_RADIUS   0

#define PROT_0_STIMULATION_H_START 0
#define PROT_0_STIMULATION_H_END 0
#define PROT_0_STIMULATION_W_START 0
#define PROT_0_STIMULATION_W_END 0

#define PROT_1_STIMULATION_H_START 1
#define PROT_1_STIMULATION_H_END 2	
#define PROT_1_STIMULATION_W_START 1
#define PROT_1_STIMULATION_W_END (W-1)

#define PROT_2_STIMULATION_H_START 1
#define PROT_2_STIMULATION_H_END (H-1)
#define PROT_2_STIMULATION_W_START 1
#define PROT_2_STIMULATION_W_END 2

// Space parameters
#define dW 0.01 // mm/node
#define dH 0.01 // mm/node
#define dZ 0.01 // mm/node
#define W 142
#define H 142

// Stimulation parameters
#define STIMULATION_AMP (-100.0) //3.0752 -100.0;//-100*Cm; // pA
#define STIMULATION_TOTAL_TIME 2.0 // 50.0 milliseconds
#define STIMULATION_BEGIN 5.0 // milliseconds

// Diffusion  parameters
#define DIFFUSION_COEF 0.001//0.00056
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

/*
double* g_pGardientX = (double*)malloc(W*H*sizeof(double));
double* g_pGardientY = (double*)malloc(W*H*sizeof(double));
double* g_pExtraCellularV = (double*)malloc(W*H*sizeof(double));
double* g_pExtraCellularRiseTimeMat = (double*)malloc(W*H*sizeof(double));
*/
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
	printf("\n");
	for(int i = 1; i < W-1; i += 4)	
	{	
		printf("\t");
		for (int j = 1; j < H-1; j += 4)
		{
			int nCurIndex = i*H + j;
			if (g_pFibroblastMat[nCurIndex] != 0.0)	
			{
				printf("  ");
			}
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

void SaveModelOutputToFile(double sim_time, int nProtocol)
{
	double sim_time_to_save = floor(sim_time + 0.5);
	char fileName[1024] = {0};
	sprintf(fileName, "Output\\ModelLogs\\Prot%d\\ModelOutput_at_%.2f.txt", nProtocol, sim_time_to_save);
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

void SaveModelOutputToFile2(double sim_time, int nProtocol)
{
	char fileName[1024] = {0};
	sprintf(fileName, "Output\\ModelLogs\\Prot%d\\ModelMeasurments.txt", nProtocol);
	FILE* pFile = sim_time == 0.0 ? fopen(fileName, "w") : fopen(fileName, "a");
	if(pFile == NULL)
	{
		return;
	}
	if(sim_time != 0.0)
	{
		fprintf(pFile, "\n");
	}
	if(nProtocol == 1)
	{
		for (int j = 1; j < H-1; j++)
		{	
			int i = (W-2-20);
			fprintf(pFile, "%.6f ", g_pStateVars[i*H + j].V);
		}				
	}
	else if(nProtocol == 2)
	{
		for (int i = 1; i < W-1; i++)
		{	
			int j = (H-2-20);
			fprintf(pFile, "%.6f ", g_pStateVars[i*H + j].V);
		}				
	}	
	fclose(pFile);
	pFile = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void SaveModelOutputToFile3(double sim_time, int nProtocol)
{
	char fileName[1024] = {0};
	sprintf(fileName, "Output\\ModelLogs\\Prot%d\\ModelOutput_crunch_time_at_%.3f.txt", nProtocol, sim_time);
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

void SaveFibroblastsToFile()
{
	char fileName[1024] = {0};
	sprintf(fileName, "Output\\TargetFibroblastMat.txt");
	FILE* pFile = fopen(fileName, "w");
	if(pFile == NULL)
	{
		return;
	}

	for (int i = 1; i < W-1; i++)
	{	
		for (int j = 1; j < H-1; j++)
		{
			fprintf(pFile, "%.3f ", g_pFibroblastMat[i*H + j]);
		}
		fprintf(pFile, "\n");
	}	
	fclose(pFile);
	pFile = NULL;
	
	sprintf(fileName, "Output\\TargetFibroblastMatReadable.txt");
	pFile = fopen(fileName, "w");
	if(pFile == NULL)
	{
		return;
	}

	for (int i = 0; i < W; i++)
	{	
		for (int j = 0; j < H; j++)
		{
			fprintf(pFile, "%d", (int)g_pFibroblastMat[i*H + j]);
		}
		fprintf(pFile, "\n");
	}	
	fclose(pFile);
	pFile = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void SaveRiseTimeToFile(int nProtocol)
{
	char fileName[1024] = {0};
	sprintf(fileName, "Output\\TargetFibroblastMatResults%d.txt", nProtocol);
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
	ReadStateVariablesInitialConditionFromFile();
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
#if FIBROBLAST_PATCHES
	int nHStart[FIBROBLAST_PATCHES] = {0};
	int nHEnd[FIBROBLAST_PATCHES] = {0};
	int nWStart[FIBROBLAST_PATCHES] = {0};
	int nWEnd[FIBROBLAST_PATCHES] = {0};

	nHStart[0] = FIBROBLAST_H_START_PATCH_0;
	nHEnd[0] = FIBROBLAST_H_END_PATCH_0;
	nWStart[0] = FIBROBLAST_W_START_PATCH_0;
	nWEnd[0] = FIBROBLAST_W_END_PATCH_0;

	nHStart[1] = FIBROBLAST_H_START_PATCH_1;
	nHEnd[1] = FIBROBLAST_H_END_PATCH_1;
	nWStart[1] = FIBROBLAST_W_START_PATCH_1;
	nWEnd[1] = FIBROBLAST_W_END_PATCH_1;

	nHStart[2] = FIBROBLAST_H_START_PATCH_2;
	nHEnd[2] = FIBROBLAST_H_END_PATCH_2;
	nWStart[2] = FIBROBLAST_W_START_PATCH_2;
	nWEnd[2] = FIBROBLAST_W_END_PATCH_2;

	nHStart[3] = FIBROBLAST_H_START_PATCH_3;
	nHEnd[3] = FIBROBLAST_H_END_PATCH_3;
	nWStart[3] = FIBROBLAST_W_START_PATCH_3;
	nWEnd[3] = FIBROBLAST_W_END_PATCH_3;

	for(int iPatch = 0; iPatch < FIBROBLAST_PATCHES; ++iPatch)
	{
		for (int i = nHStart[iPatch]; i <= nHEnd[iPatch]; ++i)
		{	
			for (int j = nWStart[iPatch]; j <= nWEnd[iPatch]; ++j)
			{
				g_pFibroblastMat[i*H + j] = 1.0;
			}
		}
	}
#endif

	if(FIBROBLAST_RADIUS > 0)
	{
		START_LOOP		
		double dRadius = sqrt(pow(double(FIBROBLAST_H_CENTER - j), int(2)) + pow(double(FIBROBLAST_W_CENTER - i), int(2)));
		if(dRadius <= FIBROBLAST_RADIUS)
		{
			g_pFibroblastMat[CUR_INDEX] = 1.0;
		}
		END_LOOP
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

void InitProtcolMat(int nProtocol)
{
	int nHStart = (nProtocol == 0 ? PROT_0_STIMULATION_H_START : (nProtocol == 1 ? PROT_1_STIMULATION_H_START : PROT_2_STIMULATION_H_START));
	int nWStart = (nProtocol == 0 ? PROT_0_STIMULATION_W_START : (nProtocol == 1 ? PROT_1_STIMULATION_W_START : PROT_2_STIMULATION_W_START));
	int nHEnd = (nProtocol == 0 ? PROT_0_STIMULATION_H_END : (nProtocol == 1 ? PROT_1_STIMULATION_H_END : PROT_2_STIMULATION_H_END));
	int nWEnd = (nProtocol == 0 ? PROT_0_STIMULATION_W_END : (nProtocol == 1 ? PROT_1_STIMULATION_W_END : PROT_2_STIMULATION_W_END));
	START_LOOP
		if(i >= nHStart && i <= nHEnd &&
			j >= nWStart && j <= nWEnd)
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
	clock_t total_program_starting_time = clock();
	for(int iProtocol = 1; iProtocol <= 2; ++iProtocol)
	{
		printf("\n*******************\nStarted protocol %d\n*******************\n\n", iProtocol);

		InitFibroblastMat();
		CreateFibroblastBorders();
		CreateFibroblastPatch();
		SaveFibroblastsToFile();
		InitStateVars();
		InitTensorMat();
		InitProtcolMat(iProtocol);
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
			SaveModelOutputToFile2(sim_time, iProtocol);

			if (sim_time >= 28.0)
			{
				SaveModelOutputToFile3(sim_time, iProtocol);
			}

			if (sim_time >= dNextTimeToSaveOutput)
			{
				clock_t endingTime = clock();
				double runningTime = (endingTime - startingTime)/double(CLOCKS_PER_SEC);
				startingTime = clock();
				dNextTimeToSaveOutput += SAVE_OUTPUT_PERIOD;
				printf("Vm mat at time: %.2f (duration: %.3f)\n", sim_time, runningTime);
				ShowVmMat();
				SaveModelOutputToFile(sim_time, iProtocol);

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
			
					if((g_pRiseTimeMat[CUR_INDEX] == INVALID_RISE_TIME) && (g_pStateVars[CUR_INDEX].V > RISE_TIME_VM_THRESHOLD))
					{				
						if(dFirstRiseTime == 0.0)
						{
							dFirstRiseTime = sim_time;
						}
						g_pRiseTimeMat[CUR_INDEX] = sim_time - dFirstRiseTime;
					}
				}
			END_LOOP

			/*
			START_LOOP
				if(j != 0 && j != (H-1))
				{
					g_pGardientX[CUR_INDEX] = (g_pStateVars[i*H + j + 1].V - g_pStateVars[i*H + j - 1].V)/2;
				}
				else if(j == 0)
				{
					g_pGardientX[CUR_INDEX] = (g_pStateVars[i*H + j + 1].V - g_pStateVars[i*H + j].V);
				}
				else // j == H-1
				{
					g_pGardientX[CUR_INDEX] = (g_pStateVars[i*H + j].V - g_pStateVars[i*H + j - 1].V);
				}

				if(i != 0 && i != (W-1))
				{
					g_pGardientY[CUR_INDEX] = (g_pStateVars[(i+1)*H + j].V - g_pStateVars[(i-1)*H + j].V)/2;
				}
				else if(i == 0)
				{
					g_pGardientY[CUR_INDEX] = (g_pStateVars[(i+1)*H + j].V - g_pStateVars[i*H + j].V);
				}
				else // i == W-1
				{
					g_pGardientY[CUR_INDEX] = (g_pStateVars[i*H + j].V - g_pStateVars[(i-1)*H + j].V);
				}
			END_LOOP

			START_LOOP
				for(int ii = 0; ii < W; ii++)
				{
					for(int jj = 0; jj < H; jj++)
					{
						double dDistanceW = double(dW*(i-ii));
						double dDistanceH = double(dH*(j-jj));
						double dDistanceR3 = pow((pow(dDistanceW,2) + pow(dDistanceH,2) + dZ*dZ), 3/2);
						double dUnipolar = (g_pGardientY[ii*H + jj]*dDistanceW + g_pGardientX[ii*H + jj]*dDistanceH)/dDistanceR3;
					}
				} 
			END_LOOP
			*/

		} // sim_time loop

		SaveRiseTimeToFile(iProtocol);
	}

	clock_t total_program_ending_time = clock();
	double total_program_running_time = (total_program_ending_time - total_program_starting_time)/double(CLOCKS_PER_SEC);
	printf("Total program running time: %.3f\n", total_program_running_time);
	getchar();
	return 0;
} // of main

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
