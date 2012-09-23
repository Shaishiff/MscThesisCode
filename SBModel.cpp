// FkModel.cpp : Defines the entry point for the console application.
//

#include "defs.h"
#include "SBModel.h"
#include "Mat.h"
#include "SBModelDefs.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

CSBModel::CSBModel()
{
	m_dDiffusion = Diffusion;
	m_dj = j_var;
	dVdh = 0.0;
	dVdw = 0.0;
	dW2 = dW*dW;
	dH2 = dH*dH;	
	mat = NULL;
	new_mat = NULL;
	s_mat = NULL;
	new_s_mat = NULL;
	f_mat = NULL;
	new_f_mat = NULL;
	InitModel();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

CSBModel::~CSBModel()
{
	DeleteModel();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CSBModel::InitModel()
{
	mat = CreateMat();
	new_mat = CreateMat();
	s_mat = CreateMat();
	new_s_mat = CreateMat();
	f_mat = CreateMat();
	new_f_mat = CreateMat();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CSBModel::DeleteModel()
{
	DestroyMat(mat);
	DestroyMat(new_mat);
	DestroyMat(s_mat);
	DestroyMat(new_s_mat);
	DestroyMat(f_mat);
	DestroyMat(new_f_mat);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CSBModel::CleanupModel()
{
	for(int iH = 0; iH < Nh_with_border; ++iH)
	{
		for(int iW = 0; iW < Nw_with_border; ++iW)
		{
			mat[iH][iW] = V_rest;			
			new_mat[iH][iW] = V_rest;
			s_mat[iH][iW] = 1.0;
			new_s_mat[iH][iW] = 1.0;
			f_mat[iH][iW] = 0.0;
			new_f_mat[iH][iW] = 0.0;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CSBModel::CalculateDer(int iH, int iW, double** inFibroblastMat)
{
	if(inFibroblastMat[iH-1][iW] == 1.0)
	{
		dVdh = (mat[iH+1][iW] - mat[iH][iW])/dH2;
	}
	else if (inFibroblastMat[iH + 1][iW] == 1.0)
	{
		dVdh = (mat[iH-1][iW] - mat[iH][iW])/dH2;
	}
	else
	{
		dVdh = (mat[iH-1][iW] + mat[iH+1][iW] - 2.0*mat[iH][iW])/dH2;
	}

	if(inFibroblastMat[iH][iW-1] == 1.0)
	{
		dVdw = (mat[iH][iW+1] - mat[iH][iW])/dW2;
	}
	else if (inFibroblastMat[iH][iW+1] == 1.0)
	{
		dVdw  = (mat[iH][iW-1] - mat[iH][iW])/dW2;
	}
	else
	{
		dVdw  = (mat[iH][iW-1] + mat[iH][iW+1] - 2.0*mat[iH][iW])/dW2;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#define EPSILON	double(0.000001)

/***************************
* Schraudolph's algorithm: *
***************************/
static union
{
	double d;
	struct
	{
		int j, i; // #ifdef LITTLE_ENDIAN
	} n;
} _eco;

#define EXP_A (1048576/M_LN2)	/* 2^20/LN(2) use 1512775 for integer version */
#define EXP_C 60801	 /* see text for choice of c values */
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)

#if 1
#define shaisExp(in_var) EXP(in_var)
/*
double shaisExp(double in_var)
{
#define N 6.0
	double dRet = 1.0 + in_var;
	double dMul = 1.0;
	double dPow = in_var;
	for(double iN = 2.0; iN <= N; iN += 1.0)
	{
		dMul = dMul*iN;
		dPow = dPow*in_var;
		dRet += dPow/dMul;
	}
	return dRet;
	//return (1.0 + in_var + (in_var*in_var)/2.0 + (in_var*in_var*in_var)/6.0);
	//return exp(in_var);
}
*/
#else 
#define shaisExp(in_var) exp(in_var)
#endif

bool heavySide(double in_var)
{
	if(in_var > 0.0)
        return true;
    else return false;
/*
    if(in_var > 0.0)
        return 1.0;
    else if (in_var < 0.0) 
        return 0.0;
    else return 0.5;
	*/
}

double calc_alpha_h(double V)
{
	if(heavySide(-V-40.0)) return 0.0;
    return 0.135*shaisExp(-(V+80.0)/6.8);
}

double calc_alpha_m(double V)
{
	double dividend = 0.32*(V+47.13);
	double divisor = 1.0-shaisExp(-0.1*(V+47.13));
	if(fabs(divisor) < EPSILON) 
	{
		//printf("EPSILON for %.6f where dividend= %.6f and divisor=%.6f\n", V, dividend, divisor);
		divisor = EPSILON;
	}
    return dividend/divisor;
}

double calc_beta_h(double V)
{
	double dRet = 0.0;
	if(heavySide(-V-40.0))
	{
		dRet += 3.56*shaisExp(0.079*V);
		dRet += 3.1*100000.0*shaisExp(0.35*V);
	}
	if(heavySide(V+40.0))
	{
		dRet += (1.0/(0.13*(1.0+shaisExp(-(V+10.66)/11.1))));
	}
	return dRet;
    //return (3.56*shaisExp(0.079*V) + 3.1*100000.0*shaisExp(0.35*V))*heavySide(-V-40.0) + heavySide(V+40.0)*(1.0/(0.13*(1.0+shaisExp(-(V+10.66)/11.1))));
}

double calc_beta_m(double V)
{
    return 0.08*shaisExp(-V/11.0);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#define INVALID_RISE_TIME	1000.0

void CSBModel::ExecuteModel(double** inFibroblastMat, double** outRiseTimeMat, const ProtocolParams& protParams, char* outputFolder)
{			
	CleanupModel();
	
	// Init the rise time mat.
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
	
	double dFirstRiseTime = 0.0;	
	// Start the temporal loop.
	for (int iT = 0; iT < Nt; ++iT)
	{   		
		// Calc the time in ms and check if we need to use the S1 protocol.
		double curTime = iT*dt;
		bool bS1 = ((curTime >= protParams.m_BeginTime) && (curTime <= (protParams.m_BeginTime+protParams.m_TotalTime)));
		
		if(outputFolder != NULL)
		{
			if(iT%200 == 0)
			{
				char modelOutput[FILE_NAME_BUFFER_SIZE] = {0};				
				sprintf(modelOutput, "%s/ModelOutput_at_%.2f.txt", outputFolder, curTime);				
				SaveMatToFileWithFullName(mat, modelOutput);
				printf("Saved file %s\n", modelOutput);
			}
		}
		
		// Start the spatial loop.
		for (int iW = 1; iW < Nw+1; ++iW)		
		{
			for (int iH = 1; iH < Nh+1; ++iH)
			{            
				if(inFibroblastMat[iH][iW] != 1)
				{										
					// Check if we need to add a stimulation current
					// according to the protocol.
					double Jstim = 0.0;
					if(bS1)
					{
						if(iH >= protParams.m_hStart &&
						   iH <= protParams.m_hEnd &&
						   iW >= protParams.m_wStart &&
						   iW <= protParams.m_wEnd)
						{
							Jstim = S1Amp*100.0;
						}
					}

					// Search for the rise time of each cell.
					double v = mat[iH][iW];
					if((outRiseTimeMat[iH][iW] == INVALID_RISE_TIME) && (v > -50.0))
					{
						if(dFirstRiseTime == 0.0)
						{
							dFirstRiseTime = iT*dt;
						}
						outRiseTimeMat[iH][iW] = iT*dt - dFirstRiseTime;
					}

					// Start calculating the ion currents.
					double s = s_mat[iH][iW];
					double f = f_mat[iH][iW];

					// Ungated  Outward  Current
					double I_na = g_na*(V_na - v);
                    double alpha_h = calc_alpha_h(v);
                    double alpha_m = calc_alpha_m(v);
                    double beta_h = calc_beta_h(v);
                    double beta_m = calc_beta_m(v);
					double tau_h = 1.0/(alpha_h + beta_h);
                    double tau_m = 1.0/(alpha_m + beta_m);
					
					CalculateDer(iH, iW, inFibroblastMat);
					                    
					// h = s
					// m = f					
					// Calculate the ODEs.
					new_mat[iH][iW] = mat[iH][iW] + dt*(m_dDiffusion*(dVdh + dVdw) + I_na*m_dj*s*(f*f*f) + Jstim);
					new_s_mat[iH][iW] = dt*((heavySide(V_h-v)-s)/tau_h) + s;
					new_f_mat[iH][iW] = dt*((heavySide(v-V_m)-f)/tau_m) + f;
									
					if(isnan(dVdh) || isnan(alpha_h) || isnan(alpha_m) || isnan(beta_h) || isnan(beta_m) ||
					   isnan(new_mat[iH][iW]) || isnan(new_s_mat[iH][iW]) || isnan(new_f_mat[iH][iW]))
					{
						printf("********* ERROR **********\n");
						printf("NAN at iH=%d, iW=%d\n", iH, iW);
						for(int xx = iH-1; xx <= iH+1; xx++)
						{
							for(int yy = iW-1; yy <= iW+1; yy++)
							{
								printf("V(%d,%d)=%.3f  ", xx, yy, mat[xx][yy]);
							}
							printf("\n");
						}
						printf("%d with %.3f: alpha_h: %.4f, alpha_m: %.4f, beta_h: %.4f, beta_m: %.4f\n", iT, v, alpha_h, alpha_m, beta_h, beta_m);
						printf("%d with %.3f: tau_h: %.4f, tau_m: %.4f, dVdh: %.4f, dVdw: %.4f\n", iT, v, tau_h, tau_m, dVdh, dVdw);
						printf("**************************\n");
						throw;
					}
				}
			}
		}
		
		// Move the current mats to be the old ones.
		double** temp = NULL;

		temp = mat;
		mat = new_mat;
		new_mat = temp;

		temp = s_mat;
		s_mat = new_s_mat;
		new_s_mat = temp;

		temp = f_mat;
		f_mat = new_f_mat;
		new_f_mat = temp;

		// Check if we can stop the simulation.
		bool bStopSim = true;
		for (double curH = protParams.m_hMeasureStart; curH <= protParams.m_hMeasureEnd; curH += dH)
		{
			for (double curW = protParams.m_wMeasureStart; curW <= protParams.m_wMeasureEnd; curW += dW)
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
			return;
		}
	}
	
	// If should get here normally. This means that
	// the simulation time wasn't long enough.
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
