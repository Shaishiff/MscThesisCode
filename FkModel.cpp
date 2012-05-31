// FkModel.cpp : Defines the entry point for the console application.
//

#include "defs.h"
#include "FkModel.h"
#include "Mat.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

CFkModel::CFkModel()
{
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
	InitFkModel();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

CFkModel::~CFkModel()
{
	DeleteFkModel();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CFkModel::InitFkModel()
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

void CFkModel::DeleteFkModel()
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

void CFkModel::CleanupFkModel()
{
	for(int iH = 0; iH < Nh_with_border; ++iH)
	{
		for(int iW = 0; iW < Nw_with_border; ++iW)
		{
			mat[iH][iW] = 0.0;			
			new_mat[iH][iW] = 0.0;
			s_mat[iH][iW] = 0.0;
			new_s_mat[iH][iW] = 0.0;
			f_mat[iH][iW] = 0.0;
			new_f_mat[iH][iW] = 0.0;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CFkModel::CalculateDer(int iH, int iW, double** inFibroblastMat)
{
	if(inFibroblastMat[iH-1][iW] == 1)
	{
		dVdh = (mat[iH+1][iW] - mat[iH][iW])/dH2;
	}
	else if (inFibroblastMat[iH + 1][iW] == 1)
	{
		dVdh = (mat[iH-1][iW] - mat[iH][iW])/dH2;
	}
	else
	{
		dVdh = (mat[iH-1][iW] + mat[iH+1][iW] - 2*mat[iH][iW])/dH2;
	}

	if(inFibroblastMat[iH][iW-1] == 1)
	{
		dVdw = (mat[iH][iW+1] - mat[iH][iW])/dW2;
	}
	else if (inFibroblastMat[iH][iW+1] == 1)
	{
		dVdw  = (mat[iH][iW-1] - mat[iH][iW])/dW2;
	}
	else
	{
		dVdw  = (mat[iH][iW-1] + mat[iH][iW+1] - 2*mat[iH][iW])/dW2;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#define INVALID_RISE_TIME	1000.0

void CFkModel::ExecuteModel(double** inFibroblastMat, double** outRiseTimeMat, const ProtocolParams& protParams, char* outputFolder)
{	
	CleanupFkModel();
	double t_ung = 0.0;
	double J_ung = 0.0;
	double s_inf = 0.0;
	double tau_s = 0.0;
	double J_slow = 0.0;
	double f_inf = 0.0;
	double tau_f = 0.0;
	double J_fast  = 0.0;
	double Jion = 0.0;
	double dFirstRiseTime = 0.0;
	
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
				//sprintf(modelOutput, "ModelOutput_at_%.2f.txt", curTime);
				//SaveMatToFile(mat, modelOutput, outputFolder);
				
				sprintf(modelOutput, "%s\\ModelOutput_at_%.2f.txt", outputFolder, curTime);				
				SaveMatToFileWithFullName(mat, modelOutput);
				printf("Saved file %s\n", modelOutput);
			}
		}
		
		// Start the spatial loop.
		for (int iH = 1; iH < Nh+1; ++iH)
		{
			for (int iW = 1; iW < Nw+1; ++iW)
			{            
				if(inFibroblastMat[iH][iW] != 1)
				{					
					CalculateDer(iH, iW, inFibroblastMat);
					
					double Jstim = 0.0;
					if(bS1)
					{
						double curH = iH*dH;
						double curW = iW*dW;
						if(curH >= protParams.m_hStart &&
						   curH <= protParams.m_hEnd &&
						   curW >= protParams.m_wStart &&
						   curW <= protParams.m_wEnd)
						{
							Jstim = S1Amp;
						}
					}

					double v = mat[iH][iW];
					if((v > 0.3) && (outRiseTimeMat[iH][iW] == INVALID_RISE_TIME))
					{
						//if(Jstim == 0.0)
						{							
							if(dFirstRiseTime == 0.0)
							{
								dFirstRiseTime = iT*dt;
							}
							outRiseTimeMat[iH][iW] = iT*dt - dFirstRiseTime;
						}
					}

					double s = s_mat[iH][iW];
					double f = f_mat[iH][iW];

					 // Ungated  Outward  Current
					if (v < V_u)
					{
						t_ung = t_ug1; // tau_ungated
						J_ung = (1 / t_ung) * v; // J_ungated
					}
					else // v > V_u
					{
						t_ung =t_ug2; // tau_ungated
						J_ung = (1 / t_ung) * 1; // J_ungated
					}

					// Slow  Gate  Current
					if (v < V_s)
					{
						s_inf = 1;
						tau_s = ts_opn;
					}
					else // v > V_s
					{
						s_inf = 0;
						tau_s = ts_cls;
					}

					// Slow  Inward  Current
					J_slow = -(s / t_slow) * (1+tanh(kwm * (v - V_slow) ) ) / 2;

					// Fast  Gate  Current
					if (v < V_f)
					{
						f_inf = 1;
						if (v < Vf_open)
							tau_f = tf_opn1;
						else // v > Vf_open
							tau_f = tf_opn2;
					}
					else // v < V_f
					{
						f_inf = 0;
						tau_f = tf_cls;
					}

					// Fast  Inward  Current
					if (v < V_f)
					{
						J_fast = 0;
					}
					else
					{
						J_fast = -(f / t_fast) * (v - V_f) * (1 - v);
					}

					Jion = J_fast + J_slow + J_ung;

					new_mat[iH][iW] = mat[iH][iW] + ((dt*sigma)/(Cm*Am))*(dVdh + dVdw) - (dt/Cm)*Jion + (dt/(Am*Cm))*Jstim;
					new_s_mat[iH][iW] = s_mat[iH][iW] + dt*(s_inf - s)/tau_s;                        
					new_f_mat[iH][iW] = f_mat[iH][iW] + dt*(f_inf - f)/tau_f;                
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

		bool bStopSim = true;
		for (int iH = 1; iH < Nh+1; ++iH)
		{
			for (int iW = 1; iW < Nw+1; ++iW)
			{            									
				double curH = iH*dH;
				double curW = iW*dW;
				if(curH >= protParams.m_hMeasureStart &&
				   curH <= protParams.m_hMeasureEnd &&
				   curW >= protParams.m_wMeasureStart &&
				   curW <= protParams.m_wMeasureEnd)
				{
					if(outRiseTimeMat[iH][iW] == INVALID_RISE_TIME)
					{
						bStopSim = false;
						break;
					}
				}
			}
		}

		if(bStopSim)
		{
			break;
		}		
	}	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
