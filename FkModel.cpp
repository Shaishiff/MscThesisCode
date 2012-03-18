// FkModel.cpp : Defines the entry point for the console application.
//

#include "defs.h"
#include "FkModel.h"

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
	mat = new double*[Nh+2];
	new_mat = new double*[Nh+2];
	s_mat = new double*[Nh+2];
	new_s_mat = new double*[Nh+2];
	f_mat = new double*[Nh+2];
	new_f_mat = new double*[Nh+2];
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		mat[iH] = new double[Nw+2];
		new_mat[iH] = new double[Nw+2];
		s_mat[iH] = new double[Nw+2];
		new_s_mat[iH] = new double[Nw+2];
		f_mat[iH] = new double[Nw+2];
		new_f_mat[iH] = new double[Nw+2];
		for(int iW = 0; iW < Nw+2; ++iW)
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

void CFkModel::DeleteFkModel()
{
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		delete [] mat[iH];
		delete [] new_mat[iH];
		delete [] s_mat[iH];
		delete [] new_s_mat[iH];
		delete [] f_mat[iH];
		delete [] new_f_mat[iH];
	}

	delete [] mat;
	delete [] new_mat;
	delete [] s_mat;
	delete [] new_s_mat;
	delete [] f_mat;
	delete [] new_f_mat;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CFkModel::CleanupFkModel()
{
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		for(int iW = 0; iW < Nw+2; ++iW)
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

bool CFkModel::SaveToInputFile(double** mat, char* fileName)
{
	char fullFilePath[FILE_NAME_BUFFER_SIZE] = {0};
	sprintf(fullFilePath, "C:\\Users\\Shai\\Documents\\MATLAB\\MSc\\FkModel\\Input\\%s", fileName);
	return SaveToFile(mat, fullFilePath);
}

bool CFkModel::SaveToOutputFile(double** mat, char* fileName)
{
	char fullFilePath[FILE_NAME_BUFFER_SIZE] = {0};
	sprintf(fullFilePath, "C:\\Users\\Shai\\Documents\\MATLAB\\MSc\\FkModel\\Output\\%s", fileName);
	return SaveToFile(mat, fullFilePath);
}

bool CFkModel::SaveToFile(double** mat, char* fileName)
{
	FILE* pFile = fopen(fileName, "w");
	if(pFile == NULL)
	{
		return false;
	}

	for (int iH = 1; iH < Nh+1; ++iH)
	{
		for (int iW = 1; iW < Nw+1; ++iW)
		{
			fprintf(pFile, "%f ", mat[iH][iW]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

const double SaveToFileInterval = 5.0;

void CFkModel::SaveToFile(int iT, int& nFileNumber, double** mat, char* nameOfVar)
{
#ifndef SAVE_FK_MODEL_TO_FILE
	return;
#endif

	if((nFileNumber*SaveToFileInterval) <= (iT*dt))
	{		
		char fileName[FILE_NAME_BUFFER_SIZE] = {0};
		int n = sprintf(fileName, "C:\\Users\\Shai\\Documents\\MATLAB\\MSc\\FkModel\\Output\\SimulationOutput_%s_%d.txt", nameOfVar, nFileNumber);
		if(SaveToFile(mat, fileName))
		{
			++nFileNumber;
			// printf("Saved data to file %d at time: %.4f\n", nFileNumber-1, (iT*dt));
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

void CFkModel::ExecuteModel(double** inFibroblastMat, double** outRiseTimeMat, const ProtocolParams& protParams)
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
	int nFibroFileNumber = 0;
	int nuFileNumber = 0;
	int nsFileNumber = 0;
	int nfFileNumber = 0;

	SaveToFile(0, nFibroFileNumber, inFibroblastMat, "Fibroblasts"); 

	// Start the temporal loop.
	for (int iT = 0; iT < Nt; ++iT)
	{   		
		// Save the current data to the file.	
		//SaveToFile(iT, nuFileNumber, mat, "u"); 
		//SaveToFile(iT, nsFileNumber, s_mat, "s"); 
		//SaveToFile(iT, nfFileNumber, f_mat, "f"); 
		
		// Calc the time in ms and check if we need to use the S1 protocol.
		double curTime = iT*dt;
		bool bS1 = ((curTime >= protParams.m_BeginTime) && (curTime <= (protParams.m_BeginTime+protParams.m_TotalTime)));
		
		// Start the spatial loop.
		for (int iH = 1; iH < Nh+1; ++iH)
		{
			for (int iW = 1; iW < Nw+1; ++iW)
			{            
				if(inFibroblastMat[iH][iW] == 1)
				{
					outRiseTimeMat[iH][iW] = TotalSimulationTime*2;
				}
				else
				{
					CalculateDer(iH, iW, inFibroblastMat);
					
					double Jstim = 0.0;
					if(bS1)
					{
						double curH = iH*dH;
						double curW = iW*dW;
						if(curH >= protParams.m_hStart && curH <= protParams.m_hEnd && curW >= protParams.m_wStart && curW <= protParams.m_wEnd)
						{
							Jstim = S1Amp;
						}
					}

					double v = mat[iH][iW];
					if((v > 0.95) && (outRiseTimeMat[iH][iW] == 0.0))
					{
						if(Jstim == 0.0)
						{
							outRiseTimeMat[iH][iW] = iT*dt;
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
	}

	//SaveToFile(Nt, nuFileNumber, mat, "u"); 
	//SaveToFile(Nt, nsFileNumber, s_mat, "s"); 
	//SaveToFile(Nt, nfFileNumber, f_mat, "f");
	
	//int n = 0;
	//SaveToFile(Nt, n, outRiseTimeMat, "r");

	// return mat;
	//return rise_time_mat;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
