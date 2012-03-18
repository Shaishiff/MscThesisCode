
#ifndef __FK__MODEL__
#define __FK__MODEL__

// #define SAVE_FK_MODEL_TO_FILE
#include "defs.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

class CFkModel
{
public:
	CFkModel();
	~CFkModel();
	void CleanupFkModel();
	double** ExecuteFkModel(double** inFibroblastMat, const ProtocolParams& protParams);
	void SaveToFile(int iT, int& nFileNumber, double** mat, char* nameOfVar);
	static bool SaveToInputFile(double** mat, char* fileName);
	static bool SaveToOutputFile(double** mat, char* fileName);
	void InitTargetFibroblastMat();
	void InitTargetFibroblastMat(int nCenterH, int nCenterW, int nHeight, int nWidth);

private:
	void InitFkModel();
	void DeleteFkModel();
	void AddFibroblasts(int hStart, int hEnd, int wStart, int wEnd);
	void CalculateDer(int iH, int iW, double** inFibroblastMat);
	static bool SaveToFile(double** mat, char* fileName);

private:
	double dVdh;
	double dVdw;
	double dW2;
	double dH2;
	double** mat;
	double** fibroMat;
	double** new_mat;
	double** s_mat;
	double** new_s_mat;
	double** f_mat;
	double** new_f_mat;
	double** rise_time_mat;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __FK__MODEL__

