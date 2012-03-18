
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
	virtual ~CFkModel();
	void ExecuteModel(double** inFibroblastMat, double** outRiseTimeMat, const ProtocolParams& protParams);
	void SaveToFile(int iT, int& nFileNumber, double** mat, char* nameOfVar);
	static bool SaveToInputFile(double** mat, char* fileName);
	static bool SaveToOutputFile(double** mat, char* fileName);
	static bool SaveToFile(double** mat, char* fileName);

private:
	void InitFkModel();
	void DeleteFkModel();
	void CleanupFkModel();
	void CalculateDer(int iH, int iW, double** inFibroblastMat);	

private:
	double dVdh;
	double dVdw;
	double dW2;
	double dH2;
	double** mat;
	double** new_mat;
	double** s_mat;
	double** new_s_mat;
	double** f_mat;
	double** new_f_mat;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __FK__MODEL__

