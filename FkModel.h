
#ifndef __FK__MODEL__
#define __FK__MODEL__

#include "defs.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

class CFkModel
{
public:
	CFkModel();
	virtual ~CFkModel();
	void ExecuteModel(double** inFibroblastMat, double** outRiseTimeMat, const ProtocolParams& protParams);

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

