#ifndef __GA_H__
#define __GA_H__

#include "defs.h"
#include "Candidate.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

class Ga
{
public:
	Ga()
	{
		nCurIteration = 0;
		pGlobalMeasurement1 = NULL;
		pGlobalMeasurement2 = NULL;
		pTempMeasurement1 = NULL;
		pTempMeasurement2 = NULL;		
		MinCost = NULL;
		Rank = NULL;		
	}
	void RunGa();
	static void Measure(double* pMeasurement1, double* pMeasurement2, Candidate* pCandidate);
	static double CalculateCost(Candidate* pCandidate);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef USE_MPI

public:
	static HANDLE ghEvents[Npop];
private:
	void CreateEvents();
	void CloseEvents();

#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

private:
	void InitGa();	
	void CalculateCosts();
	int GetMate();
	void CreateChild(Candidate* pParent1, Candidate* pParent2, Candidate* pChild);
	bool FindSimilarCandidate(int nIndex);
	void CreateMutations();

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
	
private:
	static int nCurIteration;
	static double* pGlobalMeasurement1;
	static double* pGlobalMeasurement2;
	static double** pTempMeasurement1;
	static double** pTempMeasurement2;	
	double* MinCost;
	double* Rank;
	vector<Candidate*> Population;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __GA_H__