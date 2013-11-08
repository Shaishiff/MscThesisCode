#ifndef __GA_H__
#define __GA_H__

#include <mpi.h>
#include "safeJobVector.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

class CGA
{
public:
	CGA();
	~CGA();
	void RunGA(int iAlgoIndex, double** pCombinedFibroblastMat);
	void Test();
	
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

private:
	bool ReadTargetMeasurements();
	void CreateJobs();
	void ProcessJobs();
	void ProcessResults();
	void ProcessResults(Candidate* pCandidate);
	int ChooseRandomParent();
	int Mutate(int nInValue);
	void CreateRandomPopulation();
	void CreateNextGeneration();
	Candidate* CreateChild(Candidate* pParent1, Candidate* pParent2, int nIndex);
	Candidate* CreateRandomCandidate(int nIndex);	
	double CalculateTargetCoverage(double** pBestMatch);
	int CalculateError(double** pBestMatch);
	bool IsNewChild(Candidate* pChild);
	void AddChildToPastCandidates(Candidate* pChild);
	void ClearPastCandidates();
	void ClearPopulation();
	void AddSortedToPatchesVector(FibroblastPatchVector& vecFibroblastPatchVector, FibroblastPatch newFibroblastPatch);
	double CalculatePatchDistance(const FibroblastPatch& newFibroblastPatch);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
	
private:
	int m_nNumberOfMachines;
	double* m_pTargetMeasurement1;
	double* m_pTargetMeasurement2;
	double** m_pTargetFibroblastMat;
	unsigned long int* MinCost;
	double* Rank;
	vector<Candidate*> m_population;
	CSafeJobVector m_jobVector;
	vector<char*> PastCandidates;
	bool bLogToFileOnly;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __GA_H__

