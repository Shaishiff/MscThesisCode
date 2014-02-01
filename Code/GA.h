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
	void RunGA(char* target, char* targetResults, int iAlgoIndex, int nSamplingInterval, int nNoise, double** pCombinedFibroblastMat);
	void Test();

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

private:
	bool ReadIntoArray(ifstream& myfile, double** arr);
	bool ReadTargetMeasurements(char* target, char* targetResults);
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
	int m_nNumberOfSlaveMachines;
	#define MAX_SLAVE_MACHINES	100
	MPI_Request m_mpiRequestArray[MAX_SLAVE_MACHINES];
	MPI_Status m_mpiStatusArray[MAX_SLAVE_MACHINES];
	double* m_pTargetMeasurement1;
	double* m_pTargetMeasurement2;
	double** m_pTargetFibroblastMat;
	unsigned long int* MinCost;
	double* Rank;
	vector<Candidate*> m_population;
	vector<Job> m_jobVector;
	vector<char*> PastCandidates;
	bool bLogToFileOnly;
	int m_nSamplingInterval;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __GA_H__

