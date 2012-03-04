
#include "GA.h"

#ifdef USE_MPI
	#include <mpi.h>
#else
	// For detecting memory leaks: http://msdn.microsoft.com/en-us/library/x98tx3cf.aspx
	#define _CRTDBG_MAP_ALLOC
	#define _CRT_SECURE_NO_WARNINGS
	#include <crtdbg.h>
#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef USE_MPI
#define LOG_FILE_NAME "Log.txt"
#else
#define LOG_FILE_NAME "C:\\Users\\Shai\\Documents\\MATLAB\\MSc\\FkModel\\Output\\Log.txt"
#endif

FILE* pLogFile;
char strLog[1024];
#define LOG \
	printf("%s", strLog); \
	pLogFile = fopen(LOG_FILE_NAME, "a"); \
	if(pLogFile != NULL){ fprintf(pLogFile,"%s", strLog); fclose(pLogFile); }

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int Ga::nCurIteration = 0;
double* Ga::pGlobalMeasurement1 = NULL;
double* Ga::pGlobalMeasurement2 = NULL;
double** Ga::pTempMeasurement1 = NULL;
double** Ga::pTempMeasurement2 = NULL;
#ifndef USE_MPI
HANDLE Ga::ghEvents[];
#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::InitGa()
{
	pGlobalMeasurement1 = new double[Nw];
	pGlobalMeasurement2 = new double[Nh];
	pTempMeasurement1 = new double*[Npop];
	pTempMeasurement2 = new double*[Npop];

	for(int iPop = 0; iPop < Npop; ++iPop)
	{
		#ifndef USE_MPI
		ghEvents[iPop] = NULL;
		#endif
		pTempMeasurement1[iPop] = new double[Nw];
		pTempMeasurement2[iPop] = new double[Nh];

		Candidate* pCandidate = new Candidate(iPop);
		Population.push_back(pCandidate);

		char fileName[FILE_NAME_BUFFER_SIZE] = {0};
		sprintf(fileName, "OriginalCandidate_%d.txt", iPop);
		CFkModel::SaveToInputFile(pCandidate->m_mat, fileName);
	}

	MinCost = new double[MaxIterations];
	for(int iteration=0; iteration < MaxIterations; ++iteration)
	{
		MinCost[iteration] = 0.0;
	}
	
	int nRankSum = 0;
	for(int iRankSum = 1; iRankSum <= NsurvivingPopulation; ++iRankSum)
	{
		nRankSum += iRankSum;
	}

	Rank = new double[NsurvivingPopulation];
	for (int iRank = 1 ; iRank <= NsurvivingPopulation; ++iRank)
	{
		Rank[iRank-1] = double(NsurvivingPopulation - iRank + 1)/double(nRankSum);
		if(iRank != 1)
		{
			Rank[iRank-1] = Rank[iRank-1] + Rank[iRank-2];
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::Measure(double* pMeasurement1, double* pMeasurement2, Candidate* pCandidate)
{
	S1Protocol s1;
	S2Protocol s2;
	double** pTempResult1 = NULL;
	double** pTempResult2 = NULL;
	double** pFibroblasts = NULL;

	pTempResult1 = pCandidate->m_pModel->ExecuteFkModel(pCandidate->m_mat, s1);
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		for(int iW = 0; iW < Nw+2; ++iW)
		{											
			pCandidate->m_pResult1[iH][iW] = pTempResult1[iH][iW];			
		}
	}
	for (int iW = 1; iW < Nw+1; ++iW)
	{
		pMeasurement1[iW-1] = pCandidate->m_pResult1[Nh][iW];
	}

	pTempResult2 = pCandidate->m_pModel->ExecuteFkModel(pCandidate->m_mat, s2);
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		for(int iW = 0; iW < Nw+2; ++iW)
		{											
			pCandidate->m_pResult2[iH][iW] = pTempResult2[iH][iW];			
		}
	}
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		pMeasurement2[iH-1] = pCandidate->m_pResult2[iH][Nw];
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double Ga::CalculateCost(Candidate* pCandidate)
{	
	if((nCurIteration == 0) || (pCandidate->m_nIndex != 0))
	{
		Measure(pTempMeasurement1[pCandidate->m_nIndex], pTempMeasurement2[pCandidate->m_nIndex], pCandidate);
		for (int iW = 1; iW < Nw+1; ++iW)
		{
			pCandidate->m_cost += std::abs((long int)(pTempMeasurement1[pCandidate->m_nIndex][iW-1] - pGlobalMeasurement1[iW-1]));
		}

		for (int iH = 1; iH < Nh+1; ++iH)
		{
			pCandidate->m_cost += std::abs((long int)(pTempMeasurement2[pCandidate->m_nIndex][iH-1] - pGlobalMeasurement2[iH-1]));
		}
	}
	return pCandidate->m_cost;
}

#ifndef USE_MPI

DWORD WINAPI MyControllingFunction(LPVOID pParam)
{
	Candidate* pCandidate = (Candidate*)pParam;		
	double dCost = Ga::CalculateCost(pCandidate);	
	if(!SetEvent(Ga::ghEvents[pCandidate->m_nIndex])) 
    {
       printf("ERROR - SetEvent failed (%d)\n", GetLastError());
       return 1;
    }	
	return 0;
}

void Ga::CalculateCosts()
{
	for (int nPopIndex=0; nPopIndex < Npop; ++nPopIndex)
	{	
		sprintf(strLog,"Calculating cost for candidate #%d: %s\n", Population[nPopIndex]->m_nIndex, Population[nPopIndex]->m_str);
		LOG
		HANDLE hThread = CreateThread(NULL, 0, MyControllingFunction, (LPVOID)(Population[nPopIndex]), 0, NULL);
		if(hThread == NULL)
		{
			sprintf(strLog,"CreateThread error: %d\n", GetLastError());
			LOG
			ExitProcess(0);
		}
	}

	DWORD dwWaitRet = WaitForMultipleObjects(Npop, ghEvents, TRUE, 60000);
	sprintf(strLog,"Finished waiting for threads, result: %d\n", dwWaitRet);
	LOG
	
	for (int nPopIndex=0; nPopIndex < Npop; ++nPopIndex)
	{	
		sprintf(strLog,"Calculated cost for candidate #%d: %s, Cost: %.3f\n", Population[nPopIndex]->m_nIndex, Population[nPopIndex]->m_str, Population[nPopIndex]->m_cost);
		LOG
	}
}
#else

void Ga::CalculateCosts()
{
}

#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int Ga::GetMate()
{
	double parentRand = (double)rand()/(double)RAND_MAX;        
    for (int iRank = 0 ; iRank < (NsurvivingPopulation-1); ++iRank)
	{
        if(Rank[iRank] >= parentRand)
		{
            return iRank;
		}
	}
	return (NsurvivingPopulation-1);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CreateChild(Candidate* pParent1, Candidate* pParent2, Candidate* pChild)
{
	if((pChild == pParent1) || (pChild == pParent2))
	{
		throw;
	}
		
	sprintf(strLog,"Parent1, Candidate #%d: %s\n", pParent1->m_nIndex, pParent1->m_str);
	LOG
	sprintf(strLog,"Parent2, Candidate #%d: %s\n", pParent2->m_nIndex, pParent2->m_str);
	LOG

	pChild->m_nCenterH = pParent1->m_nCenterH;
	pChild->m_nCenterW = pParent1->m_nCenterW;
	pChild->m_nHeight = pParent1->m_nHeight;
	pChild->m_nWidth = pParent1->m_nWidth;
	
	int nCrossPoint1 = rand()%4;
	switch(nCrossPoint1)
	{
		case 0:
			pChild->m_nCenterH = pParent2->m_nCenterH;	
			break;
		case 1:
			pChild->m_nCenterW = pParent2->m_nCenterW;
			break;
		case 2:
			pChild->m_nHeight = pParent2->m_nHeight;			
			break;
		case 3:
			pChild->m_nWidth = pParent2->m_nWidth;
			break;
	}

	int nCrossPoint2 = rand()%4;
	switch(nCrossPoint2)
	{
		case 0:
			pChild->m_nCenterH = pParent2->m_nCenterH;	
			break;
		case 1:
			pChild->m_nCenterW = pParent2->m_nCenterW;
			break;
		case 2:
			pChild->m_nHeight = pParent2->m_nHeight;			
			break;
		case 3:
			pChild->m_nWidth = pParent2->m_nWidth;
			break;
	}

	pChild->CreateFibroblastPatch();
	sprintf(strLog,"Child,   Candidate #%d: %s\n", pChild->m_nIndex, pChild->m_str);	
	LOG
	sprintf(strLog,"-\n");
	LOG
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool Ga::FindSimilarCandidate(int nIndex)
{
	for(int iPop = 0; iPop < Npop; ++iPop)
	{
		if(iPop != nIndex)
		{
			if(Population[iPop]->m_nCenterH == Population[nIndex]->m_nCenterH &&
			   Population[iPop]->m_nCenterW == Population[nIndex]->m_nCenterW &&
			   Population[iPop]->m_nHeight == Population[nIndex]->m_nHeight &&
			   Population[iPop]->m_nWidth == Population[nIndex]->m_nWidth)
			{
				sprintf(strLog,"Candiates are identical, #%d & #%d\n", nIndex, iPop);
				LOG
				return true;
			}
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CreateMutations()
{
	sprintf(strLog,"Starting mutations...\n");
	LOG
	sprintf(strLog,"--\n");
	LOG
	for(int iPop = 1; iPop < Npop; ++iPop)
	{
		sprintf(strLog,"Candiate #%d before mutation: %s\n", iPop, Population[iPop]->m_str);
		LOG
		Population[iPop]->Mutate();
		while(FindSimilarCandidate(iPop))
		{
			Population[iPop]->UnMutate();
			Population[iPop]->Mutate();
		}
		Population[iPop]->CreateFibroblastPatch();
		Population[iPop]->m_cost = 0;
		sprintf(strLog,"Candiate #%d after mutation : %s\n", iPop, Population[iPop]->m_str);
		LOG
		sprintf(strLog,"-\n");
		LOG
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef USE_MPI

void Ga::CreateEvents()
{
	for(int iPop = 0; iPop < Npop; ++iPop)
	{
		ghEvents[iPop] = CreateEvent( 
				NULL,   // default security attributes
				FALSE,  // auto-reset event object
				FALSE,  // initial state is nonsignaled
				NULL);  // unnamed object

		if(ghEvents[iPop] == NULL) 
		{ 
			sprintf(strLog,"CreateEvent error: %d\n", GetLastError());
			LOG
			ExitProcess(0); 
		} 
	}
}

void Ga::CloseEvents()
{
	for(int iPop = 0; iPop < Npop; ++iPop)
	{
		CloseHandle(ghEvents[iPop]);
	}
}

#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::RunGa()
{
	sprintf(strLog,"\nInit GA...\n");
	LOG
	InitGa();

	// Create and measure our target mat.
	Candidate* pTargetCandidate = new Candidate(-1);	
	Measure(pGlobalMeasurement1, pGlobalMeasurement2, pTargetCandidate);
	CFkModel::SaveToInputFile(pTargetCandidate->m_pResult1, "TargetRiseTime1.txt");
	CFkModel::SaveToInputFile(pTargetCandidate->m_pResult2, "TargetRiseTime2.txt");
	sprintf(strLog,"Target: %s\n", pTargetCandidate->m_str);
	LOG

	sprintf(strLog,"------------------\n");	
	LOG
	while(nCurIteration <= MaxIterations)
	{
		sprintf(strLog,"Starting iteration #%d\n", nCurIteration);
		LOG
		sprintf(strLog,"------------------\n");
		LOG
		clock_t iterationStartingTime = clock();
	
		#ifndef USE_MPI
		CreateEvents();
		#endif

		// Determine the cost function for each member of the population.
		CalculateCosts();		

		#ifndef USE_MPI
		CloseEvents();
		#endif

		// Sort according to the cost and then mate according to the rank.
		std::sort(Population.begin(), Population.end(), CandidateCompare);
		for(int iPop = 0; iPop < Npop; ++iPop)
		{
			Population[iPop]->m_nIndex = iPop;
		}

		char bestCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
		sprintf(bestCandidateFileName, "BestCandidate_%d.txt", nCurIteration);
		CFkModel::SaveToOutputFile(Population[0]->m_mat, bestCandidateFileName);

		double dAvgCost = 0.0;
		for(int iPop = 0; iPop < Npop; ++iPop)
		{			
			char fileName[FILE_NAME_BUFFER_SIZE] = {0};
			sprintf(fileName, "Candidate_%d_%d.txt", nCurIteration, iPop);
			CFkModel::SaveToInputFile(Population[iPop]->m_mat, fileName);

			char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_1.txt", nCurIteration, iPop);
			CFkModel::SaveToOutputFile(Population[iPop]->m_pResult1, riseTimeCandidateFileName);
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_2.txt", nCurIteration, iPop);
			CFkModel::SaveToOutputFile(Population[iPop]->m_pResult2, riseTimeCandidateFileName);

			dAvgCost += Population[iPop]->m_cost;
		}
		dAvgCost = dAvgCost/Npop;

		// Create the next generation.		
		sprintf(strLog,"--\n");
		LOG
		sprintf(strLog,"Creating the next generation...\n");
		LOG
		sprintf(strLog,"--\n");
		LOG
		int iOffspring = 0;
		while(iOffspring < Nmates)
		{	
			if(iOffspring == 0)
			{
				Candidate* Parent1 = Population[0];
				Candidate* Parent2 = Population[0];
				CreateChild(Parent1, Parent2, Population[NsurvivingPopulation + iOffspring]);
				++iOffspring;
			}
			else
			{
				int nFirstParent = GetMate();
				int nSecondParent = GetMate();
				if(nFirstParent != nSecondParent)
				{
					Candidate* Parent1 = Population[nFirstParent];
					Candidate* Parent2 = Population[nSecondParent];
				
					// Choose a cross over point.
					// Add the new offspirng to the population.
					CreateChild(Parent1, Parent2, Population[NsurvivingPopulation + iOffspring]);
					++iOffspring;
				}
			}
		}
    		         
		CreateMutations();		
		++nCurIteration;

		MinCost[nCurIteration] = Population[0]->m_cost;
		sprintf(strLog,"Avg Cost: %.3f\n", dAvgCost);
		LOG
			sprintf(strLog,"Min Cost: %.3f for %s\n", MinCost[nCurIteration], Population[0]->m_str);
		LOG
		if(MinCost[nCurIteration] <= 0.0)
		{
			break;
		}

		clock_t iterationEndingTime = clock();
		double iterationRunningTime = (iterationEndingTime - iterationStartingTime)/double(CLOCKS_PER_SEC);
		sprintf(strLog,"Iteration duration: %.3f seconds\n", iterationRunningTime);
		LOG
		sprintf(strLog,"------------------\n");
		LOG
	}

	delete pTargetCandidate;
	pTargetCandidate = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef USE_MPI
int _tmain(int argc, _TCHAR* argv[])
#else
int main(int argc, char *argv[])
#endif
{
	// add in MPI startup routines.
	// 1st: launch the MPI processes on each node.
	#ifdef USE_MPI
	MPI_Init(&argc, &argv);
	#endif
	  
	MPI_Comm_rank(MPI_COMM_WORLD, &tid);
	MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
	cpu_name = (char*)calloc(nCpuNameLen,sizeof(char));
	MPI_Get_processor_name(cpu_name,&nCpuNameLen);
	printf("hello MPI user: from process = %i on machine=%s, of NCPU=%i processes\n", tid, cpu_name, nthreads);

	/*
	srand((unsigned int)time(NULL));
	pLogFile = fopen(LOG_FILE_NAME, "w");
	if(pLogFile != NULL){ fclose(pLogFile); }
	clock_t mainStartingTime = clock();
	Ga ga;
	ga.RunGa();
	clock_t mainEndingTime = clock();
	double mainRunningTime = (mainEndingTime - mainStartingTime)/double(CLOCKS_PER_SEC);
	printf("Main duration: %.3f seconds", mainRunningTime);
	_CrtDumpMemoryLeaks();	
	getchar();
	*/

	#ifdef USE_MPI
	MPI_Finalize();
	#endif

	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////