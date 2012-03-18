
#include <pthread.h>
#include "GA.h"
#include "FkModel.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

FILE* pLogFile;
char strLogFileName[1024];
char strLog[1024];
#define LOG \
	printf("%s", strLog); \
	pLogFile = fopen(strLogFileName, "a"); \
	if(pLogFile != NULL){ fprintf(pLogFile,"%s", strLog); fclose(pLogFile); }

Ga ga;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int Ga::nNumberOfMachines = 0;
int Ga::nCurIteration = 0;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void* ThreadFunction(void* arg)
{
	int* pThreadIndex = (int*)arg;
	ga.JobProcessingThreadFunc(*pThreadIndex);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Ga::Ga()
{
	nCurIteration = 0;
	m_pTargetMeasurement1 = new double[Nw];
	m_pTargetMeasurement2 = new double[Nh];
	MinCost = NULL;
	Rank = NULL;
	for(int iThread = 0; iThread < MAX_NUMBER_OF_THREADS; iThread++)
	{
		m_threadIndexArray[iThread] = iThread;
		m_threadFlagArray[iThread] = true;
		m_threadHandleArray[iThread] = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::InitGa()
{		
	for(int iPop = 0; iPop < Npop; ++iPop)
	{		
		Candidate* pCandidate = new Candidate(iPop);
		Population.push_back(pCandidate);

		char fileName[FILE_NAME_BUFFER_SIZE] = {0};
		sprintf(fileName, "OriginalCandidate_%d.txt", iPop);
		CFkModel::SaveToFile(pCandidate->GetFibroblastMat(), fileName);
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
	
	// Create the threads which will allocate the jobs to the
	// different machines. We start from 1 (and not 0) because
	// 0 is the index of the master machine (which runs the main program).
	sprintf(strLog, "Creating %d threads to communicate with the slave processes\n", nNumberOfMachines-1);
	LOG
    for(int iThread = 1; iThread < nNumberOfMachines; iThread++)
	{		
		int nRet = pthread_create(&m_threadHandleArray[iThread], NULL, ThreadFunction, (void*)(&m_threadIndexArray[iThread]));
	}   
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CreateTargetMeasurements()
{
	S1Protocol s1;
	S2Protocol s2;
	CFkModel* pModel = new CFkModel();	
	Candidate* pTargetCandidate = new Candidate(-1);		
	
	pModel->ExecuteModel(pTargetCandidate->GetFibroblastMat(), pTargetCandidate->m_pResult1, s1);	
	for (int iW = 1; iW < Nw+1; ++iW)
	{
		m_pTargetMeasurement1[iW-1] = pTargetCandidate->m_pResult1[Nh][iW];
	}

	pModel->ExecuteModel(pTargetCandidate->GetFibroblastMat(), pTargetCandidate->m_pResult2, s2);
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		m_pTargetMeasurement2[iH-1] = pTargetCandidate->m_pResult2[iH][Nw];
	}
	
	CFkModel::SaveToFile(pTargetCandidate->m_pResult1, "TargetRiseTime1.txt");
	CFkModel::SaveToFile(pTargetCandidate->m_pResult2, "TargetRiseTime2.txt");
	sprintf(strLog, "Target: %s\n", pTargetCandidate->GetFullName());
	LOG
	
	delete pModel;
	pModel = NULL;
	delete pTargetCandidate;
	pTargetCandidate = NULL;	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::JobProcessingThreadFunc(int nThreadIndex)
{	
	printf("JobProcessingThreadFunc | Started | Thread: %d\n", nThreadIndex);				
	while(ga.GetThreadFlag(nThreadIndex))
	{
		Job* pJob = ga.GetJob();
		if(pJob != NULL)
		{	
			printf("JobProcessingThreadFunc | Read job | Thread: %d\n", nThreadIndex);
						
			// Tell the corresponding machine that we want to start a new job.
			double dFlag = MPI_FLAG_START_JOB;
			printf("JobProcessingThreadFunc | Sending start job flag | Thread: %d\n", nThreadIndex);
			MPI_Send(&dFlag, 1, MPI_DOUBLE, nThreadIndex, MPI_FLAG_MSG_TAG, MPI::COMM_WORLD);
			sleep(1);
			
			// Send the actual matrix to process.
			Candidate* pCandidate = pJob->m_pCandidate;
			double** mat = pCandidate->GetFibroblastMat();
			printf("JobProcessingThreadFunc | Processing candidate %d | Thread: %d\n", pCandidate->m_nIndex, nThreadIndex);
			
			printf("JobProcessingThreadFunc | Sending mat for processing | Thread: %d\n", nThreadIndex);
			MPI_Send(mat, (Nw+2)*(Nh+2), MPI_DOUBLE, nThreadIndex, MPI_JOB_MSG_TAG, MPI::COMM_WORLD);
			sleep(1);
		
			MPI_Status* status = new MPI_Status();
		
			// We wait for the result.
			printf("JobProcessingThreadFunc | Waiting to receive result mat 1 | Thread: %d\n", nThreadIndex);
			MPI_Recv(pCandidate->m_pResult1, (Nw+2)*(Nh+2), MPI_DOUBLE, nThreadIndex, MPI_RESULT_MSG_TAG, MPI::COMM_WORLD, status);
			printf("JobProcessingThreadFunc | Received result mat 1 | Thread: %d\n", nThreadIndex);
			sleep(1);
			
			printf("JobProcessingThreadFunc | Waiting to receive result mat 2 | Thread: %d\n", nThreadIndex);
			MPI_Recv(pCandidate->m_pResult2, (Nw+2)*(Nh+2), MPI_DOUBLE, nThreadIndex, MPI_RESULT_MSG_TAG, MPI::COMM_WORLD, status);
			printf("JobProcessingThreadFunc | Received result mat 2 | Thread: %d\n", nThreadIndex);
			sleep(1);
					
			char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_1.txt", 0, pCandidate->m_nIndex);
			CFkModel::SaveToOutputFile(pCandidate->m_pResult1, riseTimeCandidateFileName);
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_2.txt", 0, pCandidate->m_nIndex);
			CFkModel::SaveToOutputFile(pCandidate->m_pResult2, riseTimeCandidateFileName);
			
			// Calculate the cost for this candidate.
			pCandidate->m_cost = 0;
			
			// Calculate cost from result 1.
			for (int iW = 1; iW < Nw+1; ++iW)
			{
				pCandidate->m_cost += std::abs((long int)(pCandidate->m_pResult1[Nh][iW] - m_pTargetMeasurement1[iW]));
			}

			// Calculate cost from result 2.
			for (int iH = 1; iH < Nh+1; ++iH)
			{
				pCandidate->m_cost += std::abs((long int)(pCandidate->m_pResult2[iH][Nw] - m_pTargetMeasurement2[iH]));
			}
			
			printf("JobProcessingThreadFunc | calculated cost: %.2f, candidate: %d, thread: %d\n", pCandidate->m_cost, pCandidate->m_nIndex, nThreadIndex);
			
			delete status;
			status = NULL;
			
			// Delete this job and move on to the next one.
			delete pJob;
			pJob = NULL;
		}
		usleep(500);
	}
	
	// No more jobs to process.
	printf("JobProcessingThreadFunc | No more jobs, exiting | Thread: %d\n", nThreadIndex);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CalculateCosts()
{
	sprintf(strLog, "CalculatingCosts started\n"); 
	LOG
		
	// Create the jobs.
	for (int nPopIndex = 0; nPopIndex < Npop; ++nPopIndex)
	{			
		// We don't need to recalculate the first candidate because
		// it hasn't been mutated.
		bool bAddJob = (Population[nPopIndex]->m_nIndex != 0);
		
		// But we do need to calculate is cost if this is 
		// the first iteration ever.
		bAddJob |= (nCurIteration == 0);
		
		if(bAddJob)
		{	
			sprintf(strLog, "Creating a job for candidate #%d: %s\n", Population[nPopIndex]->m_nIndex, Population[nPopIndex]->GetFullName());
			LOG
		
			Job* pJob = new Job();
			pJob->m_pCandidate = Population[nPopIndex];
			m_jobVector.AddJob(pJob);
		}
	}
		
	sprintf(strLog, "CalculatingCosts finished\n"); 
	LOG
}

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
		
	sprintf(strLog,"Parent1, Candidate #%d: %s\n", pParent1->m_nIndex, pParent1->GetFullName());
	LOG
	sprintf(strLog,"Parent2, Candidate #%d: %s\n", pParent2->m_nIndex, pParent2->GetFullName());
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

	pChild->CreateFibroblasts();
	sprintf(strLog,"Child,   Candidate #%d: %s\n", pChild->m_nIndex, pChild->GetFullName());	
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
		sprintf(strLog,"Candiate #%d before mutation: %s\n", iPop, Population[iPop]->GetFullName());
		LOG
		Population[iPop]->Mutate();
		while(FindSimilarCandidate(iPop))
		{
			Population[iPop]->UnMutate();
			Population[iPop]->Mutate();
		}
		Population[iPop]->CreateFibroblasts();
		Population[iPop]->m_cost = 0.0;
		sprintf(strLog,"Candiate #%d after mutation : %s\n", iPop, Population[iPop]->GetFullName());
		LOG
		sprintf(strLog,"-\n");
		LOG
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::RunGa()
{
	sprintf(strLog, "\nInit GA...\n");
	LOG
	InitGa();

	// Create and measure our target mat.
	CreateTargetMeasurements();

	sprintf(strLog, "------------------\n");	
	LOG
	while(nCurIteration <= MaxIterations)
	{
		sprintf(strLog, "Starting iteration #%d\n", nCurIteration);
		LOG
		sprintf(strLog, "------------------\n");
		LOG
		clock_t iterationStartingTime = clock();
	
		// Determine the cost function for each member of the population.
		CalculateCosts();

		// Sort according to the cost and then mate according to the rank.
		sprintf(strLog, "Sorting the population...\n");
		LOG
		std::sort(Population.begin(), Population.end(), CandidateCompare);
		for(int iPop = 0; iPop < Npop; ++iPop)
		{
			Population[iPop]->m_nIndex = iPop;
		}

		char bestCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
		sprintf(bestCandidateFileName, "BestCandidate_%d.txt", nCurIteration);
		CFkModel::SaveToFile(Population[0]->GetFibroblastMat(), bestCandidateFileName);

		double dAvgCost = 0.0;
		for(int iPop = 0; iPop < Npop; ++iPop)
		{			
			char fileName[FILE_NAME_BUFFER_SIZE] = {0};
			sprintf(fileName, "Candidate_%d_%d.txt", nCurIteration, iPop);
			CFkModel::SaveToFile(Population[iPop]->GetFibroblastMat(), fileName);

			char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_1.txt", nCurIteration, iPop);
			CFkModel::SaveToFile(Population[iPop]->m_pResult1, riseTimeCandidateFileName);
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_2.txt", nCurIteration, iPop);
			CFkModel::SaveToFile(Population[iPop]->m_pResult2, riseTimeCandidateFileName);

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
		sprintf(strLog,"Min Cost: %.3f for %s\n", MinCost[nCurIteration], Population[0]->GetFullName());
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
		
    // Terminate the threads.
	for(int iThread = 1; iThread < nNumberOfMachines; iThread++)
	{		
		m_threadFlagArray[iThread] = false;
	}		
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartMainProcess(char* sMachineName)
{
	sprintf(strLog, "Starting the main process on %s\n", sMachineName);
	LOG
	
	// This is for the controlling master process.
	srand((unsigned int)time(NULL));
	
	// Start timing.
	clock_t mainStartingTime = clock();
		
	// Run main program.
	ga.RunGa();

	for(int iProcess = 1; iProcess < Ga::nNumberOfMachines; ++iProcess)
	{
		double dFlag = MPI_FLAG_QUIT;
		MPI_Send(&dFlag, 1, MPI_DOUBLE, iProcess, MPI_FLAG_MSG_TAG, MPI::COMM_WORLD);
		sprintf(strLog, "Sending quit flag to process %d\n", iProcess);
		LOG
	}
	
	// End timing.
	clock_t mainEndingTime = clock();
	double mainRunningTime = (mainEndingTime - mainStartingTime)/double(CLOCKS_PER_SEC);
	sprintf(strLog, "Main duration: %.3f seconds\n", mainRunningTime);
	LOG
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartSlaveProcess(int nProcess, char* sMachineName)
{
	sprintf(strLog, "Process %i on %s | Starting\n", nProcess, sMachineName);
	LOG
		
	// This is for all the other processes which are not the master.		
	// Create all the vars we will use for data transfer.
	double dFlag = 0.0;
	MPI_Status* status = new MPI_Status; // can save resources by using the predefined constant MPI_STATUS_IGNORE as a special value for the status argument.	
	S1Protocol s1;
	S2Protocol s2;
	CFkModel* pModel = new CFkModel();	
	double** mat = new double*[Nh+2];
	double** result_mat = new double*[Nh+2];
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		mat[iH] = new double[Nw+2];		
		result_mat[iH] = new double[Nw+2];		
		for(int iW = 0; iW < Nw+2; ++iW)
		{								
			mat[iH][iW] = 0;
			result_mat[iH][iW] = 0;
		}
	}

	// Start the infinite loop (until the master tells us to quit).
	sprintf(strLog, "Process %i on %s | Starting loop\n", nProcess, sMachineName);
	LOG
	while(true)
	{
		sprintf(strLog, "Process %i on %s | Waiting for flag\n", nProcess, sMachineName);
		LOG
		MPI_Recv(&dFlag, 1, MPI_DOUBLE, MPI_MASTER, MPI_FLAG_MSG_TAG, MPI::COMM_WORLD, status);
		if(dFlag == MPI_FLAG_QUIT)
		{
			sprintf(strLog, "Process %i on %s | Got a flag to quit\n", nProcess, sMachineName);
			LOG
			break;
		}
		else if (dFlag == MPI_FLAG_START_JOB)
		{					
			sprintf(strLog, "Process %i on %s | Got a flag to start job\n", nProcess, sMachineName);
			LOG			
			sprintf(strLog, "Process %i on %s | Waiting to receive mat\n", nProcess, sMachineName);
			LOG			
			MPI_Recv(mat, (Nw+2)*(Nh+2), MPI_DOUBLE, MPI_MASTER, MPI_JOB_MSG_TAG, MPI::COMM_WORLD, status);
			sprintf(strLog, "Process %i on %s | Mat received\n", nProcess, sMachineName);
			LOG
			sleep(1);
			
			// First protocol.
			sprintf(strLog, "Process %i on %s | Executing 1st protocol\n", nProcess, sMachineName);
			LOG
			//result_mat = pModel->ExecuteModel(mat, s1);
			sprintf(strLog, "Process %i on %s | Finished 1st protocol, sending results\n", nProcess, sMachineName);
			LOG
			MPI_Send(result_mat, (Nw+2)*(Nh+2), MPI_DOUBLE, MPI_MASTER, MPI_RESULT_MSG_TAG, MPI::COMM_WORLD);
			sprintf(strLog, "Process %i on %s | 1st protocol results were sent\n", nProcess, sMachineName);
			LOG		
			sleep(1);
			
			// Second protocol.
			sprintf(strLog, "Process %i on %s | Executing 2nd protocol\n", nProcess, sMachineName);
			LOG
			//result_mat = pModel->ExecuteModel(mat, s2);
			sprintf(strLog, "Process %i on %s | Finished 2nd protocol, sending results\n", nProcess, sMachineName);
			LOG
			MPI_Send(result_mat, (Nw+2)*(Nh+2), MPI_DOUBLE, MPI_MASTER, MPI_RESULT_MSG_TAG, MPI::COMM_WORLD);
			sprintf(strLog, "Process %i on %s | 2nd protocol results were sent\n", nProcess, sMachineName);
			LOG		
			sleep(1);
		}
		else
		{
			sprintf(strLog, "Got an invalid flag. Process = %i on %s\n", nProcess, sMachineName);
			LOG
		}
	} // End of loop.
	
	// Clear up the matrix we use for data transfer.
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		delete [] mat[iH];
		delete [] result_mat[iH];
	}
	delete [] mat;
	delete [] result_mat;
	delete status;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);	
	srand((unsigned int)time(NULL));	
		
	// Getting general info about MPI.
	int nCpuNameLen = MPI_MAX_PROCESSOR_NAME;
	char sMachineName[MPI_MAX_PROCESSOR_NAME] = {0};
	MPI_Get_processor_name(sMachineName, &nCpuNameLen);
	int nCurProcess = MPI::COMM_WORLD.Get_rank(); // same as MPI_Comm_rank(MPI_COMM_WORLD, &tid)
	Ga::nNumberOfMachines = MPI::COMM_WORLD.Get_size(); // same as MPI_Comm_size(MPI_COMM_WORLD, &nthreads)
	
	// Create the log file name.
	sprintf(strLogFileName, "Log_%i_on_%s.txt", nCurProcess, sMachineName);
	
	// Clear the current log file.
	pLogFile = fopen(strLogFileName, "w");
	if(pLogFile != NULL) { fclose(pLogFile); }
		
	sprintf(strLog, "Starting process = %i on %s, out of %i processes\n", nCurProcess, sMachineName, Ga::nNumberOfMachines);
	LOG
	
	if(nCurProcess == MPI_MASTER)
	{		
		StartMainProcess(sMachineName);		
	}
	else
	{		
		StartSlaveProcess(nCurProcess, sMachineName);
	}
	
	printf("Ending process = %i on %s, out of %i processes\n", nCurProcess, sMachineName, Ga::nNumberOfMachines);
	MPI_Finalize();
		
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// TODO:
// Create threads the same number as the number of slave processes.
// Each thread will take a job from the queue and process it.
// Processing involves: Executing model, saving results and measurements and updating the cost
// of the specific jobs candidate.
// The thread will exit when the queue is empty.
// On exit the thread will signal that it's done and this is how the master process
// will know he can continue on.
// Don't forget that when master process exists he needs to signal all other slave
// processes to exit as well by sending them the appropriate flag.