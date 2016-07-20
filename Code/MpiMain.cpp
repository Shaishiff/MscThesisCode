
#include "GA.h"
#include "Log.h"
#include "Protocols.h"
#include "SBModel.h"
#include "Mat.h"
#include <iostream>
#include <fstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

//#define TESTING

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#ifdef TESTING
void StartTestingProcess()
{
	bool bLogToFileOnly = false;
	CGA ga;
	ga.Test();
	LOG("Ending process");
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartMainGaProcess()
{
	bool bLogToFileOnly = false;
	LOG("Starting main process");

	// This is for the controlling master process.
	srand((unsigned int)time(NULL));

	// Start timing.
	clock_t mainStartingTime = clock();

	// Run main program.
	vector<char*> vecTargets;
	char targetName[FILE_NAME_BUFFER_SIZE] = {0};

	//#define NO_NOISE
	//#define NO_DOWN_SAMPLING
	//sprintf(targetName, "TargetFibroblast_(66_66)-(85_85)"); vecTargets.push_back(targetName); // Single square.
	//sprintf(targetName, "TargetFibroblast_(36_36)-(55_55)_+_(86_86)-(105_105)"); vecTargets.push_back(targetName); // Two squares.
	//sprintf(targetName, "TargetFibroblast_(36_36)-(55_55)_+_(36_76)-(55_115)_+_(80_35)-(100_65)"); vecTargets.push_back(targetName); // Three squares.
	//sprintf(targetName, "TargetFibroblast_(90_90)-R=15"); vecTargets.push_back(targetName); // Disc.

	// Target #1
	//sprintf(targetName, "TargetFibroblast_(60_35)-(80_60)_+_(90_90)-(110_110)"); vecTargets.push_back(targetName);

	// Target #2
	sprintf(targetName, "TargetFibroblast_(60_60)-R=20"); vecTargets.push_back(targetName);

	// Target #3
	//sprintf(targetName, "TargetFibroblast_(61_61)-(80_80)_+_(79_79)-(98_98)"); vecTargets.push_back(targetName);

	for(int iVecTargets = 0; iVecTargets < (int)vecTargets.size(); ++iVecTargets)
	{
		#ifndef NO_NOISE
		vector<int> vecNoise;
		//vecNoise.push_back(-21);
		//vecNoise.push_back(-22);
		//vecNoise.push_back(-23);
		//vecNoise.push_back(-24);
		//vecNoise.push_back(-25);
		vecNoise.push_back(-26);
		//for(int iVecNoise = -22; iVecNoise <= -26; iVecNoise -= 2) vecNoise.push_back(iVecNoise);
		//for(int iVecNoise = 30; iVecNoise <= 50; iVecNoise += 5) vecNoise.push_back(iVecNoise);
		for(int iVecNoise = 0; iVecNoise < (int)vecNoise.size(); ++iVecNoise)
		{
			int nNoise = vecNoise[iVecNoise];
		#else
		{
			int nNoise = 100;
		#endif

			char targetFullPath[FILE_NAME_BUFFER_SIZE] = {0};
			#ifndef NO_NOISE
				sprintf(targetFullPath, "_TargetFibroblastMatResults_withNoise=%d_SNRdb_prot", vecNoise[iVecNoise]);
			#else
				sprintf(targetFullPath, "TargetFibroblastMatResults");
			#endif

			#ifndef NO_DOWN_SAMPLING
			vector<int> vecSamplingIntervals;
			vecSamplingIntervals.push_back(16);
			for(int iVecSample = 0; iVecSample < (int)vecSamplingIntervals.size(); iVecSample++)
			{
				int nSampleInterval = vecSamplingIntervals[iVecSample];
			#else
			{
				int nSampleInterval = 1;
			#endif

				double** pCombinedFibroblastMat = CreateMat();
				for(int iAlgoIndex = 30; iAlgoIndex <= 40; ++iAlgoIndex)
				{
					CGA ga;
					ga.RunGA(vecTargets[iVecTargets], targetFullPath, iAlgoIndex, nSampleInterval, nNoise, pCombinedFibroblastMat);
				}
				//char combinedFibroblastFileName[FILE_NAME_BUFFER_SIZE] = {0};
				//sprintf(combinedFibroblastFileName, "%s/CombinedFibroblasts_samp_%d.txt", LOG_FOLDER, vecSamplingIntervals[iVec]);
				//SaveMatToFileWithFullNameIntFormat(pCombinedFibroblastMat, combinedFibroblastFileName);
				//DestroyMat(pCombinedFibroblastMat);
			}
		}
	}

	int nNumberOfMachines = MPI::COMM_WORLD.Get_size(); // same as MPI_Comm_size(MPI_COMM_WORLD, &nthreads)
	for(int iProcess = 1; iProcess < nNumberOfMachines; ++iProcess)
	{
		MPI_Send(NULL, 0, MPI_DOUBLE, iProcess, MPI_DIE_TAG, MPI::COMM_WORLD);
		LOG1("Sending quit flag to process %d", iProcess);
	}

	// End timing.
	clock_t mainEndingTime = clock();
	double mainRunningTime = (mainEndingTime - mainStartingTime)/double(CLOCKS_PER_SEC);
	LOG1("Main duration: %.3f seconds", mainRunningTime);
	LOG("Ending process");
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartSlaveGaProcess()
{
	bool bLogToFileOnly = true;
	LOG("Starting slave process");

	// This is for all the other processes which are not the master.
	// Create all the vars we will use for data transfer.
	MPI_Status status; // can save resources by using the predefined constant MPI_STATUS_IGNORE as a special value for the status argument.
	S1Protocol s1;
	S2Protocol s2;

	// Start the infinite loop (until the master tells us to quit).
	LOG("Starting loop");
	while(true)
	{
		//CFkModel* pModel = new CFkModel();
		CSBModel* pModel = new CSBModel();
		double** fibroblast_mat = CreateMat();
		double** result_mat = CreateMat();

		LOG("Waiting to receive mat");
		int nRet = MPI_Recv(&(fibroblast_mat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, MPI_MASTER, MPI_ANY_TAG, MPI::COMM_WORLD, &status);
		LOG1("Mat received, res: %d", nRet);

		// Start timing.
		clock_t startingTime = clock();

		if (status.MPI_TAG == MPI_JOB_1_TAG)
		{
			LOG("Executing 1st protocol");
			pModel->ExecuteModel(fibroblast_mat, result_mat, s1);
		}
		else if (status.MPI_TAG == MPI_JOB_2_TAG)
		{
			LOG("Executing 2nd protocol");
			pModel->ExecuteModel(fibroblast_mat, result_mat, s2);
		}
		else if (status.MPI_TAG == MPI_DIE_TAG)
		{
			LOG("Got the die tag. Exiting...");
			break;
		}
		else
		{
			LOG("Got an invalid tag. Aborting !");
			break;
		}

		// End timing.
		clock_t endingTime = clock();
		double runningTime = (endingTime - startingTime)/double(CLOCKS_PER_SEC);

		LOG1("Finished executing protocol after %.3f seconds, sending results.", runningTime);
		MPI_Send(&(result_mat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, MPI_MASTER, MPI_RESULT_TAG, MPI::COMM_WORLD);
		LOG("Results were sent.");

		// Clear up the matrix we use for data transfer.
		delete pModel;
		pModel = NULL;
		DestroyMat(fibroblast_mat);
		DestroyMat(result_mat);

	} // End of loop.
	LOG("Ending process");
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	// Init stuff.
	int nProvided = -1;
	int nMpiIntRet = MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &nProvided);
	srand((unsigned int)time(NULL));

	// Getting general info about MPI.
	int nCpuNameLen = MPI_MAX_PROCESSOR_NAME;
	char sMachineName[MPI_MAX_PROCESSOR_NAME] = {0};
	MPI_Get_processor_name(sMachineName, &nCpuNameLen);
	int nCurProcess = MPI::COMM_WORLD.Get_rank(); // same as MPI_Comm_rank(MPI_COMM_WORLD, &tid)

	// Start the appropriate process: master/slave
	if(MPI::COMM_WORLD.Get_rank() == MPI_MASTER)
	{
		CreateLogFile(nCurProcess, sMachineName);
	#ifdef TESTING
		StartTestingProcess();
	}
	#else
		StartMainGaProcess();
	}
	else
	{
		CreateLogFile(nCurProcess, sMachineName);
		StartSlaveGaProcess();
	}
	#endif // TESTING

	MPI_Finalize();

	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
