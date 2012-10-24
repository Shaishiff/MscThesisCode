
#include "Log.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

FILE* pLogFile;
char strLogSourceName[1024];
char strLogFileName[1024];
char strLog[1024];

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CreateLogFile(int nCurProcess, char* sMachineName)
{
	// Create the log file name.
	if(nCurProcess != 0)
	{
		sprintf(strLogFileName, "%s/Log_%i_on_%s.txt", LOG_FOLDER, nCurProcess, sMachineName);
	}
	else
	{
		sprintf(strLogFileName, "%s/Log_master_on_%s.txt", LOG_FOLDER, sMachineName);
	}
	sprintf(strLogSourceName, "Process %i on %s | ", nCurProcess, sMachineName);
	
	// Clear the current log file.
	pLogFile = fopen(strLogFileName, "w");
	if(pLogFile != NULL) { fclose(pLogFile); }	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
