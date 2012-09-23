
#ifndef __LOG__H_
#define __LOG__H_

#include "defs.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#define LOG_FOLDER "/a/home/cc/students/enginer/shaishif/Logs"

#define PRINTLOG printf("%s%s\n", strLogSourceName, strLog); fflush(stdout); \
	pLogFile = fopen(strLogFileName, "a"); \
	if(pLogFile != NULL) \
	{ fprintf(pLogFile,"%d | %s%s\n", clock(), strLogSourceName, strLog); fclose(pLogFile); }

#define LOG(log) 								sprintf(strLog, log); \
												PRINTLOG
#define LOG1(log, var1) 						sprintf(strLog, log, var1); \
												PRINTLOG
#define LOG2(log, var1, var2)					sprintf(strLog, log, var1, var2); \
												PRINTLOG
#define LOG3(log, var1, var2, var3)				sprintf(strLog, log, var1, var2, var3); \
												PRINTLOG										
#define LOG4(log, var1, var2, var3, var4)		sprintf(strLog, log, var1, var2, var3, var4); \
												PRINTLOG
#define LOG5(log, var1, var2, var3, var4, var5)	sprintf(strLog, log, var1, var2, var3, var4, var5); \
												PRINTLOG
											
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
										
#endif // __LOG__H_
