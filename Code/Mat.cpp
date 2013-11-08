
#include "Mat.h"
#include "Log.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double** CreateMat()
{
	double* data = (double *)malloc(Nh_with_border*Nw_with_border*sizeof(double));
    double** mat = (double **)malloc(Nh_with_border*sizeof(double*));
    for (int i = 0; i < Nh_with_border; i++)
	{
        mat[i] = &(data[Nw_with_border*i]);
	}
	
	for (int iH = 0; iH < Nh_with_border; ++iH)
	{
		for (int iW = 0; iW < Nw_with_border; ++iW)
		{
			mat[iH][iW] = 0.0;
		}
	}
	
	return mat;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void DestroyMat(double** mat)
{
	free(mat[0]);
	free(mat);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void PrintMat(double** mat)
{
	printf("Printing mat\n"); fflush(stdout);
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		for (int iW = 1; iW < Nw+1; ++iW)
		{
			printf("%d,%d: %f\n", iH, iW, mat[iH][iW]); 
			fflush(stdout);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool SaveMatToFileWithFullNameIntFormat(double** mat, char* fullFileName)
{
	FILE* pFile = fopen(fullFileName, "w");
	if(pFile == NULL)
	{
		return false;
	}
	
	for (int iH = 0; iH <= Nh+1; ++iH)
	{
		for (int iW = 0; iW <= Nw+1; ++iW)
		{
			fprintf(pFile, "%d", (int)ceil(mat[iH][iW]));
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
	return true;
}

bool SaveMatToFileWithFullName(double** mat, char* fullFileName)
{
	FILE* pFile = fopen(fullFileName, "w");
	if(pFile == NULL)
	{
		return false;
	}
	
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		for (int iW = 1; iW < Nw+1; ++iW)
		{
			fprintf(pFile, "%4.3f ", mat[iH][iW]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
	return true;
}

bool SaveMatToFile(double** mat, char* fileName)
{
	return SaveMatToFile(mat, fileName, NULL);
}

bool SaveMatToFile(double** mat, char* fileName, char* folder)
{
	char fullFileName[1024] = {0};
	if(folder == NULL)
	{
		sprintf(fullFileName, "%s/%s", LOG_FOLDER, fileName);
	}
	else
	{
		sprintf(fullFileName, "%s/%s", folder, fileName);
	}
	return SaveMatToFileWithFullName(mat, fullFileName);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
