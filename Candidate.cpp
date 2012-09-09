
#include "defs.h"
#include "Candidate.h"
#include "Mat.h"
#include <algorithm>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate::Candidate(int nIndex, int nHStart, int nWStart, int nHEnd, int nWEnd)
{
	// Create the mats.
	m_pFibroblastMat = CreateMat();	
	m_pResult1 = CreateMat();
	m_pResult2 = CreateMat();
		
	// Init the vars.
	m_nIndex = nIndex;	
	m_nHStart = nHStart;
	m_nWStart = nWStart;
	m_nHEnd = nHEnd;
	m_nWEnd = nWEnd;
	m_cost = 0;
	/*
	m_nParam1 = 0;
	m_nParam2 = 0;
	m_nVal1 = 0;
	m_nVal2 = 0;
	*/
	memset(m_cCandidateFullName, 0, sizeof(m_cCandidateFullName));
	
	// Now create the fibroblats.
	ClearMat();
	CreateFibroblastMat();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate::~Candidate()
{
	DestroyMat(m_pFibroblastMat);
	DestroyMat(m_pResult1);
	DestroyMat(m_pResult2);
}
	
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::ClearMat()
{
	for (int iH = 0; iH < Nh_with_border; ++iH)
	{
		for (int iW = 0; iW < Nw_with_border; ++iW)
		{								
			if((iH == 0) || (iH == Nh+1) || (iW == 0) || (iW == Nw+1))
			{
				m_pFibroblastMat[iH][iW] = 1.0;
			}
			else
			{
				m_pFibroblastMat[iH][iW] = 0.0;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::CreateFibroblastMat()
{
	//ClearMat();
	for (int iH = m_nHStart; iH <= m_nHEnd; ++iH)
	{	
		for (int iW = m_nWStart; iW <= m_nWEnd; ++iW)
		{
			m_pFibroblastMat[iH][iW] = 1.0;
		}
	}	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
/*
void Candidate::Mutate()
{
	m_nParam1 = rand()%4;
	m_nVal1 = rand()%9 - 4;
	m_nParam2 = rand()%4;
	while(m_nParam2 == m_nParam1) { m_nParam2 = rand()%4; }
	m_nVal2 = rand()%9 - 4;

	switch(m_nParam1)
	{
		case 0:
			m_nHStart += m_nVal1;
			break;
		case 1:
			m_nWStart += m_nVal1;
			break;
		case 2:
			m_nHEnd += m_nVal1;
			break;
		case 3:
			m_nWEnd += m_nVal1;
			break;
	}		
	switch(m_nParam2)
	{
		case 0:
			m_nCenterH += m_nVal2;
			break;
		case 1:
			m_nCenterW += m_nVal2;
			break;
		case 2:
			m_nHeight += m_nVal2;
			break;
		case 3:
			m_nWidth += m_nVal2;
			break;
	}

	m_nCenterH = min(max(m_nCenterH, 1), Nh);
	m_nCenterW = min(max(m_nCenterW, 1), Nw);
	m_nHeight = min(max(m_nHeight, 1), Nh);
	m_nWidth = min(max(m_nWidth, 1), Nw);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::UnMutate()
{
	switch(m_nParam1)
	{
		case 0:
			m_nCenterH -= m_nVal1;
			break;
		case 1:
			m_nCenterW -= m_nVal1;
			break;
		case 2:
			m_nHeight -= m_nVal1;
			break;
		case 3:
			m_nWidth -= m_nVal1;
			break;
	}		
	switch(m_nParam2)
	{
		case 0:
			m_nCenterH -= m_nVal2;
			break;
		case 1:
			m_nCenterW -= m_nVal2;
			break;
		case 2:
			m_nHeight -= m_nVal2;
			break;
		case 3:
			m_nWidth -= m_nVal2;
			break;
	}
}
*/
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

char* Candidate::GetFullName()
{ 
	sprintf(m_cCandidateFullName,"(%d,%d) -> (%d,%d)", m_nHStart, m_nWStart, m_nHEnd, m_nWEnd);
	return m_cCandidateFullName; 
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CandidateCompare(Candidate* pFirstCandidate, Candidate* pSecondCandidate)
{
	return (pFirstCandidate->m_cost < pSecondCandidate->m_cost);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

