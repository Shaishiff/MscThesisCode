
#include "defs.h"
#include "Candidate.h"
#include <algorithm>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate::Candidate(int nIndex)
{	
	// Create the mats.
	m_pFibroblastMat = CreateMat();
	m_pResult1 = CreateMat();
	m_pResult2 = CreateMat();
	ClearMat();
	
	// Init the vars.
	m_nIndex = nIndex;
	m_cost = 0;
	m_nCenterH = 0;
	m_nCenterW = 0;
	m_nHeight = 0;
	m_nWidth = 0;	
	
	// Now create the fibroblats.
	if(nIndex == -1)
	{
		// Use the target fibroblast mat.
		m_nCenterH = TargetCenterH;
		m_nCenterW = TargetCenterW;
		m_nHeight = TargetHeight;
		m_nWidth = TargetWidth;	
	}
	else
	{
		// Create a random fibroblast mat.
		m_nCenterH = rand()%int(Nh) + 1;
		m_nCenterW = rand()%int(Nw) + 1;
		m_nHeight = rand()%int(Nh) + 1;
		m_nWidth = rand()%int(Nw) + 1;
	}
	CreateFibroblasts();
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

void Candidate::CreateFibroblasts()
{
	ClearMat();
	int nhStart = max((m_nCenterH),1);
	int nhEnd = min((m_nCenterH+m_nHeight),Nh);
	int nwStart = max((m_nCenterW),1);
	int nwEnd = min((m_nCenterW+m_nWidth),Nw);
	for (int iH=nhStart; iH < nhEnd; ++iH)
	{	
		for (int iW=nwStart; iW < nwEnd; ++iW)
		{
			m_pFibroblastMat[iH][iW] = 1.0;
		}
	}
	sprintf(m_cCandidateFullName,"(%d,%d) %dX%d", m_nCenterH, m_nCenterW, m_nHeight, m_nWidth);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

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
			m_nCenterH += m_nVal1;
			break;
		case 1:
			m_nCenterW += m_nVal1;
			break;
		case 2:
			m_nHeight += m_nVal1;
			break;
		case 3:
			m_nWidth += m_nVal1;
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

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CandidateCompare(Candidate* pFirstCandidate, Candidate* pSecondCandidate)
{
	return (pFirstCandidate->m_cost < pSecondCandidate->m_cost);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

