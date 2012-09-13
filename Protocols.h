
#ifndef __PROTOCOLS__
#define __PROTOCOLS__

//////////////////////////////////////////////////////////////////////////////////////////
// Protocols
//////////////////////////////////////////////////////////////////////////////////////////

struct ProtocolParams
{
	double m_Amp;
	double m_TotalTime;
	double m_BeginTime;
	double m_hStart;
	double m_hEnd;
	double m_wStart;
	double m_wEnd;
	double m_hMeasureStart;
	double m_hMeasureEnd;
	double m_wMeasureStart;
	double m_wMeasureEnd;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// Protocol parameters:
const double S1Amp = 1.5;//5000.0/3000.0;//500.0; // microA/cm^3
const double S1TotalTime = 5; // milliseconds
const double S1BeginTime = 5.0; // milliseconds
//const double S1hStart = 0.0;
//const double S1hEnd = dH*2;
const double S1hStart = dH*MeasurementMarginIndexes;
const double S1hEnd = dH*(MeasurementMarginIndexes+2);
const double S1wStart = 0.0;
const double S1wEnd = W;
const double S1hMeasureStart = H-dH*MeasurementMarginIndexes;
const double S1hMeasureEnd = H-dH*MeasurementMarginIndexes;
const double S1wMeasureStart = 0.0;
const double S1wMeasureEnd = W;

struct S1Protocol : public ProtocolParams
{
	S1Protocol()
	{
		m_Amp = S1Amp;
		m_TotalTime = S1TotalTime;
		m_BeginTime = S1BeginTime;
		m_hStart = S1hStart;
		m_hEnd = S1hEnd;
		m_wStart = S1wStart;
		m_wEnd = S1wEnd;
		m_hMeasureStart = S1hMeasureStart;
		m_hMeasureEnd = S1hMeasureEnd;
		m_wMeasureStart = S1wMeasureStart;
		m_wMeasureEnd = S1wMeasureEnd;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

const double S2Amp = 1.5;//5000.0/3000.0;//500.0; // microA/cm^3
const double S2TotalTime = 5; // milliseconds
const double S2BeginTime = 5.0; // milliseconds
const double S2hStart = 0.0;
const double S2hEnd = H;
//const double S2wStart = 0.0;
//const double S2wEnd = dW*2;
const double S2wStart = dW*MeasurementMarginIndexes;
const double S2wEnd = dW*(MeasurementMarginIndexes+2);
const double S2hMeasureStart = 0.0;
const double S2hMeasureEnd = H;
const double S2wMeasureStart = W-dW*MeasurementMarginIndexes;
const double S2wMeasureEnd = W-dW*MeasurementMarginIndexes;

struct S2Protocol : public ProtocolParams
{
	S2Protocol()
	{
		m_Amp = S2Amp;
		m_TotalTime = S2TotalTime;
		m_BeginTime = S2BeginTime;
		m_hStart = S2hStart;
		m_hEnd = S2hEnd;
		m_wStart = S2wStart;
		m_wEnd = S2wEnd;
		m_hMeasureStart = S2hMeasureStart;
		m_hMeasureEnd = S2hMeasureEnd;
		m_wMeasureStart = S2wMeasureStart;
		m_wMeasureEnd = S2wMeasureEnd;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif //__PROTOCOLS__
