
#ifndef __FK_MODEL_DEFS__
#define __FK_MODEL_DEFS__

//////////////////////////////////////////////////////////////////////////////////////////
// Physiological consts
//////////////////////////////////////////////////////////////////////////////////////////

const double Cm = 1.0; // membrance capacitance - microF/cm^2
const double Am = 3000.0; // 1/cm
const double sigma = 0.5;//8.3;//7.75; // 0.3; 3.0*10^-4; mS/cm

//////////////////////////////////////////////////////////////////////////////////////////
// APD & CV Restitution
//////////////////////////////////////////////////////////////////////////////////////////

const double t_fast  = 0.249; // tau_fast
const double tf_opn1 = 40; // tau_f_open_1
const double tf_opn2 = 82.5; // tau_f_open_2
const double tf_cls = 5.75*1.6; // tau_f_close
const double t_slow = 226.9; // tau_slow
const double ts_opn = 100*0.2; // tau_s_open
const double ts_cls = 300*1.6; // tau_s_close
const double t_ug1 = 64.7*1.15; // tau_ung_1
const double t_ug2 = 222.9*1.15; // tau_ung_2
const double Vrest = -81.2; // in  mV
const double Vpeak = 3.6; // in  mV
const double V_f = 0.16;  
const double V_s = V_f;
const double V_u = V_f;
const double Vf_open = 0.04;
const double V_slow  =  0.85;
const double kwm = 10;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif //__FK_MODEL_DEFS__
