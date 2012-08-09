/* constants for the control human action potential*/
/* As defined in the Nattel paper ? FULL REF ? */
/* Nconst_chol_la.h, for LA*/

double Cm = 100.0; /* area of cell decreased by 4 from 100*100 to 50*50 um^2, was 100.0 pF */
double dt = 0.005;              /*ms*/
double R = 8.3143;             /*J/mol*K*/
double T = 310.0;              /*K*/
double F = 96.4867;            /*C/mmol*/
double RTF = (R*T/F);	       /* A widely used combination */
double Vcell = 20100.0;        /*um^3*/
double Vi = 13668.0;           /*um^3*/
double Vup = 1109.52;          /*um^3*/
double Vrel = 96.48;           /*um^3*/
double Ko = 5.4;               /*mM*/
double Nao = 140;              /*mM*/
double Cao = 1.8;              /*mM*/
double gna = 7.8*Cm;           /*nS/pF*/

//double gk1 = 0.58*0.09*Cm;          /*nS/pF  multiply by 0.58 for PV  multiply by 4 as ytest*/
//double gto = 0.5*0.75*0.1652*Cm;                /*nS/pF  multiply by 0.75 for PV  multiplied by 0.5 for CHF from Li et al Circulation 2000  */
//double gkr = 1.5*0.0294*Cm;                /*nS/pF  try to multiply by 1.5 for PV  */
//double gks = 0.7*1.6*0.129*Cm;         /*nS/pF   multiply by 1.6 for PV, multiplied by 0.7 for CHF from Li et al Circulation 2000   */
//double gcal = 0.7*0.35*0.1238*Cm;               /*nS/pF multiply by 0.7 for PV  multiplied by 0.7 for CHF from Li et al Circulation 2000 */

// double gk1 = 0.09*Cm;          /*nS/pF  multiply by 0.58 for PV  multiply by 4 as ytest*/
// double gto = 0.1652*Cm;                /*nS/pF  multiply by 0.75 for PV  multiplied by 0.5 for CHF from Li et al Circulation 2000  */
// double gkr = 0.0294*Cm;                /*nS/pF  try to multiply by 1.5 for PV  */
// double gks = 0.129*Cm;         /*nS/pF   multiply by 1.6 for PV, multiplied by 0.7 for CHF from Li et al Circulation 2000   */
// double gcal = 0.1238*Cm;               /*nS/pF multiply by 0.7 for PV  multiplied by 0.7 for CHF from Li et al Circulation 2000 */




double gk1 = 0.09*Cm;          /*nS/pF*/
double gto = 0.1652*Cm;                /*nS/pF multiplied by 0.5 for CHF from Li et al Circulation 2000 */
double gkr = 0.0294*Cm;                /*nS/pF*/
double gks = 0.129*Cm;         /*nS/pF multiplied by 0.7 for CHF from Li et al Circulation 2000 */
double gcal = 0.1238*Cm;               /*nS/pF multiplied by 0.7 for CHF from Li et al Circulation 2000 */


double gbca = 0.00113*Cm;      /*nS/pF*/
double gbna = 0.000674*Cm;     /*nS/pF*/
double gkur_sub = 1.0*Cm;     /*nS/pF*/
/*the value of actual gkur is defined in the main program*/
double Inakmax = 0.6*Cm;               /*pA/pF*/
//double Inacamax = 1600.0*Cm;           /*pA/pF multiplied by 1.45 for CHF from Li et al Circulation 2000 */

double Inacamax = 1600.0*Cm;           /*pA/pF multiplied by 1.45 for CHF from Li et al Circulation 2000 */
double Ipcamax = 0.275*Cm;     /*pA/pF*/
double Iupmax = 0.005;         /*mM/ms*/
double Kq10 = 3.0;
double lambda = 0.35;
double Kmnai = 10.0;           /*mM*/
double Kmko = 1.5;             /*mM*/
double Kmna = 87.5;            /*mM*/
double Kmca = 1.38;            /*mM*/
double ksat = 0.1;
double krel = 30.0;            /*ms^-1*/
double kup = 0.00092;          /*mM*/
double Caupmax = 15.0;         /*mM*/
double Cmdnmax = 0.05;         /*mM*/
double Trpnmax = 0.07;         /*mM*/
double Csqnmax = 10.0;         /*mM*/
double KmCmdn = 0.00238;       /*mM*/
double KmTrpn = 0.0005;        /*mM*/
double KmCsqn = 0.8;           /*mM*/
double taufca = 2.0;
double tautr = 180.0;
double tauu = 8.0;
double Ach = 0.01;              /*uM*/
double gkach = 2.0*Cm;           /*nS/pF*/
