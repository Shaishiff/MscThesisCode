/* conconst stants for the control human action potential*/
/* As defined in the Nattel paper ? FULL REF ? */
/* Nconst_chol_la.h, for LA*/

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

const double Cm = 100.0; /* area of cell decreased by 4 from 100*100 to 50*50 um^2, was 100.0 pF */
const double dt = 0.005; /*ms*/
const double R = 8.3143; /*J/mol*K*/
const double T = 310.0; /*K*/
const double F = 96.4867; /*C/mmol*/
const double RTF = (R*T/F); /* A widely used combination */
const double Vcell = 20100.0; /*um^3*/
const double Vi = 13668.0; /*um^3*/
const double Vup = 1109.52; /*um^3*/
const double Vrel = 96.48; /*um^3*/
const double Ko = 5.4; /*mM*/
const double Nao = 140; /*mM*/
const double Cao = 1.8; /*mM*/
const double gna = 7.8*Cm; /*nS/pF*/
const double gk1 = 0.09*Cm; /*nS/pF*/
const double gto = 0.1652*Cm; /*nS/pF multiplied by 0.5 for CHF from Li et al Circulation 2000 */
const double gkr = 0.0294*Cm; /*nS/pF*/
const double gks = 0.129*Cm; /*nS/pF multiplied by 0.7 for CHF from Li et al Circulation 2000 */
const double gcal = 0.1238*Cm; /*nS/pF multiplied by 0.7 for CHF from Li et al Circulation 2000 */
const double gbca = 0.00113*Cm; /*nS/pF*/
const double gbna = 0.000674*Cm; /*nS/pF*/
const double gkur_sub = 1.0*Cm; /*nS/pF*/ /*the value of actual gkur is defined in the main program*/
const double Inakmax = 0.6*Cm; /*pA/pF*/
const double Inacamax = 1600.0*Cm; /*pA/pF multiplied by 1.45 for CHF from Li et al Circulation 2000 */
const double Ipcamax = 0.275*Cm; /*pA/pF*/
const double Iupmax = 0.005; /*mM/ms*/
const double Kq10 = 3.0;
const double lambda = 0.35;
const double Kmnai = 10.0; /*mM*/
const double Kmko = 1.5; /*mM*/
const double Kmna = 87.5; /*mM*/
const double Kmca = 1.38; /*mM*/
const double ksat = 0.1;
const double krel = 30.0; /*ms^-1*/
const double kup = 0.00092; /*mM*/
const double Caupmax = 15.0; /*mM*/
const double Cmdnmax = 0.05; /*mM*/
const double Trpnmax = 0.07; /*mM*/
const double Csqnmax = 10.0; /*mM*/
const double KmCmdn = 0.00238; /*mM*/
const double KmTrpn = 0.0005; /*mM*/
const double KmCsqn = 0.8; /*mM*/
const double taufca = 2.0;
const double tautr = 180.0;
const double tauu = 8.0;
const double Ach = 0.01; /*uM*/
const double gkach = 2.0*Cm; /*nS/pF*/

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////