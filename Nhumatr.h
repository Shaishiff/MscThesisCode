 
/********************************************************************************
Action Potential Model of a Human Atrial Cell

humatr.c

Modified to humatr.h to be included in 2d_atria.cc by Omer Berenfeld
on 10/11/02

Last modified by Sandeep V. Pandit
on 10/05/02

Fixed step integration method

Based on Program Written By Joseph Tranquillo, BME, Duke University

From The Journal Article : "Ionic Mechanisms Underlying Human Atrial Action
        Properties: Insights from a Mathematical Model"
Authors: Courtemanche,M. , Ramirez,R. and Nattel,S.
Am. J. Physiol. 275 (Heart Circ. Physiol. 44); pages H301-H321. 1998.

NOTES
1. The constants are in the file "hum_*.h", whose path needs to be set
2. Remember to change the time you want to write to the output file;
   By default it the program writes only the last 1000 seconds of simulation
3. Added Ist to Ki
4. Write the last file at the start of the last second (usually at the 12th second)
5. Sample data to reduce the output file size, also write only AP
6. Change shift time to 20 msec.
7. Also includes equation for modified IK1 FOR kir2.3
	
********************************************************************************/
#include <stdio.h>
#include <math.h>
#include "Nconst_chol_la.h" /* File of constants must be in the same folder */


#define Initial_condition_fname "..\\msc_thesis_c_code\\Ninit_chol_la.d"
const int count_init = 21;


/*Global Variables */
double state[count_init];

/* List of state variables (variables of ODE) */
class state_variables
{
public:
	double V;
	double h;
	double d;
	double xr;
	double Nai;
	double Ki;
	double Carel;
	double si;
	double ui;
	double v;
	double m;
	double j;
	double f;
	double xs;
	double Cai;
	double Caup;
	double sa;
	double ua;
	double fca;
	double u;
	double w;
}; /* of class state_variables */

// Function prototypes

void	Get_state_variables_initial_condition(void);
void	Assign_node_initial_condition(state_variables &rN);
double	Total_transmembrane_currents(state_variables &rN, double &Ist);

FILE *initcond;

/*********************************************************************/

void	Get_state_variables_initial_condition()
{

/*      Reading initial conditions*/
initcond = fopen(Initial_condition_fname,"r");
for(int jj=0; jj<count_init; jj++)
	fscanf(initcond,"%le ",&state[jj]);

fclose(initcond);
} /* of Get_state_variables_initial_condition

/*********************************************************************/

void Assign_node_initial_condition(state_variables &rN)
{

/* Assign Initial Conditions to variables*/
rN.V = state[0];
rN.h = state[1];
rN.d = state[2];
rN.xr = state[3];
rN.Nai = state[4];
rN.Ki = state[5];
rN.Carel = state[6];
rN.si = state[7];
rN.ui = state[8];
rN.v = state[9];
rN.m = state[10];
rN.j = state[11];
rN.f = state[12];
rN.xs = state[13];
rN.Cai = state[14];
rN.Caup = state[15];
rN.sa = state[16];
rN.ua = state[17];
rN.fca = state[18];
rN.u = state[19];
rN.w = state[20];

} /* of Assign_node_initial_condition */


/*********************************************************************/

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/***************************
* Schraudolph's algorithm: *
***************************/
static union
{
	double d;
	struct
	{
		int j, i; // #ifdef LITTLE_ENDIAN
	} n;
} _eco;

#define EXP_A (1048576/M_LN2)	/* 2^20/LN(2) use 1512775 for integer version */
#define EXP_C 60801	 /* see text for choice of c values */
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#if 0
#define shaisExp(in_var) EXP(in_var)
#else 
#define shaisExp(in_var) exp(in_var)
#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

class Model
{
	double Iion,Ina,Ik1,Ito,Ikur,Ikr,Iks,Ical,Iup;
	double Ipca,Inak,Inaca,Ibna,Ibca,Iupleak,Irel,Itr;
	double B1,B2,Ikach,ach2;
	double ach0;
	double ach1;
	double am,bm,aj,bj,ah,bh,taum,tauh,tauj,infm,infh,infj;
	double asa,bsa,asi,bsi,tausa,infsa,tausi,infsi;
	double aui,bui,aua,bua,tauui,infui,tauua,infua;
	double axr,bxr,infxr,tauxr,axs,bxs,tauxs,infxs;
	double taud,infd,tauf,inff,inffca,fnak,sigma;
	double infu,tauv,infv,tauw,infw,Fn;
	double Ena,Ek,Eca;
	double gkur;

public:
	Model()
	{
		ach0 = pow(Ach,0.477811);
		ach1 = gkach*10.0*ach0/(ach0 + 9.13652); 
	}

	double Total_transmembrane_currents(const state_variables &rN, state_variables &outN, double &Ist)
	{
		Ena = RTF*log(Nao/rN.Nai);
		Ek = RTF*log(Ko/rN.Ki);
		Eca = RTF*0.5*log(Cao/rN.Cai);

		/* update currents */
		Ina = gna*rN.m*rN.m*rN.m*rN.h*rN.j*(rN.V-Ena);

		Ik1 = gk1*(rN.V-Ek)/(1+shaisExp(0.07*(rN.V+80.0)));

		/* FOR Kir2.3: Ik1 = gk1*(rN.V-Ek)/(1+shaisExp(0.07*(rN.V+80.0))) + 15.0/(1+exp((rN.V+40.0)/-7.0)); */

		Ito = gto*rN.sa*rN.sa*rN.sa*rN.si*(rN.V-Ek);
		Ikur = gkur*rN.ua*rN.ua*rN.ua*rN.ui*(rN.V-Ek);
		Ikr = gkr*rN.xr*(rN.V-Ek)/(1.0 + exp((rN.V+15.0)/22.4));
		Iks = gks*rN.xs*rN.xs*(rN.V-Ek);
		Ical = gcal*rN.d*rN.f*rN.fca*(rN.V-65.0);
		Inak = Inakmax*fnak*(1.0/(1.0+ pow(Kmnai/rN.Nai,1.5))) * Ko/(Ko+Kmko);
		Inaca = Inacamax*(exp(lambda*rN.V/(RTF))*rN.Nai*rN.Nai*rN.Nai*Cao - exp((lambda-1)*rN.V/(RTF))*Nao*Nao*Nao*rN.Cai)/
		((Kmna*Kmna*Kmna + Nao*Nao*Nao)*(Kmca + Cao)*(1 + ksat*exp((lambda-1)*rN.V/(RTF))));
		Ibca = gbca*(rN.V-Eca);
		Ibna = gbna*(rN.V-Ena);
		Ipca = Ipcamax*rN.Cai/(0.0005 + rN.Cai);

		Irel = krel*rN.u*rN.u*rN.v*rN.w*(rN.Carel - rN.Cai);
		Itr = (rN.Caup-rN.Carel)/tautr;
		Iup = Iupmax/(1.0+(kup/rN.Cai));
		Iupleak = rN.Caup*Iupmax/Caupmax;

		/*Add new Ikach, from Kneller et al., 2002*/
		//ach1 =  5.5/(1.0 + 9.13652/(pow(Ach,0.477811))); 
		// ach1 =  10.0/(1.0 + 9.13652/(pow(Ach,0.477811)));  Calculated in declaration area and includes gkach
		ach2 = (0.0517 + 0.4516/(1.0 + exp((rN.V+59.53)/17.18)));	// for [ACH]=1uM only!  Kneller 2002
 

		Ikach = 0.0;//ach1*ach2*(rN.V-Ek); // gkach is included in ach1

		Iion = Ina+Ik1+Ito+Ikur+Ikr+Iks+Ical+Ipca+Inak+Inaca+Ibna+Ibca+Ikach;
		Iion = Iion/Cm+Ist; // Return values in A/F units

		/* update variables */

		Fn = 1e-12*Vrel*Irel - (5e-13/F)*(0.5*Ical - 0.2*Inaca);

		if(rN.V==-47.13){
		am = 3.2;
		}else{
		am = 0.32*(rN.V+47.13)/(1-exp(-0.1*(rN.V+47.13)));
		}
		bm = 0.08*exp(-rN.V/11.0);
		taum = 1.0/(am+bm);
		infm = am*taum;

		if(rN.V>=-40.0){
		ah = 0.0;
		}else{
		ah = 0.135*exp( -1.0*(rN.V+80.0)/6.8);
		}
		if(rN.V>=-40.0){
		bh = 1.0/(0.13*(1+exp(-1.0*(rN.V+10.66)/11.1)));
		}else{
		bh = 3.56*exp(0.079*rN.V) + (3.1e5)*exp(0.35*rN.V);
		}
		tauh = 1.0/(ah+bh);
		infh = ah*tauh;

		if(rN.V>=-40.0){
		aj = 0.0;
		}else{
		aj = ( -127140.0*exp(0.2444*rN.V) - (3.474e-5)*exp(-0.04391*rN.V) )*(rN.V+37.78)/
				(1.0+exp(0.311*(rN.V+79.23)));
		}
		if(rN.V>=-40.0){
		bj = 0.3*exp((-2.535e-7 )*rN.V)/(1.0+exp(-0.1*(rN.V+32.0)));
		}else{
		bj = 0.1212*exp(-0.01052*rN.V)/(1.0+exp(-0.1378*(rN.V+40.14)));
		}
		tauj = 1.0/(aj+bj);
		infj = aj*tauj;

		asa = 0.65/( (exp(-1.0*(rN.V+10.0)/8.5)) + (exp(-1.0*(rN.V-30.0)/59.0)) );
		bsa = 0.65/(2.5 + exp((rN.V+82.0)/17.0) );
		infsa = 1.0/(1.0 + exp(-1.0*(rN.V+20.47)/17.54));
		tausa = (1.0/(asa+bsa))/Kq10;

		asi = 1.0/(18.53 + exp((rN.V+113.7)/10.95));
		bsi = 1.0/(35.56 + exp(-1.0*(rN.V+1.26)/7.44));
		infsi = 1.0/(1.0 + exp((rN.V+43.1)/5.3));
		tausi = (1.0/(asi+bsi))/Kq10;

		gkur = gkur_sub*(0.005 + 0.05/(1.0+exp(-1.0*(rN.V-15.0)/13.0)));
		aua = 0.65/(exp(-1.0*(rN.V+10.0)/8.5) + exp(-1.0*(rN.V-30.0)/59.0));
		bua = 0.65/(2.5+exp((rN.V+82.0)/17.0));
		infua = 1.0/(1.0 + exp(-1.0*(rN.V+30.3)/9.6)); 
		tauua = (1.0/(aua+bua))/Kq10;

		aui = 1.0/(21.0 + exp(-1.0*(rN.V-185.0)/28.0));
		bui = exp((rN.V-158.0)/16.0);
		infui = 1.0/(1.0 + exp((rN.V-99.45)/27.48)); 
		tauui = (1.0/(aui+bui))/Kq10;

		axr = 0.0003*(rN.V+14.1)/(1.0 - exp(-1.0*(rN.V+14.1)/5.0));
		bxr = (7.3898e-5)*(rN.V-3.3328)/(exp((rN.V-3.3328)/5.1237)-1.0);
		infxr = 1.0/(1.0 + exp(-1.0*(rN.V+14.1)/6.5)); 
		/*tauxr = (1.0/(axr+bxr))/Kq10;*/
		tauxr = (1.0/(axr+bxr));

		axs = (4e-5)*(rN.V-19.9)/(1.0 - exp(-1.0*(rN.V-19.9)/17.0));
		bxs = (3.5e-5)*(rN.V-19.9)/(exp((rN.V-19.9)/9.0) - 1.0);
		tauxs = 0.5/(axs+bxs);
		infxs = 1.0/sqrt(1.0 + exp(-1.0*(rN.V-19.9)/12.7));

		taud = (1.0 - exp(-1.0*(rN.V+10.0)/6.24))/(0.035*(rN.V+10.0)*(1.0+exp(-1.0*(rN.V+10.0)/6.24)));
		infd = 1.0/(1.0 + exp(-1.0*(rN.V+10.0)/8.0));

		tauf = 9.0/(0.0197*exp(-1.0*0.0337*0.0337*(rN.V+10.0)*(rN.V+10.0)) + 0.02);
		inff = 1.0/(1.0 + exp((rN.V+28.0)/6.9));

		inffca = 1.0/(1.0+rN.Cai/0.00035);

		sigma = (1.0/7.0)*(exp(Nao/67.3) - 1.0);
		fnak = 1.0/(1.0 + 0.1245*exp(-0.1*rN.V/(RTF)) + 0.0365*sigma*exp(-1.0*rN.V/(RTF)));

		infu = 1.0/(1.0+ exp(-1.0*(Fn - 3.4175e-13)/13.67e-16));

		tauv = 1.91 + 2.09/(1.0 + exp(-1.0*(Fn - 3.4175e-13)/13.67e-16));
		infv = 1.0 - 1.0/(1.0 + exp(-1.0*(Fn - 6.835e-14)/13.67e-16));

		tauw = 6.0*(1.0 - exp(-1.0*(rN.V-7.9)/5.0))/( (1.0 + 0.3*exp(-1.0*(rN.V-7.9)/5.0))*(rN.V-7.9) );
		infw = 1.0 - 1.0/(1.0 + exp(-1.0*(rN.V-40.0)/17.0));


		/* update differential equations */
		outN.Nai += dt*( (-3.0*Inak - 3.0*Inaca - Ibna - Ina)/(F*Vi) );
		outN.Ki += dt*( (2.0*Inak-Ik1-Ito-Ikur-Ikr-Iks-Ist)/(F*Vi) );
		B1 = (2.0*Inaca - Ipca - Ical - Ibca)/(2.0*F*Vi) + (Vup*(Iupleak-Iup) + Irel*Vrel)/Vi;
		B2 = 1.0 + Trpnmax*KmTrpn/((rN.Cai +KmTrpn)*(rN.Cai +KmTrpn)) + Cmdnmax*KmCmdn/((rN.Cai+KmCmdn)*(rN.Cai+KmCmdn));
		outN.Cai += dt*(B1/B2);
		outN.Caup += dt*(Iup-Iupleak - Itr*Vrel/Vup);
		outN.Carel += dt*( (Itr-Irel)*(1.0/(1.0 + (Csqnmax*KmCsqn)/((rN.Carel+KmCsqn)*(rN.Carel+KmCsqn)))));

		outN.m += dt*(infm - rN.m)/taum;
		outN.h += dt*(infh - rN.h)/tauh;
		outN.j += dt*(infj - rN.j)/tauj;
		outN.sa += dt*(infsa - rN.sa)/tausa;
		outN.si += dt*(infsi - rN.si)/tausi;
		outN.ua += dt*(infua - rN.ua)/tauua;
		outN.ui += dt*(infui - rN.ui)/tauui;
		outN.xr += dt*(infxr - rN.xr)/tauxr;
		outN.xs += dt*(infxs - rN.xs)/tauxs;
		outN.d += dt*(infd - rN.d)/taud;
		outN.f += dt*(inff - rN.f)/tauf;
		outN.fca += dt*(inffca - rN.fca)/taufca;
		outN.u += dt*(infu - rN.u)/tauu;
		outN.v += dt*(infv - rN.v)/tauv;
		outN.w += dt*(infw - rN.w)/tauw;


		/* rN.V += dt*-1.0*Iion/Cm; UPDATED IN THE MAIN PROGRAM TO INCLUDE COUPLING CURRENTS */

		return Iion;

		} // of Total_transmembrane_currents

};

