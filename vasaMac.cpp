
/*
  (dec07)
  Template fuer main() Programm; Ort:
  
  ~/versionedProjects/lib/templates/templateMain.cpp
  
  makefile dazu:
  
  ~/versionedProjects/lib/templates/makefile
  
  brauche nur das templaeMain.cpp und makefile; alles andere wird
  beim Compilieren dazugeholt!

  Achtung: Bei Klassendeklaration .h File mit Strichpunkt enden!!
  sonst moeglicherweise Fehlermeldung ohne jede Aussage!

  Achtung! Auch ohne .h File muss man bei $OBJECTS immer auch das 
  File mit der Main-Methode dazunehmen!
  (sonst "ld: undefined reference to main"

  vectors:

  vector<double>xtab;
  vector<double>xtabWithInit(ninit, 0);
  vector<double> xtabValInit={val1,val2,val3};
  vector<vehicle> vehicles;
  vehicles.push_back(oneVehicle); //out of bounds when vehicles[n]=oneVehicle
  vehicles[i];
  vehicles.size(); // size() elements

  2 ways of Calling an API function expecting arrays such as intp:

  intpextp (&xValsVec[0], &yValsVec[0], xValsVec.size(), x);
  intpextp (xValsVec.data(), yValsVec.data(), xValsVec.size(), x);

  Why? The elements of a vector are guaranteed to be contiguous.

*/


// c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//alternatively there is <cstdio> which declares
//everything in namespace std
//but the explicit "using namespace std;" puts
//everything in global namespace


// c++ 
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

// own
#include "general.h"

#include "Statistics.h"
#include "RandomUtils.h" // contains, e.g.,  myRand()
#include "InOut.h"


// constants

static const int NDATA_MAX=10000;// max. number of data points
static const int MAXSTR=500;// max. string length

// ################################################################
// central control variables and parameters
// ################################################################
class SimulationInput{
public:
  
  SimulationInput(){;}

  // spacetimes
  
  double tmax;
  double dt;
  double xmax;
  double dx;

  // parameters from SimulationInput

  double rhomax;
  double w=-1.4;        // wave speed for all groups and all gradients
  double mu=0.05;       // friction coefficient (!!wind drag not considered)
  double sigmarel=0.15; // relative in-class dispersion in free flow 

  // general control from SimulationInput

  int dntout;
  int dnxout;


  SimulationInput(const char projectName[]){
    char fname[MAXSTR];
    InOut inout;
    double data[NDATA_MAX];
    int nData;
    sprintf(fname,"%s.proj",projectName);
    inout.get_col(fname,1, nData, data);

    tmax=data[0];
    dt=data[1];
    xmax=data[2];
    dx=data[3];

    rhomax=data[4];
    w=data[5];
    mu=data[6];
    sigmarel=data[7];

    dntout=int(data[8]);
    dnxout=int(data[9]);
    
  }

private:
  double localVar;

};   // end HelperClass: don't forget semikolon!

SimulationInput sim;


// ############################################
// (1) Skier and geometry input
// ############################################

int nClasses;
double kTab[]={0,1,2,3,4,5,6,7,8,9,10};



// #participants in the groups

vector<int>n;
//int n[]={324,470,853,1154,1351,1751,1665,1897,2256,2136,2899};

// median pace (inverse speed, s/m) on level terrain
// (median because then v_k=1/invv_k)

vector<double>invv;
//double invv[]={0.15, 0.18, 0.21, 0.24, 0.27,0.30, 0.33, 0.36, 0.39, 0.42, 0.42};

// Median arrival times at starting line

vector<double>tau;
//double tau[]={0, 14.5, 33.1, 60.4, 84.5,120.6, 166.7, 211.0, 276, 358.5, 400};

// staggered wavestart (if no, all values are zero)
vector<double>tau_wave;



// from a GPX trajectory (file xxx) (jumps 1.5 km->1.7 km, 3.2 km->3.4 km
// 4.9 km -> 5.1 km 106m

vector<double>xTab; // defined in main from file
vector<double>nlaneTab;
vector<double>elevTab;

double nLanes( double x){
  return intpextp (xTab.data(), nlaneTab.data(), xTab.size(), x);
}


double elev(double x){
  return intpextp (xTab.data(), elevTab.data(), xTab.size(), x);
}


double gradient(double x){
  double dx=20;  // "smoothing width"
  return (elev(x+0.5*dx)-elev(x-0.5*dx))/dx;
}




// ############################################
// (3) Free flow (assumed dispersive, wave speed w globally const)
// ############################################



// free-flow pace 1/v as a function of the gradient

double factGradient(double gradient){return 1/(1+max(gradient,0.)/sim.mu);}
double invvgrad(int k, double gradient){
  return invv[k]/factGradient(gradient);
}

// warped time
// assumes relative slowdown due to gradient const for all classes

vector<double> initializeTimewarpfactorTab(){
  double timeAct=0;  
  double timeHoriz=0;  // time passed at given x for given v0=1
  vector<double> timewarpfactorTab;
  timewarpfactorTab.push_back(1.);

  for(int ix=0; ix<(int)(sim.xmax/sim.dx); ix++){
    cout<<"ix="<<ix<<endl;
    double x=(ix+0.5)*sim.dx;
    timeHoriz+=sim.dx; // v0=1
    timeAct += sim.dx/factGradient(gradient(x));
    timewarpfactorTab.push_back(timeAct/timeHoriz);
  }
  return timewarpfactorTab;
}

vector<double> timewarpfactorTab; // !! must be initialized in main()

double twarpFact(double x){
  return intpextp (timewarpfactorTab.data(), timewarpfactorTab.size(),
		   x, 0, sim.xmax);
}


// free-flow as a function of x and t (all rest of (2) old)

// universal relative pace distribution function
// (pace invv=1/v is Gaussian, not speed!)

double murel=1;
double ftilde(double Trel){
  return 1./sqrt(2*PI*sim.sigmarel*sim.sigmarel)
    *exp(-pow(Trel-murel,2)/(2*sim.sigmarel*sim.sigmarel));
}


double Qtotfree(int k, double x, double t){
  double teff=tau[k]+x*invv[k]*twarpFact(x);
  return (t>tau_wave[k])
    ? n[k]/teff * ftilde((t-tau_wave[k])/teff)
    :0;
}
  
double Qtotfree(double x, double t){
  double sum=0;
  for(int  k=0; k<nClasses; k++){
    sum+=Qtotfree(k,x,t);
  }
  return sum;
}

double rhototfree(int k,double x,double t){
  return Qtotfree(k,x,t)
    *invv[k]/factGradient(gradient(x));
}
  
double rhototfree(double x, double t){
  double sum=0;
  for(int  k=0; k<nClasses; k++){
    sum+=rhototfree(k,x,t);
  }
  return sum;
}

double Qfree(double x, double t){
  return Qtotfree(x,t)/nLanes(x);}
double rhofree(double x, double t){
  return rhototfree(x,t)/nLanes(x);}


// ############################################
// (4) congested flow: Single-class LWR
// with class-depenend bottleneck strength
// ############################################


// outflow from jam per lane

double Qmax(double V0){return -sim.w*sim.rhomax/(1-sim.w/V0);}
double QmaxPace(double invvloc){return -sim.w*sim.rhomax/(1-sim.w*invvloc);}

double capacClass(int k, double x){
  return nLanes(x)*QmaxPace(invvgrad(k,gradient(x)));
}

//!!! Quick hack to estimate the approx class median at (x,t)
// (later on by providing and updating vector kmedTab[i] and kmed(x)
// for the actual simulation time t)

double capacHack(double x, double t){
  double tJamFree=0.4*t; // w/o jams, the class arriving at time t would be
  // here at time 0.4*t
  double vHorizEff=x/tJamFree * twarpFact(x);
  int k=0;
  bool success=(vHorizEff>1/invv[0]);
  for(k=1; ((!success)&&(k<nClasses)); k++){
    success=(vHorizEff>1/invv[k]);
  }
  if(k>0){k--;}
  cout<<" (k="<<k<<") ";
  return nLanes(x)*QmaxPace(invvgrad(k,gradient(x)));
}

// median free speed on horizontal terrain
// as a function of (continuous) group index k

double V0levelfun(double k){
  return 1./intpextp(invv.data(), nClasses, k, 0, nClasses);
}

// estimated average group index as f(cumulated #participants)
// ncum<=n[0]/2: group 0; ncum>=N-n[nClasses-1]/2: group nClasses

double k_groupfun(double ncum){
  vector<double>ncumTab;
  ncumTab.push_back(0.5*n[0]);
  for(int k=1; k<nClasses; k++){
    ncumTab.push_back(ncumTab[k-1]+0.5*(n[k-1]+n[k]));
  }
  return intpextp(ncumTab.data(), kTab, nClasses, ncum);
}



// trajectory settings: get trajectories for all athletes in index_traj
// initialize xtraj with -1: "this athlete not yet crossed start"
vector<int>index_traj={100,500,1500,2500,4000,6000,8000,10000,13000,16000};
vector<double>xtraj(int(index_traj.size()),-1); 








// ###########################################
// if we need Gaussians
// ###########################################


double norm(double x){return 0.5*(1+erf(x/sqrt(2.)));}

//!!! define invnorm(p) more efficiently if needed!

double invnorm(double p){
  if(p<0.5){return -invnorm(1-p);}
  double n=100;
  double xmax=5;
  bool success=false;
  int i;
  for(i=0; (i<=n+1e-6)&&(!success); i++){
    //cout<<"i="<<i<<" p="<<p<<" norm(i/n*xmax)="<<norm(i/n*xmax)<<endl;
    success=(p<norm(i/n*xmax));
  }
  i--;
  double dp=norm(i/n*xmax)-norm((i-1)/n*xmax);
  double w=(p-norm((i-1)/n*xmax))/dp;
  cout<<"i="<<i<<" w="<<w<<endl;
  return (i==n) ? xmax : xmax*(i-1+w)/n;
}

// whatch out! reference ofstream& !!

void write_timestep(int it, ofstream& outfile,
		    vector<double> rhotot,
		    vector<double> Qtot, vector<double> capac,
		    vector<double> V0, vector<double> ncum,
		    vector<double> k_group,
		    vector<bool> flag_cong, vector<bool> flag_disp){
  if(it/sim.dntout<=1){
    outfile<<"#t\tx\trhotot\tQtot\tcapac\tV0\tncum\tk_group\tcong(0=false)\tdisp\n";
  }
  double t=it*sim.dt;
  // for i=int(rhotot.size()) some quantities are undef/zero=> decrement by 1
  for(int i=0; i<int(rhotot.size())-1; i+=sim.dnxout){
    double x=(i+0.5)*sim.dx;
    outfile<< std::fixed<<setprecision(2)<<t<<"\t"<<x<<"\t"
	   <<rhotot[i]<<"\t"
	   <<Qtot[i]<<"\t"
 	   <<capac[i]<<"\t"
 	   <<V0[i]<<"\t"<<setprecision(1)
 	   <<ncum[i]<<"\t"
 	   <<k_group[i]<<"\t"
	   <<setprecision(0)<<flag_cong[i]<<"\t\t"
	   <<flag_disp[i]<<endl;
  }
  outfile<<endl;
}
   
// whatch out! reference ofstream& !!

void write_trajectories(int it, ofstream& outfile,
			vector<int> index_traj, vector<double> xtraj){

  if(it<=1){
    outfile<<"#t\t";
    for(int i=0; i<int(index_traj.size()); i++){
      outfile<<"n="<<index_traj[i]<<"\t";
    }
    outfile<<endl;
  }
  outfile<<fixed<<setprecision(2);
  outfile<<it*sim.dt<<"\t";
  for(int i=0; i<int(xtraj.size()); i++){
    outfile<<xtraj[i]<<"\t";
  }
  outfile<<endl;
}

// whatch out! reference ofstream& !!

void write_SDD(int it, int ix, ofstream& outfile,
		    vector<double> Qtot,
		    vector<double> rhotot){
  if(it/sim.dntout<=1){
    outfile<<"#t[s]\tQtot[skiers/h]\trhotot[skiers/m]\n";
  }
  outfile<< std::fixed<<setprecision(0)<<it*sim.dt<<"\t"
	 <<setprecision(0)<<3600*Qtot[ix]<<"\t\t"
	 <<setprecision(2)<<rhotot[ix]<<endl;
}

// prepare HoegstaPunkten_splitTimes.csv to calibrate with outfile_SDD data

void write_calibHoechstaPunkten(string str_infile, string str_calibfile){

  InOut inout;
  int nSplitData=inout.getNumberOfDataLines(str_infile.c_str());
  
  int nt=ceil(sim.tmax/sim.dt);
  int nIntervals=nt/sim.dntout;
  
  int classSplitData[nSplitData];
  double timeSplitData[nSplitData];
  
  inout.get_col(str_infile.c_str(),2, nSplitData, classSplitData);
  inout.get_col(str_infile.c_str(),3, nSplitData, timeSplitData);

  ofstream outfile(str_calibfile, ios::out);
  outfile<<"#t[s]\tnPassed\tflow[1/h]"<<endl;

  vector<double>flowSplitData;
  int iSplit=0;
  for(int i=1; i<nIntervals; i++){
    int nPassed=0;
    while( (timeSplitData[iSplit]<(i+1)*sim.dt*sim.dntout)
	   &&(iSplit<nSplitData)){
      iSplit++;
      nPassed++;
    }
    flowSplitData.push_back(nPassed/sim.dntout);
    outfile<< ((i+0)*sim.dt*sim.dntout) <<"\t"
	   <<nPassed <<"\t"
	   <<(3600*nPassed/(sim.dt*sim.dntout))<<endl;
    if(false){
      cout<<"HoegstaPunkten: t="<<(i+0)*sim.dt*sim.dntout
	  <<" iSplit="<<iSplit
	  <<" timeSplitData[iSplit]="<<timeSplitData[iSplit]
	  <<" nPassed="<<nPassed
	  <<" flow[invh]="<<3600*nPassed/sim.dntout<<endl;
    }
  }
  cout<<"counted "<<iSplit<<" skiers at Hoechsta Punkten"<<endl;
}

  

//#####################################################
//#####################################################

int main(int argc, char* argv[]) {

  
  
  //###############################
  // Input
  // ##############################

  char   projName[MAXSTR];

  if (argc!=2){ // argc=number of cmdline params + 1
    cerr <<"\nCalling syntax: vasaMac projName w/o extension\n";
    cerr <<"Example: vasaMac sim1";
    exit(-1);
  }
  sprintf(projName,"%s",argv[1]);
  sim=SimulationInput(projName);
  if(true){
    cout<<"sim.xmax="<<sim.xmax<<endl
	<<"sim.dx="<<sim.dx<<endl
	<<"sim.tmax="<<sim.tmax<<endl
	<<"sim.dt="<<sim.dt<<endl;
  }

  int nx=ceil(sim.xmax/sim.dx);  // cell i in [i*dx, (i+1)*dx], i=0..nx-1
  int nt=ceil(sim.tmax/sim.dt); 
  int nxSDD=min(nx-1,int(2980/sim.dx));
  
  char   outName_macro[MAXSTR+15];
  char   outName_traj[MAXSTR+15];
  char   outName_SDD[MAXSTR+15];
  char   inName_trackParams[MAXSTR+15];
  char   inName_skiers[MAXSTR+15];
  string str_indataHoechstaPunkten="HoegstaPunkten_splitTimes.csv";
  string str_calibHoechstaPunkten=string(projName)+".HoegstaPunktenCalib";

  
  sprintf(outName_macro,"%s.%s",projName,"macro");
  sprintf(outName_traj,"%s.%s",projName,"traj");
  sprintf(outName_SDD,"%s.x%i",projName,int(sim.dx*nxSDD));
  sprintf(inName_trackParams,"%s.%s",projName,"trackParams");
  sprintf(inName_skiers,"%s.%s",projName,"skiers");

  ofstream outfile_macro(outName_macro, ios::out);
  ofstream outfile_traj(outName_traj, ios::out);
  ofstream outfile_SDD(outName_SDD, ios::out);
  if(!outfile_macro){
    cerr<<"vasaMac.main: Error: cannot open files for writing"<<endl;
    exit(-1);
  }

  outfile_macro<<"#vasaMac with parameters rhomax="<<sim.rhomax<<" w="<<sim.w
	       <<" mu="<<sim.mu<<" sigmarel="<<sim.sigmarel<<endl;
  outfile_traj<<"#vasaMac with parameters rhomax="<<sim.rhomax<<" w="<<sim.w
	       <<" mu="<<sim.mu<<" sigmarel="<<sim.sigmarel<<endl;
  outfile_SDD<<"#vasaMac with parameters rhomax="<<sim.rhomax<<" w="<<sim.w
	       <<" mu="<<sim.mu<<" sigmarel="<<sim.sigmarel<<endl
	     <<"#Virtual macroscopic detector data at x="
	     <<(nxSDD*sim.dx)<<endl;


      
// #########################################################
// track attributes
// #########################################################

  InOut inout;
  int nTrackData=inout.getNumberOfDataLines(inName_trackParams);
  double data1[nTrackData]; 
  double data2[nTrackData]; 
  double data3[nTrackData]; 
  inout.get_col(inName_trackParams,1, nTrackData, data1);
  inout.get_col(inName_trackParams,2, nTrackData, data2);
  inout.get_col(inName_trackParams,3, nTrackData, data3);


  // no easily understandable way to transform arrays to vectors
  // -> do it the simple way
  
  for(int i=0; i<nTrackData; i++){
    xTab.push_back(data1[i]);
    nlaneTab.push_back(data2[i]);		
    elevTab.push_back(data3[i]);
  }

  for(int i=0; i<nTrackData; i++){
    cout<<"i="<<i<<" x="<<xTab[i]<<" nlaneTab="<<nlaneTab[i]
	<<" elevTab="<<elevTab[i]<<endl;
  }


  // #########################################################
  // Skier attributes
  // #########################################################

  nClasses=inout.getNumberOfDataLines(inName_skiers);
  double skidata2[nClasses]; 
  double skidata3[nClasses]; 
  double skidata4[nClasses]; 
  double skidata5[nClasses];
  
  inout.get_col(inName_skiers,2, nClasses, skidata2);
  inout.get_col(inName_skiers,3, nClasses, skidata3);
  inout.get_col(inName_skiers,4, nClasses, skidata4);
  inout.get_col(inName_skiers,5, nClasses, skidata5);

  
  for(int i=0; i<nClasses; i++){
    n.push_back(skidata2[i]);		
    invv.push_back(skidata3[i]);
    tau.push_back(skidata4[i]);
    tau_wave.push_back(skidata5[i]);
  }
  

  for(int k=0; k<nClasses; k++){
    cout<<"k="<<k<<" n[k]="<<n[k]
      //<<" invv[k]="<<invv[k]
	<<endl;
  }
  
  // ###############################################################


  
  // need to initialize this variable defined before main !!
  
  timewarpfactorTab=initializeTimewarpfactorTab();

  // #athletes

  int ntot=0;
  for(int k=0; k<nClasses; k++){ntot+=n[k];}

  
  // test of flows and timewarps
  
  int k=10;
  cout<<" k="<<k<<endl;

  for(int i=0; i<81; i++){
    double x=i*50;
    cout<<"x="<<x<<" nLanes="<<nLanes(x)<<" gradient="<<gradient(x)
	<<" factGradient="<<factGradient(gradient(x))
        <<" vMed(k=5)="<<1/invvgrad(k,gradient(x))
	<<" twarpFact="<<twarpFact(x)
        <<" capacClass(k=5)="<<capacClass(k,x)
	<<endl;
  }

  for(int i=9; i<25; i++){
    double x=i*100;
    cout<<endl<<"x="<<x<<" capacClass(3,x)="<<capacClass(3,x)<<endl;
    for(int j=5; j<60; j++){
      double t=j*60;
      cout<<"  x="<<x<<" t="<<t
	  <<" Qtotfree(x,t)="<<Qtotfree(x,t)
	  <<" capacHack(x,t)="<<capacHack(x,t)
	  <<endl;
    }
    double ntotsum=0;
    for(int j=0; j<100; j++){
      double t=(j+0.5)*60;
      ntotsum+=60*Qtotfree(x,t);
    }
    cout<<" Test athlete conservation: ntotsum="<<ntotsum
	<<" ntot="<<ntot<<endl;
    
  }

  /*
  // space test
  for(int i=0; i<12; i++){
    double x=(i+0.5)*100;
    for(int j=0; j<25; j++){
      double t=j*30;
      cout<<"  x="<<x<<" t="<<t
	  <<" Qtotfree(x,t)="<<Qtotfree(x,t)
	  <<" capacHack(x,t)="<<capacHack(x,t)
	  <<endl;
    }
  }
  */

  /* 
  // time test 
  for(int it=0; it<20; it++){
    double dt=60;
    double sim.dx=50;
    double t=it*dt;
    cout<<"\nt="<<t<<endl;
    for(int i=0; i<22; i++){
      double x=(i+0.5)*sim.dx;
      cout<<" x="<<x
	  <<" Qtotfree(x,t)="<<Qtotfree(x,t)
	  <<" rhototfree(x,t)="<<Qtotfree(x,t)*t/(x+60)
	  <<" capacHack(x,t)="<<capacHack(x,t)
	  <<endl;
    }
    double ncum=0;
    for(int i=0; i<100; i++){
      double x=(i+0.5)*sim.dx;
      double v=x/(t+0.1);
      ncum+=sim.dx*Qtotfree(x,t)/v;
    }
    cout<<"ncum(x=0..5000)="<<ncum<<" ntot="<<ntot<<endl;
    
  }
  */
  
  //test norm(x)=0.5*(1+erf(x)) and own invnorm(p)

  if(false){
    cout<<"Test: norm0)="<<norm(0)<<" norm(1.69)="<<norm(1.69)
        <<" norm(-1.96)="<<norm(-1.96)<<endl;
    cout<<"invnorm(norm(-0.2))="<<invnorm(norm(-0.2))<<endl;
    cout<<"invnorm(norm(0))="<<invnorm(norm(0))<<endl;
    cout<<"invnorm(norm(0.001))="<<invnorm(norm(0.001))<<endl;
    cout<<"invnorm(norm(3))="<<invnorm(norm(3))<<endl;
    cout<<"invnorm(norm(0.02))="<<invnorm(norm(0.02))<<endl;
    cout<<"invnorm(norm(-10))="<<invnorm(norm(-10))<<endl;
  }

  // test V0 and k_group

  

  


  //###############################
  // Simulation
  // ##############################
  



  // state variables for actual simulation (for a given time)

  vector<double>Qtot(nx,0);  // (initially) nx elemenets with value 0
  vector<double>Qtotold(nx,0);  // !! needed for simultaneous update
  vector<double>rhotot(nx,0);  
  vector<double>rhoctot(nx,0);  
  vector<double>capac(nx,0);  // local capacity (dep. on x, group, grade)
  vector<double>V0(nx,0);     // local max speed (dep. on x, group, grade)
  vector<double>V(nx,-1); //only need local speed for trajectories
  vector<bool>flag_cong(nx,false);
  vector<bool>flag_disp(nx,true);
  vector<double>ncum(nx,0); // cum #athletes  int_0^x rhotot dx'
  vector<double>k_group(nx,0); // avg starting group with smooth transit.

  

  //#########################################################
  // actual main simulation loop
  // from time step (it-1)*sim.dt (initial state vars)
  // to it*sim.dt (updated state vars)
  //#########################################################

  double t_onlyDisp=400.; // only dispersion for t<t_onlyDisp s 

  
 


  
  for(int it=1; it<nt; it++){
    double t=it*sim.dt;
    if(it%10==0){cout<<"it="<<it<<" t="<<t<<endl;}



    // update local capacities, history for simultaneous update
    // and flags for dispersive region and congestion
    
    for(int i=0; i<nx; i++){
      Qtotold[i]=Qtot[i];
      double x=(i+0.5)*sim.dx;
      k_group[i]=k_groupfun(ncum[i]);
      double V0level=V0levelfun(k_group[i]);
      V0[i]=V0level*factGradient(gradient(x));
      capac[i]=nLanes(x)*Qmax(V0[i]);
      rhoctot[i]=capac[i]/V0[i];

      flag_disp[i]=(t<=t_onlyDisp);
      if(t>t_onlyDisp){flag_cong[i]=(rhotot[i]>rhoctot[i]);}
 
   }

    // boundary conditions

    Qtot[0]=Qtotfree(0.5*sim.dx,t); // dispersive FT at upstream boundary
    rhotot[0]=Qtot[0]/V0[0];
    Qtot[nx-1]=Qtot[nx-2];  // Homogeneous Von Neumann downstream BC
    rhotot[nx-1]=rhotot[nx-2];


    // debug

    int itdebug=680; int imin=104; int imax=112;
    if(it==itdebug){
      cout <<"debug: begin of timestep it="<<it<<" t="<<t<<endl;
      for(int i=imin; i<=imax; i++){
	double x=(i+0.5)*sim.dx;
	cout <<"i="<<i<<" x="<<x
	     <<" rhotot[i]="<<rhotot[i]
	     <<" rhomaxtot="<<sim.rhomax*nLanes(x)
	     <<" Qtot[i]="<<Qtot[i]<<" capac[i]="<<capac[i]
	     <<endl;
      }
    }
    
    // main update for t>sim.dt*nt_onlyDisp

    bool sequentialUpdate=false; // true: use Qtot[i]; false: use Qtotold[i]
    
    for(int i=1; i<nx-1; i++){
      double x=(i+0.5)*sim.dx;
    
      // update as for it<nt_onlyDispt if in dispersion regime

      if(flag_disp[i]){
        Qtot[i]=Qtotfree(x,t);
        rhotot[i]=Qtot[i]/V0[i];
      }

      // update if in supply-demand regime

      else{

	// upstream cell boundary flow as min of supply and demand
	// flag_cong determined for all cells before -> parallel update

	double Qup=(sequentialUpdate) ? Qtot[i-1] : Qtotold[i-1];
	double Qcenter=(sequentialUpdate) ? Qtot[i] : Qtotold[i];
	double Qdown=(sequentialUpdate) ? Qtot[i+1] : Qtotold[i+1];
	
	double D_up=(flag_cong[i-1]) ? capac[i-1] : Qup;
	double S_up=(flag_cong[i]) ? Qcenter : capac[i];
	double Q_up=min(S_up, D_up);

	// downstream cell boundary flow as min of supply and demand
	
	double D_down=(flag_cong[i]) ? capac[i] : Qcenter;
	double S_down=(flag_cong[i+1]) ? Qdown : capac[i+1];
	double Q_down=min(S_down, D_down);

	//debug
	if((it==itdebug)&&(i==108)){
	  for(int is=107; is<=109; is++){
	    double rhomaxtot=sim.rhomax*nLanes((is+0.5)*sim.dx);
	    cerr<<"\nanalysis: before state update: is="<<is
		<<" flag_cong[is]="<<flag_cong[is]
		<<" capac[is]="<<capac[is]<<endl
		<<"    Qtotold[is]="<<Qtotold[is]
		<<" Qtot_update="<<((rhotot[is]<rhoctot[is]) 
				    ? rhotot[is]*V0[is] : -sim.w*(rhomaxtot-rhotot[is]))
		<<" Q_up="<<Q_up<<" Q_down="<<Q_down
		<<" rhototold[is]="<<rhotot[is]
		<<" rhotot_update="<<rhotot[is]+sim.dt/sim.dx*(Q_up-Q_down)
		<<" rhomaxtot="<<rhomaxtot
		<<endl;
	  }
	}
	
	// state update

	double rhomaxtot=sim.rhomax*nLanes(x);
	rhotot[i] +=sim.dt/sim.dx*(Q_up-Q_down);
	Qtot[i] = (rhotot[i]<rhoctot[i]) // triang FD
	  ? rhotot[i]*V0[i] : max(0.,-sim.w*(rhomaxtot-rhotot[i]));


	
	if((rhotot[i]>rhomaxtot)||(Qtot[i]>1.0001*capac[i])){
	  cerr<<"\n\nerror: i="<<i<<" x="<<x<<" it="<<it<<" t="<<t
	      <<" either rhotot[i]="<<rhotot[i]<<">rhomaxtot="<<rhomaxtot
	      <<" or Qtot[i]="<<Qtot[i]<<">capac[i]="<<capac[i]
	      <<endl;
	  exit(0);
	}
      } // end supply-demand update

      
      ncum[i]+=Qtot[i]*sim.dt; // update passed athletes
    }



   // ############################################################
    // generation of the artificial trajectories
    // ############################################################

    int ixstart=1;  // xstart=30;
    double xstart=(ixstart+0.5)*sim.dx;// need a distance where athlete number OK
    for(int i=0; i<nx; i++){
       V[i]=((rhotot[i]<1e-6) ? V0[i] : Qtot[i]/rhotot[i]);
    }
    

    for(int itraj=0; itraj<int(index_traj.size()); itraj++){
      if(xtraj[itraj]<-0.9999){ // skier not yet crossed
	if(ncum[ixstart]>=index_traj[itraj]){// let this skier start
	  cout<<"skier "<<index_traj[itraj]<<" started at x="<<xstart<<endl;
	  xtraj[itraj]=xstart;} 
      }
      
      if((xtraj[itraj]>ixstart-0.0001)&&(xtraj[itraj]<sim.xmax)){
	double Vloc=intp(V.data(), nx, xtraj[itraj], 0.5*sim.dx, (nx-0.5)*sim.dx);
	xtraj[itraj]+=Vloc*sim.dt;
      }
      else{//nothing
      }

      if(false){
	cout<<"t="<<t<<" xtraj[itraj]="<<xtraj[itraj]
	    <<" ncum[ixstart]="<<ncum[ixstart]
	    <<" index_traj[itraj]="<<index_traj[itraj]
	    <<endl;
      }
	
    }

    
    
    // test output

    if(false){
      //if(it%60==0){
      cout<<"\n time t="<<t<<endl;
      for(int i=0; i<nx; i++){
        double x=(i+0.5)*sim.dx;
	cout<<" x="<<x
	    <<" Qtot[i]="<<Qtot[i]
	    <<" rhotot[i]="<<rhotot[i]
	    <<" V0[i]="<<V0[i]
	    <<" flag_disp[i]="<<( (flag_disp[i]) ? "true" : "false")
	    <<" flag_cong[i]="<<( (flag_cong[i]) ? "true" : "false")
	    <<endl;
      }
    }
    
    // writing output (inside time loop)

    if(it%sim.dntout==0){
      write_timestep(it,outfile_macro,rhotot,Qtot,capac,
		     V0, ncum,k_group,flag_cong,flag_disp);
    }
    
    write_trajectories(it, outfile_traj, index_traj, xtraj);
    
    if(it%sim.dntout==0){
      write_SDD(it, nxSDD, outfile_SDD, Qtot, rhotot);
    }


  } // end sim

  cout<<"Simulated skiers at x="<<(nxSDD*sim.dx)
      <<" (no dropout): "<<ncum[nxSDD]<<endl;
  
  write_calibHoechstaPunkten(str_indataHoechstaPunkten,
			    str_calibHoechstaPunkten);

  cout<<"\nwrote"<<endl
      <<outName_macro<<endl
      <<outName_traj<<endl
      <<outName_SDD<<endl
      <<str_calibHoechstaPunkten
      <<endl;
  cout<<"\nvasaMac call finished\n\n";

 
  return(0);
}

// #########################################################


/*

  // Statistics stat;


  //###############################
  // Example random numbers (#include "RandomUtils.h" written by Arne
  // ##############################

  cout <<"\nrandom numbers: myRand() starts with fixed seed by default"<<endl;
  double rand= myRand()-0.5; // rand sim G(-1/2,1/2)
  cout <<"\nrandom number rand sim G(-1/2,1/2): rand="<<rand<<endl;
  rand= myRand()-0.5; // rand sim G(-1/2,1/2)
  cout <<"random number rand sim G(-1/2,1/2): rand="<<rand<<endl;

  // Random seed:
  // srand(seed); e.g., srand(42);
  // Arne's function setRandomSeed(); works only OUTSIDE of scripts;
  // otherwise, obviously the starting time of script used!

  //#####################################################
  */

