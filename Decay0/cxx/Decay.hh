// Translated from Fortran  16.03.2010. Aleksandra Bialek <abialek@ualberta.ca>
// Brief instruction on how to use DECAY as a standalone program:
// 1. Compile 
// 2. Run DECAY and answer the questions; in particular, select an option to
//    write generated events in file. After, you can read these events
//    and simulate them with GEANT, EGS, MCNP or other simulation package.
// Please refer to: O.A.Ponkratenko, V.I.Tretyak, Yu.G.Zdesenko,
// "Event Generator DECAY4 for Simulating Double-Beta 
// Processes and Decays of Radioactive Nuclei",
// Phys. At. Nucl. 63 (2000) 1282 (nucl-ex/0104018).
// This paper describes mainly DECAY4 generator, and DECAY0 was mentioned
// as initial version of the code.
// Please send questions, comments and proposals to: 
// tretyak@kinr.kiev.ua or tretyak@lngs.infn.it.
/////////////////////////////////////////////////////
//  Generation of events of decay of natural radioactive isotopes and
//  various modes of double beta decay.
//  DECAY units: energy and moment - MeV and MeV/c;
//               time              - sec.

#ifndef __DecayH__
#define __DecayH__

#include <TRandom1.h>

#include <complex>
#include "TF1.h"
const double emass=0.51099906;

using namespace std;

class Decay  {

  public:
   bool gop;
   char chfile[40], chn[16],chdspin[4], chnuclide[16];
   char chart[10][80],chmodebb[80];
   int icurrent, fNbPart, nevents, ievstart, irndmst, iwrfile, i2bbs;
   float  ebb1, ebb2, pmass, phi,phi1,phi2, toallevents,spmax;
   int allevents, modebb0 ;
   int npgeant[100];
   float fPmoment[3][100], fPtime[100], spthe1[4300],spthe2[4300];
   Double_t elow, ehigh;
   float El, e0, e1, e2,m, t, levelE, p1,p2,p4;
   int fLevelE;
   int nartparts;
   int fArtZd[10];
   double fSl[50];
   int fSlSize;
   float artemin[10], artemax[10], arttmin[10], arttmax[10],
          artfmin[10], artfmax[10], artQb[10], arttdelay[10],
          artthalf[10];
    int nartnpg[10];
    int fi2bbs,filevel,ilevel,fModebb,istart,fistart,fier,ier,itrans02;
    float fQbb,fEK,Zd,fZdbb;
    double fQbeta; //beta energy endpoint (MeV; Qbeta>50 eV);
    int fZdtr;     //atomic number of daughter nucleus (Zdtr>0 for e- and Zdtr<0 for e+ particles
    float fEbinde,fEbindeK,fEbindeL,fEbindeM,conve,convp,tevst,tdnuc1;
    int fStartbb;
    float fTdlev;  //time of decay of level (sec);
    float fTclev;  //time of creation of level from which particle will be emitted (sec);
    float fThlev;  //level halflife (sec).
    double fThnuc; //nucleus halflife (sec);
    float fTdnuc;  //time of decay of nucleus (sec);
    float fEgamma; //gamma-ray energy (MeV) [=difference in state energies];
    float fTevst;  //time of event's start (sec);
    float fC1,fC2,fC3,fC4, fKf; 

  public:
   Decay(){};
   ~Decay(){}; 
   void GENBBsub();
   void GENBBdia();
   float GetRandom();
   void particle(int np, float E1, float E2, float teta1, float teta2,
                 float phi1, float phi2, float tclev, float thlev);
   void pair(float Epair);
   void bb();
   float fe1_mod(float &e1);
   float fe2_mod(float &e2);
   float fermi(const float &Z,const float &E);
   void beta(float Qbeta,float tcnuc,float thnuc);
   void beta1f(float Qbeta,float tcnuc,float thnuc,float c1,float c2,float c3,float c4);
   void beta1fu(float Qbeta,float tcnuc, float thnuc,float c1,float c2,float c3,float c4);
   void beta2f(float Qbeta,float tcnuc,float thnuc,int kf, float c1,float c2,float c3,float c4);
   void tgold(float a,float b,TF1 &f,float eps, int nmin,float&xextt,float &fextr);
   /************************************************/
    Double_t funbeta(Double_t *x, Double_t *par);
    Double_t funbeta1fu(Double_t *x, Double_t *par);
    Double_t funbeta1f(Double_t *x, Double_t *par);
    Double_t funbeta2f(Double_t *x, Double_t *par);
    /************************************************/ 
   
   void Ti48low();
   void Fe58low();
   void Se76low();
   void Ge74low();
   void Kr82low();
   void Mo96low();
   void Zr92low();
   void Ru100low();
   void Pd106low();
   void Sn116low();
   void Cd112low();
   void Te124low();
   void Xe130low();
   void Ba136low();
   void Sm148low();
   void Sm150low();
   void As79(float tcnuc);
   void At214(float tcnuc);
   void Ac228(float tcnuc);
   void Bi207(float tcnuc);
   void Bi210(float tcnuc);
   void Bi212(float tcnuc);
   void Bi214(float tcnuc);
   void Co60(float tcnuc);
   void Cs136(float tcnuc);
   void Eu147(float tcnuc);
   void Eu152(float tcnuc);
   void Eu154(float tcnuc);
   void Gd146(float tcnuc);
   void I126(float tcnuc);
   void I133(float tcnuc);
   void I134(float tcnuc);
   void I135(float tcnuc);
   void K40(float tcnuc);
   void K42(float tcnuc);
   void Pa234(float tcnuc);
   void Pb211(float tcnuc);
   void Pb212(float tcnuc);
   void Pb214(float tcnuc);
   void Po214(float tcnuc);
   void Rn218(float tcnuc);
   void Ra222(float tcnuc);
   void Ra228(float tcnuc);
   void Rh106(float tcnuc);
   void Sb125(float tcnuc);
   void Sb126(float tcnuc);
   void Sb133(float tcnuc);
   void Sc48(float tcnuc);
   void Ta182(float tcnuc);
   void Te133(float tcnuc);
   void Te133m(float tcnuc);
   void Te134(float tcnuc);
   void Th234(float tcnuc);
   void Tl207(float tcnuc);
   void Tl208(float tcnuc);
   void Xe133(float tcnuc);
   void Xe135(float tcnuc);
   void Y88(float tcnuc);
   void Zn65(float tcnuc);
   void Nb96(float tcnuc);
   void PbAtShell(int KLMenergy);
   void nucltransK(float Egamma,float Ebinde,float conve,float convp);
   void nucltransKL(float Egamma,float EbindeK,float conveK,float EbindeL,float conveL,float convp);
   void nucltransKLM(float Egamma,float EbindeK,float conveK,float EbindeL,float conveL,float EbindeM,float conveM, float convp);
   void nucltransKLM_Pb(float Egamma,float EbindeK,float conveK,float EbindeL,float conveL,float EbindeM,float conveM, float convp);

   complex<double> cgamma(complex<double> z);
   double divdif (double xtab[], double xval);
   FILE*file;
   TRandom1 *rnd ;   

   double operator() (double *x, double *par)
   {
     //to use with GaussLegendreIntegrator
     // function implementation using class data members
      Double_t fe1mod=0.;
      float xx ;//e2
     
     if (fModebb!=10){
        xx = x[0];//e2
        //   float yy = x[1];//e1
        par[0]=e1;
        par[1]=e0;
        par[2]=emass;
        par[3]=fZdbb;
      }

      if (fModebb==4 && xx<(par[1]-par[0])) 
         fe1mod=(par[0]+par[2])*sqrt(par[0]*(par[0]+2.*par[2]))*fermi(par[3],par[0])*
                (xx+par[2])*sqrt(xx*(xx+2.*par[2]))*fermi(par[3],xx)*
                pow(par[1]-par[0]-xx,5);
      
      if (fModebb==5 && xx<(par[1]-par[0])) 
         fe1mod=(par[0]+par[2])*sqrt(par[0]*(par[0]+2.*par[2]))*fermi(par[3],par[0])*
                (xx+par[2])*sqrt(xx*(xx+2.*par[2]))*fermi(par[3],xx)*
                (par[1]-par[0]-xx);
      
      if (fModebb==6 && xx<(par[1]-par[0])) 
         fe1mod=(par[0]+par[2])*sqrt(par[0]*(par[0]+2.*par[2]))*fermi(par[3],par[0])*
                (xx+par[2])*sqrt(xx*(xx+2.*par[2]))*fermi(par[3],xx)*
                pow(par[1]-par[0]-xx,3);
      
      if (fModebb==8 && xx<(par[1]-par[0])) 
         fe1mod=(par[0]+par[2])*sqrt(par[0]*(par[0]+2.*par[2]))*fermi(par[3],par[0])*
                (xx+par[2])*sqrt(xx*(xx+2.*par[2]))*fermi(par[3],xx)*
                pow(par[1]-par[0]-xx,7)* pow(par[0]-xx,2);
      
      if (fModebb==13 && xx<(par[1]-par[0])) 
         fe1mod=(par[0]+par[2])*sqrt(par[0]*(par[0]+2.*par[2]))*fermi(par[3],par[0])*
                (xx+par[2])*sqrt(xx*(xx+2.*par[2]))*fermi(par[3],xx)*
                pow(par[1]-par[0]-xx,7);
      
      if (fModebb==14 && xx<(par[1]-par[0])) 
         fe1mod=(par[0]+par[2])*sqrt(par[0]*(par[0]+2.*par[2]))*fermi(par[3],par[0])*
                (xx+par[2])*sqrt(xx*(xx+2.*par[2]))*fermi(par[3],xx)*
                pow(par[1]-par[0]-xx,2);
      
      if (fModebb==15 && xx<(par[1]-par[0])) 
         fe1mod=(par[0]+par[2])*sqrt(par[0]*(par[0]+2.*par[2]))*fermi(par[3],par[0])*
                (xx+par[2])*sqrt(xx*(xx+2.*par[2]))*fermi(par[3],xx)*
                pow(par[1]-par[0]-xx,5)*(9*pow(par[1]-par[0]-xx,2)+21*pow(xx-par[0],2));
      
      if (fModebb==16 && xx<(par[1]-par[0])) 
         fe1mod=(par[0]+par[2])*sqrt(par[0]*(par[0]+2.*par[2]))*fermi(par[3],par[0])*
                (xx+par[2])*sqrt(xx*(xx+2.*par[2]))*fermi(par[3],xx)*
                pow(par[1]-par[0]-xx,5)*pow(xx-par[0],2);
      
      return fe1mod;
   }
};

#endif
