#include <iostream>
#include <cstring>
#include <cmath>
#include "Decay.hh"
#include <stdio.h>
#include <algorithm>
#include <complex>
#include "TF2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AllIntegrationTypes.h"

#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"

#include "Math/IFunction.h"
#include "Math/WrappedParamFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/IFunctionfwd.h"
#include "Math/GaussLegendreIntegrator.h"
#include <unistd.h>
#define pi 3.1415927
#define twopi 2*pi

/***************
 Particle masses (MeV) - in accordance with GEANT 3.21 manual of October, 1994:
 1 - gamma         2 - positron     3 - electron
 4 - neutrino      5 - muon+        6 - muon-
 7 - pion0         8 - pion+        9 - pion-
10 - kaon0 long   11 - kaon+       12 - kaon-
13 - neutron      14 - proton      15 - antiproton
16 - kaon0 short  17 - eta         18 - lambda
19 - sigma+       20 - sigma0      21 - sigma-
22 - xi0          23 - xi-         24 - omega
25 - antineutron  26 - antilambda  27 - antisigma-
28 - antisigma0   29 - antisigma+  30 - antixi0
31 - antixi+      32 - antiomega+  45 - deuteron
46 - tritium      47 - alpha       48 - geantino
49 - He3          50 - Cerenkov
0-dummy
**************/

float datamass[51]={0,    0., 0.51099906,    0.51099906,        0.,     105.658389, 
            105.658389,   134.9764,      139.5700,  139.5700,        497.672, 
               493.677,    493.677,     939.56563, 938.27231,      938.27231, 
               497.672,     547.45,      1115.684,   1189.37,        1192.55, 
              1197.436,     1314.9,       1321.32,   1672.45,      939.56563, 
              1115.684,    1189.37,       1192.55,  1197.436,         1314.9, 
               1321.32,    1672.45,            0.,        0.,             0.,       
                    0.,         0.,            0.,        0.,             0.,       
                    0.,         0.,            0.,        0.,       1875.613,    
               2809.25,   3727.417,            0.,   2809.23,              0.};

               
// Values of natural logarithm of the standard values of momentum  (in units m_e*c) from: 
// H.Behrens, J.Janecke, "Numerical tables for beta-decay and electron capture", Berlin, Springer-Verlag, 1969.
// Range of momenta correspond to kinetic energy range from 2.55 keV to 25.0 MeV.
double plog69[]={-2.302585,  -1.609438,  -1.203973,  -0.9162907, -0.6931472, -0.5108256, -0.3566750, -0.2231435, -0.1053605,  0.0000000, 0.1823216,  0.3364722,  0.4700036,  0.5877866,  0.6931472, 0.7884574,  0.8754688,  0.9555114,  1.029619,   1.098612, 1.163151,   1.223776,   1.280934,   1.335001,   1.386294, 1.504077,   1.609438,   1.704748,   1.791759,   1.871802, 1.945910,   2.014903,   2.079442,   2.197225,   2.302585, 2.397895,   2.484907,   2.564949,   2.639057,   2.772589, 2.890372,   2.995732,   3.218876,   3.401197,  3.555348, 3.688879,   3.806663,   3.912023 };               
               
Double_t Decay::funbeta(Double_t *x, Double_t *par)
{
  Double_t xx =x[0];
  Double_t funbeta=0;
  if(xx>0.) funbeta=sqrt(xx*(xx+2.*par[1]))*(xx+par[1])*pow((par[0]-xx),2)*fermi(par[2],xx);
  return funbeta;
}

Double_t Decay::funbeta1f(Double_t *x, Double_t *par)
{
  Double_t xx =x[0];
  Double_t funbeta1f=0.; 
  if (xx>0){
     float all=sqrt(xx*(xx+2.*par[1]))*(xx + par[1])*pow(par[0]-xx,2)*fermi(par[2],xx);
     float w=xx/par[1]+1.;
     float cf=1.+ fC1/w+ fC2*w+ fC3*w*w+ fC4*w*w*w;
     funbeta1f=all*cf;
  }
     return funbeta1f;
}


Double_t Decay::funbeta1fu(Double_t *x, Double_t *par)
{
  Double_t xx =x[0];
  Double_t funbeta1fu=0.;
  
  if (xx>0){ //allowed spectrum
     double all=sqrt(xx*(xx+2.*par[1]))*(xx + par[1])*pow(par[0]-xx,2)*fermi(par[2],xx);
     double w=xx/par[1]+1.;
     double pel=TMath::Sqrt(w*w-1.);
     double pnu=(par[0]-xx)/par[1];
     // calculation of the screened lambda2 value by interpolation 
     // of the table  with the help 
     // of the divdif CERN function for logarithms of p
     double pellog=log(pel); 

     double scrl2=divdif(plog69,pellog);
     double cf1=pnu*pnu+scrl2*pel*pel;
 
     //  correction factor 2 (empirical)
     double cf2=1.+fC1/w+fC2*w+fC3*w*w+fC4*w*w*w;
     // spectrum with correction
     funbeta1fu=all*cf1*cf2;
  }
     return funbeta1fu;
}  

Double_t Decay::funbeta2f(Double_t *x, Double_t *par)
{
  Double_t xx =x[0];
  Double_t funbeta2f=0.;
  if (xx>0){
     //allowed spectrum
     float all=sqrt(xx*(xx+2.*par[1]))*(xx + par[1])*pow(par[0]-xx,2)*fermi(par[2],xx);
     //correction factor 
     float w=xx/par[1]+1.;
     float pel=sqrt(w*w-1.);
     float pnu=(par[0]-xx)/par[1];
     float cf=1.;
     if(fKf==1) cf=pel*pel+fC1*pnu*pnu;
     if(fKf==2) cf=pow(pel,4)+fC1*pel*pel*pnu*pnu+fC2*pow(pnu,4);
     if(fKf==3) cf=pow(pel,6)+fC1*pow(pel,4)*pnu*pnu+
                  fC2*pel*pel*pow(pnu,4)+fC3*pow(pnu,6);
     if(fKf==4) cf=pow(pel,8)+fC1*pow(pel,6)*pnu*pnu+
                  fC2*pow(pel,4)*pow(pnu,4)+fC3*pel*pel*pow(pnu,6)+
                  fC4*pow(pnu,8);
     //spectrum with correction 
     funbeta2f=all*cf;
  } 
     return funbeta2f; 
}
               
void Decay::GENBBdia()
{
   printf("DECAY0: Generation of events of decay of natural radioactive\n");
   printf("isotopes and various modes of double beta decay\n");
   printf("DECAY units: \n");
   printf("   energy    - MeV\n   momentum  - MeV/c \n   time      - sec\n   angle     - degree\n");
   printf("Which type of events do you want to generate:\n");
   printf(" 1. double beta processes\n");
   printf(" 2. internal or external background or calibration sources");
   printf("\n?");
   scanf("%d",&i2bbs);
   if (i2bbs==1){
      printf("Double beta nuclides:\n");
      printf("  Ca48\n  Ni58\n  Zn64   Zn70\n  Ge76\n  Se74   Se82\n");
      printf("  Sr84\n  Zr94   Zr96\n  Mo92   Mo100\n  Ru96   Ru104\n");
      printf("  Cd106  Cd108  Cd114  Cd116\n  Sn112  Sn122  Sn124\n");
      printf("  Te120  Te128  Te130\n  Xe136\n  Ce136  Ce138  Ce142\n");
      printf("  Nd148  Nd150\n  W180   W186\n  Bi214+At214\n  Pb214+Po214\n");
      printf("  Po218+Rn218+Po214\n  Rn222+Ra222+Rn218+Po214\n");
			printf("Double beta nuclides(update Fieguth/Althueser):\n");
			printf("  Xe124\n");
      printf("\n?");
      scanf("%s",chn);
      sprintf(chnuclide,"%s",chn);
      if (strncmp(chn,"Ca48",4)==0 || strncmp(chn,"CA48",4)==0 || strncmp(chn,"ca48",4)==0){
         do{
            printf("48-Ti level:    0. 0+ (gs)     {0 MeV}\n");
            printf("                1. 2+ (1)  {0.984 MeV}\n");
            printf("                2. 2+ (2)  {2.421 MeV}");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>2);
      }
     else if(strncmp(chn,"Ni58",4)==0 || strncmp(chn,"NI58",4)==0 || strncmp(chn,"ni58",4)==0){
         do{
            printf("58-Fe level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.811 MeV}\n");
            printf("                2. 2+ (2)  {1.675 MeV}");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>2);
     }
     else if(strncmp(chn,"Zn64",4)==0 || strncmp(chn,"ZN64",4)==0 || strncmp(chn,"zn64",4)==0){
            printf("transition is possible only to 64-Ni 0+ (gs) level\n");
            ilevel=0;
     }
     else if(strncmp(chn,"Zn70",4)==0 || strncmp(chn,"ZN70",4)==0 || strncmp(chn,"zn70",4)==0){
            printf("transition is possible only to 70-Ge 0+ (gs) level\n");
            ilevel=0;
     }
     else if(strncmp(chn,"Ge76",4)==0 || strncmp(chn,"GE76",4)==0 || strncmp(chn,"ge76",4)==0){
         do{
            printf("76-Se level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.559 MeV}\n");
            printf("                2. 0+ (1)  {1.122 MeV}\n");
            printf("                3. 2+ (2)  {1.216 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>3);
     }
     else if(strncmp(chn,"Se74",4)==0 || strncmp(chn,"SE74",4)==0 || strncmp(chn,"se74",4)==0){
         do{
            printf("74-Ge level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.596 MeV}\n");
            printf("                2. 2+ (2)  {1.204 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>2);
     }
     else if(strncmp(chn,"Se82",4)==0 || strncmp(chn,"SE82",4)==0 || strncmp(chn,"se82",4)==0){
         do{
            printf("82-Kr level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.776 MeV}\n");
            printf("                2. 2+ (2)  {1.475 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>2);
     }
     else if(strncmp(chn,"Sr84",4)==0 || strncmp(chn,"SR84",4)==0 || strncmp(chn,"sr84",4)==0){
         do{
            printf("84-Kr level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.882 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>1);
     }
     else if(strncmp(chn,"Zr94",4)==0 || strncmp(chn,"ZR94",4)==0 || strncmp(chn,"zr94",4)==0){
         do{
            printf("94-Mo level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.871 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>1);
     }
     else if(strncmp(chn,"Zr96",4)==0 || strncmp(chn,"ZR96",4)==0 || strncmp(chn,"zr96",4)==0){
         do{
            printf("96-Mo level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.778 MeV}\n");
            printf("                2. 0+ (1)  {1.148 MeV}\n");
            printf("                3. 2+ (2)  {1.498 MeV}\n");
            printf("                4. 2+ (3)  {1.626 MeV}\n");
            printf("                5. 2+ (4)  {2.096 MeV}\n");
            printf("                6. 2+ (5)  {2.426 MeV}\n");
            printf("                7. 0+ (2)  {2.623 MeV}\n");
            printf("                8. 2+ (6)  {2.700 MeV}\n");
            printf("                9. 2+?(7)  {2.713 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>9);
     }
     else if(strncmp(chn,"Mo92",4)==0 || strncmp(chn,"MO92",4)==0 || strncmp(chn,"mo92",4)==0){
         do{
            printf("92-Zr level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.934 MeV}\n");
            printf("                2. 0+ (1)  {1.383 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>2);
     }
     else if(strncmp(chn,"Mo100",5)==0 || strncmp(chn,"MO100",5)==0 || strncmp(chn,"mo100",5)==0){
         do{
            printf("100-Ru level:    0. 0+ (gs) {0 MeV}\n");
            printf("                 1. 2+ (1)  {0.540 MeV}\n");
            printf("                 2. 0+ (1)  {1.130 MeV}\n");
            printf("                 3. 2+ (2)  {1.362 MeV}\n");
            printf("                 4. 0+ (2)  {1.741 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>4);
     }
     else if(strncmp(chn,"Ru96",4)==0 || strncmp(chn,"RU96",4)==0 || strncmp(chn,"ru96",4)==0){
         do{
            printf("96-Mo level:    0. 0+ (gs) {0 MeV}\n");
            printf("                1. 2+ (1)  {0.778 MeV}\n");
            printf("                3. 2+ (2)  {1.498 MeV}\n");
            printf("                4. 2+ (3)  {1.626 MeV}\n");
            printf("                5. 2+ (4)  {2.096 MeV}\n");
            printf("                6. 2+ (5)  {2.426 MeV}\n");
	    printf("                7. 0+ (2)  {2.623 MeV}\n");
	    printf("                8. 2+ (6)  {2.700 MeV}\n");
            printf("                9. 2+?(7)  {2.713 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>9);
     }
     else if(strncmp(chn,"Ru104",5)==0 || strncmp(chn,"RU104",5)==0 || strncmp(chn,"ru104",5)==0){
         do{
            printf("104-Pd level:    0. 0+ (gs) {0 MeV}\n");
            printf("                 1. 2+ (1)  {0.556 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>1);
      }
     else if(strncmp(chn,"Cd106",5)==0 || strncmp(chn,"CD106",5)==0 || strncmp(chn,"cd106",5)==0){
         do{
            printf("106-Pd level:    0. 0+ (gs) {0 MeV}\n");
            printf("                 1. 2+ (1)  {0.512 MeV}\n");
            printf("                 2. 2+ (2)  {1.128 MeV}\n");
            printf("                 3. 0+ (1)  {1.134 MeV}\n");
	    printf("                 4. 2+ (3)  {1.562 MeV}\n");
	    printf("                 5. 0+ (2)  {1.706 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>5);
      }
     else if(strncmp(chn,"Cd108",5)==0 || strncmp(chn,"CD108",5)==0 || strncmp(chn,"cd108",5)==0){
            printf("transition is possible only to 108-Pd 0+(gs) level\n");
            ilevel=0;
   }
     else if(strncmp(chn,"Cd114",5)==0 || strncmp(chn,"CD114",5)==0 || strncmp(chn,"cd114",5)==0){
            printf("transition is possible only to 114-Sn 0+(gs) level\n");
            ilevel=0;
   }
     else if(strncmp(chn,"Cd116",5)==0 || strncmp(chn,"CD116",5)==0 || strncmp(chn,"cd116",5)==0){
         do{
            printf("116-Sn level:   0. 0+ (gs)     {0 MeV}\n");
	    printf("                1. 2+ (1)  {1.294 MeV}\n");
	    printf("                2. 0+ (1)  {1.757 MeV}\n");
	    printf("                3. 0+ (2)  {2.027 MeV}\n");
	    printf("                4. 2+ (2)  {2.112 MeV}\n");
      	    printf("                5. 2+ (3)  {2.225 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>5);
      }
      else if(strncmp(chn,"Sn112",5)==0 || strncmp(chn,"SN112",5)==0 || strncmp(chn,"sn112",5)==0){
         do{
            printf("112-Cd level:   0. 0+ (gs)     {0 MeV}\n");
	    printf("                1. 2+ (1)  {0.618 MeV}\n");
            printf("                2. 0+ (1)  {1.224 MeV}\n");
	    printf("                3. 2+ (2)  {1.312 MeV}\n");
	    printf("                4. 0+ (2)  {1.433 MeV}\n");
	    printf("                5. 2+ (3)  {1.469 MeV}\n");
	    printf("                6. 0+ (3)  {1.871 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>6);
      }
      else if(strncmp(chn,"Sn122",5)==0 || strncmp(chn,"SN122",5)==0 || strncmp(chn,"sn122",5)==0){
            printf("transition is possible only to 122-Te 0+(gs) level");
	    ilevel=0;
       }
       else if(strncmp(chn,"Sn124",5)==0 || strncmp(chn,"SN124",5)==0 || strncmp(chn,"sn124",5)==0){
         do{
            printf("124-Te level:   0. 0+ (gs)     {0 MeV}\n");
	    printf("                1. 2+ (1)  {0.603 MeV}\n");
	    printf("                2. 2+ (2)  {1.326 MeV}\n");
	    printf("                3. 0+ (1)  {1.657 MeV}\n");
	    printf("                4. 0+ (2)  {1.833 MeV}\n");
	    printf("                5. 2+ (3)  {2.039 MeV}\n");
	    printf("                6. 2+ (4)  {2.092 MeV}\n");
	    printf("                7. 0+ (3)  {2.153 MeV}\n");
	    printf("                8. 2+ (5)  {2.182 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>8);
      }
      else if(strncmp(chn,"Te120",5)==0 || strncmp(chn,"TE120",5)==0 || strncmp(chn,"te120",5)==0){
         do{
            printf("120-Sn level:   0. 0+ (gs)     {0 MeV}\n");
	    printf("                1. 2+ (1)  {1.171 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>1);
      }
      else if(strncmp(chn,"Te128",5)==0 || strncmp(chn,"TE128",5)==0 || strncmp(chn,"te128",5)==0){
         do{
            printf("128-Xe level:   0. 0+ (gs)     {0 MeV}\n");
	    printf("                1. 2+ (1)  {0.443 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>1);
      }
      else if(strncmp(chn,"Te130",5)==0 || strncmp(chn,"TE130",5)==0 || strncmp(chn,"te130",5)==0){
         do{
            printf("130-Xe level:   0. 0+ (gs)     {0 MeV}\n");
	    printf("                1. 2+ (1)  {0.536 MeV}\n");
	    printf("                2. 2+ (2)  {1.122 MeV}\n");
	    printf("                3. 0+ (1)  {1.794 MeV}\n");
            printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>3);
      }
			else if(strncmp(chn,"Xe124",5)==0 || strncmp(chn,"XE124",5)==0 || strncmp(chn,"xe124",5)==0){
         do{
             printf("124-Te level:   0. 0+ (gs)     {0 MeV}\n");
			 printf("                1. 2+ (1)  {0.819 MeV}\n");
			 printf("                2. 2+ (2)  {1.551 MeV}\n");
			 printf("                3. 0+ (1)  {1.579 MeV}\n");
	     //printf("                4. 2+ (3)  {2.080 MeV}\n");
	     //printf("                5. 2+ (4)  {2.129 MeV}\n");
	     //printf("                6. 0+ (2)  {2.141 MeV}\n");
	     //printf("                7. 2+ (5)  {2.223 MeV}\n");
	     //printf("                8. 0+ (3)  {2.315 MeV}\n");
	     //printf("                9. 2+ (6)  {2.400 MeV}\n");  
             printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>3);
      }
      else if(strncmp(chn,"Xe136",5)==0 || strncmp(chn,"XE136",5)==0 || strncmp(chn,"xe136",5)==0){
         do{
             printf("136-Ba level:   0. 0+ (gs)     {0 MeV}\n");
	     printf("                1. 2+ (1)  {0.819 MeV}\n");
	     printf("                2. 2+ (2)  {1.551 MeV}\n");
	     printf("                3. 0+ (1)  {1.579 MeV}\n");
	     printf("                4. 2+ (3)  {2.080 MeV}\n");
	     printf("                5. 2+ (4)  {2.129 MeV}\n");
	     printf("                6. 0+ (2)  {2.141 MeV}\n");
	     printf("                7. 2+ (5)  {2.223 MeV}\n");
	     printf("                8. 0+ (3)  {2.315 MeV}\n");
	     printf("                9. 2+ (6)  {2.400 MeV}\n");  
             printf("\n?");
            scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>9);
      }
      else if(strncmp(chn,"Ce136",5)==0 || strncmp(chn,"CE136",5)==0 || strncmp(chn,"ce136",5)==0){
         do{
             printf("136-Ba level:   0. 0+ (gs)     {0 MeV}\n");
             printf("                1. 2+ (1)  {0.819 MeV}\n");
             printf("                2. 2+ (2)  {1.551 MeV}\n");
             printf("                3. 0+ (1)  {1.579 MeV}\n");
             printf("                4. 2+ (3)  {2.080 MeV}\n");
             printf("                5. 2+ (4)  {2.129 MeV}\n");
             printf("                6. 0+ (2)  {2.141 MeV}\n");
             printf("                7. 2+ (5)  {2.223 MeV}\n");
             printf("                8. 0+ (3)  {2.315 MeV}\n");
             printf("                9. 2+ (6)  {2.400 MeV}\n");
             printf("\n?");
             scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>9);
      }
      else if(strncmp(chn,"Ce138",5)==0 || strncmp(chn,"CE138",5)==0 || strncmp(chn,"ce138",5)==0){
             printf("transition is possible only to 138-Ba 0+(gs) level\n");
	     ilevel=0;
      } 
      else if(strncmp(chn,"Ce142",5)==0 || strncmp(chn,"CE142",5)==0 || strncmp(chn,"ce142",5)==0){
             printf("transition is possible only to 142-Nd 0+(gs) level\n");
	     ilevel=0;
      }
      else if(strncmp(chn,"Nd148",5)==0 || strncmp(chn,"ND148",5)==0 || strncmp(chn,"nd148",5)==0){
         do{
             printf("148-Sm level:    0. 0+ (gs)     {0 MeV}\n");
	     printf("                 1. 2+ (1)  {0.550 MeV}\n");
	     printf("                 2. 2+ (2)  {1.455 MeV}\n");
             printf("\n?");
             scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>2);
      }
      else if(strncmp(chn,"Nd150",5)==0 || strncmp(chn,"ND150",5)==0 || strncmp(chn,"nd150",5)==0){
         do{
             printf("150-Sm level:   0. 0+ (gs)     {0 MeV}\n");
	     printf("                1. 2+ (1)  {0.334 MeV}\n");
	     printf("                2. 0+ (1)  {0.740 MeV}\n");
	     printf("                3. 2+ (2)  {1.046 MeV}\n");
	     printf("                4. 2+ (3)  {1.194 MeV}\n");
	     printf("                5. 0+ (2)  {1.256 MeV}\n");
             printf("\n?");
             scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>5);
      }
      else if(strncmp(chn,"W180",4)==0 ||  strncmp(chn,"w180",4)==0){
             printf("transition is possible only to 180-Hf 0+(gs) level\n");
	     ilevel=0;
      }
      else if(strncmp(chn,"W186",4)==0 ||  strncmp(chn,"w186",4)==0){
           do{
              printf("186-Os level:   0. 0+ (gs)     {0 MeV}\n");
	      printf("                1. 2+ (1)  {0.137 MeV}\n");
              printf("\n?");
              scanf("%d",&ilevel);
         } while(ilevel<0 || ilevel>1);
      }
      else if(strncmp(chn,"Bi214",5)==0 || strncmp(chn,"BI214",5)==0 || strncmp(chn,"bi214",5)==0){
              printf("transition is possible only to 214-At 1-(gs) level\n");
	      ilevel=0;
      }
      else if(strncmp(chn,"Pb214",5)==0 || strncmp(chn,"PB214",5)==0 || strncmp(chn,"pb214",5)==0){
              printf("transition is possible only to 214-Po 0+(gs) level\n");
	      ilevel=0;
      }
      else if(strncmp(chn,"Po218",5)==0 || strncmp(chn,"PO218",5)==0 || strncmp(chn,"po218",5)==0){
              printf("transition is possible only to 218-Rn 0+(gs) level\n");
	      ilevel=0;
      }
      else if(strncmp(chn,"Rn222",5)==0 || strncmp(chn,"RN222",5)==0 || strncmp(chn,"rn222",5)==0){
              printf("transition is possible only to 222-Ra 0+(gs) level\n");
	      ilevel=0;
      }
      else printf(" unknown double beta nuclide\n");
 
      
      do{
    	printf(" mode of bb-decay:\n");
	printf("      1. 0nubb(mn) 0+ -> 0+     {2n}\n");
    	printf("      2. 0nubb(rc) 0+ -> 0+     {2n}\n");
    	printf("      3. 0nubb(rc) 0+ -> 0+, 2+ {N*}\n");
    	printf("      4. 2nubb     0+ -> 0+     {2n}\n");
    	printf("      5. 0nubbM1   0+ -> 0+     {2n}\n");
    	printf("      6. 0nubbM2   0+ -> 0+     (2n}\n");
    	printf("      7. 0nubbM3   0+ -> 0+     {2n}\n");
    	printf("      8. 0nubbM7   0+ -> 0+     {2n}\n");
    	printf("      9. 0nubb(rc) 0+ -> 2+     {2n}\n");
    	printf("     10. 2nubb     0+ -> 2+     {2n}, {N*}\n");
    	printf("     11. 0nuKb+    0+ -> 0+, 2+\n");
    	printf("     12. 2nuKb+    0+ -> 0+, 2+\n");
    	printf("     13. 0nu2K     0+ -> 0+, 2+\n");
    	printf("     14. 2nu2K     0+ -> 0+, 2+\n");
    	printf("     15. 2nubb     0+ -> 0+ with bosonic neutrinos\n");
    	printf("     16. 2nubb     0+ -> 2+ with bosonic neutrinos\n");
    	printf("\n   5-8: Majoron(s) with spectral index SI:\n");
    	printf("     SI=1 - old M of Gelmini-Roncadelli\n");
    	printf("     SI=2 - bulk M of Mohapatra\n");
    	printf("     SI=3 - double M, vector M, charged M\n");
    	printf("     SI=7\n");
    	printf("\n?");
     
    	scanf("%d",&modebb0);
    	if (modebb0>=1 && modebb0<=5)    fModebb=modebb0;
    	if (modebb0==6)                  fModebb=14;
    	if (modebb0==7)                  fModebb=6;
    	if (modebb0==8)                  fModebb=13;
    	if (modebb0>=9 && modebb0<=14)   fModebb=modebb0-2;
    	if (modebb0==15 || modebb0==16)  fModebb=modebb0;
    	if (modebb0>16 || modebb0<1) printf("\n Unknown mode try again\n");
      } while (modebb0<1||modebb0>16);

      ebb1=0.;
      ebb2=4.3;

      if (fModebb== 4 || fModebb== 5 || fModebb== 6|| 
          fModebb== 8 || fModebb==10|| fModebb==13 || fModebb==14){
          printf("Do you want to restrict energy range for generated particles?Y/N\n");
          scanf("%s",chn);
          if (strncmp(chn,"y",1)==0|| strncmp(chn,"Y",1)==0){
             printf(" range for sum of e-/e+ energies (MeV): \n");
             scanf("%f %f",&ebb1,&ebb2);
          }
      }
   }
   if (i2bbs==2){
      printf(" Background & sources:\n");
      printf(" Ac228          Hf182          Rh106\n");
      printf(" Ar39           I126           Sb125\n");    
      printf(" Ar42           I133           Sb126\n");    
      printf(" As79+Se79m     I134           Sb133\n");    
      printf(" Bi207+Pb207m   I135           Sr90\n");     
      printf(" Bi208          K40            Ta182\n");    
      printf(" Bi210          K42            Te133\n");    
      printf(" Bi212+Po212    Kr81           Te133m\n");         
      printf(" Bi214+Po214    Kr85           Te134\n");          
      printf(" C14            Mn54           Th234\n");          
      printf(" Ca48+Sc48      Na22           Tl207\n");          
      printf(" Cd113          P32            Tl208\n");         
      printf(" Co60           Pa234m         Xe129m\n");         
      printf(" Cs136          Pb210          Xe131m\n");          
      printf(" Cs137+Ba137m   Pb211          Xe133\n");          
      printf(" Eu147          Pb212          Xe135\n");
      printf(" Eu152          Pb214          Y88\n");
      printf(" Eu154          Ra228          Y90\n");
      printf(" Gd146          Rb87           Zn65\n");
      printf("                           Zr96+Nb96\n");
      printf(" (Art)- Artificial event (beta decay + internal e+e- + GEANT particles)\n");      
      printf(" (Com)- Compton effect \n Moller scattering\n");
      printf(" Pair creation from external gamma quanta\n");
      printf("\n?");
      scanf("%s",chn);
      sprintf(chnuclide,"%s",chn);
      printf("%s",chn);
      if (strncmp(chn,"Art",3)==0|| strncmp(chn,"ART",3)==0 || strncmp(chn,"art",3)==0){
         printf(" Emission of up to 10 beta [b], \n");
         printf(" e+e- internal conversion pairs [p],\n");
         printf(" and any of GEANT particles [G] \n"); 
         printf(" in region of energies and angles with time delay and halflife\n");
         printf(" GEANT particles: \n");
         printf(" 1 - gamma         2 - positron     3 - electron\n");
         printf(" 4 - neutrino      5 - muon+        6 - muon- \n");
         printf(" 7 - pion0         8 - pion+        9 - pion- \n");
         printf("10 - kaon0 long   11 - kaon+       12 - kaon- \n");
         printf("13 - neutron      14 - proton      15 - antiproton \n");
         printf("16 - kaon0 short  17 - eta         18 - lambda \n");
         printf("19 - sigma+       20 - sigma0      21 - sigma- \n");
         printf("22 - xi0          23 - xi-         24 - omega  \n");
         printf("25 - antineutron  26 - antilambda  27 - antisigma- \n");
         printf("28 - antisigma0   29 - antisigma+  30 - antixi0 \n");
         printf("31 - antixi+      32 - antiomega+  45 - deuteron \n");
         printf("46 - tritium      47 - alpha       48 - geantino \n");
         printf("49 - He3          50 - Cerenkov \n");
         printf(" number of parts in artificial event: \n");
         scanf("%d",&nartparts);
         nartparts=min(10,nartparts);
         for (int i=0;i<nartparts;i++){
             printf("%2d: identifier (b/p/G): ",i);
             scanf("%s",chn);
             if (strncmp(chn,"B",1)==0|| strncmp(chn,"b",1)==0) sprintf(chart[i],"Be");
             if (strncmp(chn,"G",1)==0|| strncmp(chn,"g",1)==0) sprintf(chart[i],"GP"); 
             if (strncmp(chn,"P",1)==0|| strncmp(chn,"p",1)==0) sprintf(chart[i],"Pi");
             if (strncmp(chn," ",1)==0) printf("unknown particle\n");
             else if (strncmp(chart[i],"Be",2)==0){
                printf("Qbeta and Z of daughter nucleus\n");
                printf(" Z>0 for beta- and  Z<0 for beta+):\n");
                scanf("%f %d",&artQb[i],&fArtZd[i]);     
                printf(" tdelay, thalf: ");
                scanf("%f %f",&arttdelay[i],&artthalf[i]);
             } 
             else if (strncmp(chart[i],"Pi",2)==0){
                printf("Epair, tdelay, thalf: \n");
                scanf("%f %f %f",&artQb[i],&arttdelay[i],&artthalf[i]);
             }
             else{
                int l;
                printf("GEANT particle number, Emin, Emax (MeV):\n");
                scanf("%d %f %f",&l,&artemin[i],&artemax[i]);
                if((l<1 || l>50 )|| (l>32 && l<45))printf("unknown particle\n");
                else nartnpg[i]=l;

                printf("tetamin, tetamax, phimin, phimax, tdelay, thalf:\n");
                scanf("%f %f %f %f %f %f",&arttmin[i],&arttmax[i],&artfmin[i],&artfmax[i],&arttdelay[i],&artthalf[i]); 
             } 
         }
      } 
      else if (strncmp(chn,"Com",3)==0|| strncmp(chn,"COM",3)==0 || strncmp(chn,"com",3)==0){
         printf("Ranges for energy and angles of initial gammas");
         printf("Emin, Emax (MeV): ");
         scanf("%f %f",&artemin[1],&artemax[1]); 
         printf("tetamin, tetamax, phimin, phimax:");
         scanf("%f %f %f %f",&arttmin[1],&arttmax[1],&artfmin[1],&artfmax[1]);
      }
      else if (strncmp(chn,"Mol",3)==0|| strncmp(chn,"MOL",3)==0 || strncmp(chn,"mol",3)==0){
         printf("Ranges for energy and angles of initial electron\n");
         printf("Emin, Emax (MeV): ");
         scanf("%f %f",&artemin[1],&artemax[1]);
         printf("tetamin, tetamax, phimin, phimax:");
         scanf("%f %f %f %f",&arttmin[1],&arttmax[1],&artfmin[1],&artfmax[1]);
         printf("lower energy threshold for emitted delta rays (MeV): ");
         scanf("%f ",&artQb[1]);
      }
      else if (strncmp(chn,"Pai",3)==0|| strncmp(chn,"PAI",3)==0 || strncmp(chn,"pai",3)==0){
         do{
             printf("Ranges for energy and angles of initial gammas\n");
             printf("Emin, Emax (MeV): ");
             scanf("%f %f",&artemin[1],&artemax[1]); 
         } while (min(artemin[1],artemax[1])>1.022);
         printf("tetamin, tetamax, phimin, phimax:");
         scanf("%f %f %f %f",&arttmin[1],&arttmax[1],&artfmin[1],&artfmax[1]);
         do{
            printf(" Z of target: ");
            scanf("%d ",&fArtZd[1]);
         } while(fArtZd[1]>1); 
          
      }
   }
   printf(" number of events to generate: ");
   scanf("%d",&nevents);
   printf(" number of first event [1]:    ");
   scanf("%d",&ievstart);
   if(ievstart<=0) ievstart=1;
   if(ievstart!=1){
     printf(" initial random integer for seed  ");
     scanf("%d",&irndmst);
   }
   sprintf(chfile,"no file");
   iwrfile=0;
   printf(" save events in file ? (Y/N)   ");
   scanf("%s",chn);
   if (strncmp(chn,"Y",1)==0||strncmp(chn,"y",1)==0){
      iwrfile=1;
      printf(" name of file:               ");
      scanf("%s",chfile); 
   }   
   istart=-1;
}


void Decay::GENBBsub()
{
  //GENBBsub generates the events of decay of natural radioactive 
  //nuclides and various modes of double beta decay. 
  //GENBB units: energy and moment - MeV and MeV/c; time - sec.
  //Energy release in double beta decay - in accordance with  
  //   G.Audi and A.H.Wapstra, Nucl. Phys. A 595(1995)409.
  ier=0;
  if (istart!=1) {
     sprintf(chn,"%s",chnuclide); 
     printf("in GENBBsub: chn= %s,chnuclide= %s, i2bbs=%d\n",chn,chnuclide,i2bbs);
     if (i2bbs==1){
        if (strncmp(chn,"Ca48",4)==0 || strncmp(chn,"CA48",4)==0 || strncmp(chn,"ca48",4)==0){
           sprintf(chnuclide,"Ca48");
           fQbb=4.272;
           fZdbb=22.0;
           fEK=0.0;
           if (ilevel<0||ilevel>2){ 
              printf("GENBBsub: illegal Ti48 level:%d",ilevel);
              ier=1;
              return;
           }
           if (ilevel==0) {fLevelE=0;   itrans02=0;}
           if (ilevel==1) {fLevelE=984; itrans02=2;}
           if (ilevel==2) {fLevelE=2421;itrans02=2;}
        }
        else if (strncmp(chn,"Ni58",4)==0 || strncmp(chn,"NI58",4)==0 || strncmp(chn,"ni58",4)==0){ 
           sprintf(chnuclide,"Ni58");
           fQbb=1.926;
           fZdbb=-26.0;
           fEK=0.007;
           if (ilevel<0|| ilevel>2){ 
              printf(" GENBBsub: illegal Fe58 level: %d",ilevel);
              ier=1;
              return;
           }
           if (ilevel==0) { fLevelE=0;itrans02=0;}
           if (ilevel==1) { fLevelE=811;itrans02=2;}
           if (ilevel==2) { fLevelE=1675;itrans02=2;} 
        }
        else if(strncmp(chn,"Zn64",4)==0 || strncmp(chn,"ZN64",4)==0 || strncmp(chn,"zn64",4)==0){
           sprintf(chnuclide,"Zn64");
           fQbb=1.096;
           fZdbb=-28.0;
           fEK=0.008;
           if (ilevel!=0){ 
              printf("GENBBsub: illegal Ni64 level: %d",ilevel);
              ier=1;
              return;
           }
           fLevelE=0;
           itrans02=0; 
        }
        else if(strncmp(chn,"Zn70",4)==0 || strncmp(chn,"ZN70",4)==0 || strncmp(chn,"zn70",4)==0){
           sprintf(chnuclide,"Zn70");
           fQbb=1.001;
           fZdbb=32.0;
           fEK=0.;
           if (ilevel!=0){
              printf(" GENBBsub: illegal Ge70 level: %d",ilevel);
              ier=1;
              return;
           }
	   fLevelE=0;
           itrans02=0;
        }
        else if (strncmp(chn,"Ge76",4)==0 || strncmp(chn,"GE76",4)==0 || strncmp(chn,"ge76",4)==0){
           sprintf(chnuclide,"Ge76");
           fQbb=2.039;
           fZdbb=34.;
           fEK=0.;
           if (ilevel<0|| ilevel>3){ 
	      printf(" GENBBsub: illegal Se76 level:%d",ilevel);
	      ier=1;
              return;
           } 
	   if (ilevel==0){ fLevelE=0;itrans02=0;}
           if (ilevel==1){ fLevelE=559;itrans02=2;}
	   if (ilevel==2){ fLevelE=1122;itrans02=0;}
	   if (ilevel==3){ fLevelE=1216;itrans02=2;}
        }
        else if (strncmp(chn,"Se74",4)==0 || strncmp(chn,"SE74",4)==0 || strncmp(chn,"se74",4)==0){
           sprintf(chnuclide,"Se74");
           fQbb=1.209;
           fZdbb=-32.;
           fEK=0.011;
           if (ilevel<0|| ilevel>2){
	      printf(" GENBBsub: illegal Ge74 level: %d",ilevel);
	      ier=1;
              return;
           }
           if (ilevel==0) { fLevelE=0;itrans02=0;}
	   if (ilevel==1) { fLevelE=596;itrans02=2;}
	   if (ilevel==2) { fLevelE=1204;itrans02=2;}
        }
        else if (strncmp(chn,"Se82",4)==0 || strncmp(chn,"SE82",4)==0 || strncmp(chn,"se82",4)==0){
           sprintf(chnuclide,"Se82");
           fQbb=2.995;
	   fZdbb=36.;
           fEK=0.;
           if (ilevel<0|| ilevel>2){
              printf(" GENBBsub: illegal Kr82 level: %d",ilevel);
              ier=1;
              return;
           }
	   if (ilevel==0) { fLevelE=0;itrans02=0;}
           if(ilevel==1) { fLevelE=776;itrans02=2;}
           if(ilevel==2) { fLevelE=1475;itrans02=2;}
        }
        else if(strncmp(chn,"Sr84",4)==0 || strncmp(chn,"SR84",4)==0 || strncmp(chn,"sr84",4)==0){
           sprintf(chnuclide,"Sr84");
           fQbb=1.787;
           fZdbb=-36.;
           fEK=0.014;
           if (ilevel<0|| ilevel>1){ 
  	      printf(" GENBBsub: illegal Kr84 level: %d",ilevel);
              ier=1;
              return;
           }
	   if (ilevel==0) { fLevelE=0;itrans02=0;}
  	   if (ilevel==1) { fLevelE=882;itrans02=2;}
        }
        else if(strncmp(chn,"Zr94",4)==0 || strncmp(chn,"ZR94",4)==0 || strncmp(chn,"zr94",4)==0){
           sprintf(chnuclide,"Zr94");
           fQbb=1.144;
           fZdbb=42.;
  	   fEK=0.;
           if (ilevel<0 || ilevel>1){ 
              printf(" GENBBsub: illegal Mo94 level: %d",ilevel);
	      ier=1;
              return;
           }
	   if (ilevel==0) { fLevelE=0;itrans02=0;}
           if (ilevel==1) { fLevelE=871;itrans02=2;}
        }
        else if (strncmp(chn,"Zr96",4)==0 || strncmp(chn,"ZR96",4)==0 || strncmp(chn,"zr96",4)==0){
           sprintf(chnuclide,"Zr96");
           fQbb=3.350;
           fZdbb=42.;
           fEK=0.;
           if (ilevel<0 || ilevel>9){ 
              printf(" GENBBsub: illegal Mo96 level: %d",ilevel);
              ier=1;
              return;
           }
           if (ilevel==0) { fLevelE=0;    itrans02=0;}
           if (ilevel==1) { fLevelE=778;  itrans02=2;}
           if (ilevel==2) { fLevelE=1148; itrans02=0;}
           if (ilevel==3) { fLevelE=1498; itrans02=2;}
           if (ilevel==4) { fLevelE=1626; itrans02=2;}
           if (ilevel==5) { fLevelE=2096; itrans02=2;}
           if (ilevel==6) { fLevelE=2426; itrans02=2;}
           if (ilevel==7) { fLevelE=2623; itrans02=0;}
           if (ilevel==8) { fLevelE=2700; itrans02=2;}
           if (ilevel==9) { fLevelE=2713; itrans02=2;}
        }
        else if (strncmp(chn,"Mo92",4)==0 || strncmp(chn,"MO92",4)==0 || strncmp(chn,"mo92",4)==0){
         sprintf(chnuclide,"Mo92");
         fQbb=1.649;
	 fZdbb=-40.;
	 fEK=0.018;
         if (ilevel<0 || ilevel>2){ 
            printf(" GENBBsub: illegal Zr92 level: %d",ilevel);
            ier=1;
            return;
         }
         if (ilevel==0) { fLevelE=0;   itrans02=0;}
         if (ilevel==1) { fLevelE=934; itrans02=2;}
         if (ilevel==2) { fLevelE=1383;itrans02=0;}
      }
      else if(strncmp(chn,"Mo100",5)==0 || strncmp(chn,"MO100",5)==0 || strncmp(chn,"mo100",5)==0){
         sprintf(chnuclide,"Mo100");
         fQbb=3.034;
	 fZdbb=44.;
         fEK=0.;
	 if (ilevel<0 || ilevel>4){
            printf(" GENBBsub: illegal Ru100 level: %d",ilevel);
            ier=1;
            return;
         }
	 if (ilevel==0) { fLevelE=0;   itrans02=0;}
         if (ilevel==1) { fLevelE=540; itrans02=2;}
         if (ilevel==2) { fLevelE=1130;itrans02=0;}
         if (ilevel==3) { fLevelE=1362;itrans02=2;}
         if (ilevel==4) { fLevelE=1741;itrans02=0;}
      } 
      else if(strncmp(chn,"Ru96",4)==0 || strncmp(chn,"RU96",4)==0 || strncmp(chn,"ru96",4)==0){
         sprintf(chnuclide,"Ru96");
         fQbb=2.718;
         fZdbb=-42.;
         fEK=0.020;
         if (ilevel<0 || ilevel>9){
	    printf(" GENBBsub: illegal Mo96 level %d",ilevel);
            ier=1;
            return;
         } 
	 if (ilevel==0) { fLevelE=0;    itrans02=0;}
         if (ilevel==1) { fLevelE=778;  itrans02=2;}
         if (ilevel==2) { fLevelE=1148; itrans02=0;}
         if (ilevel==3) { fLevelE=1498; itrans02=2;}
         if (ilevel==4) { fLevelE=1626; itrans02=2;}
         if (ilevel==5) { fLevelE=2096; itrans02=2;}
	 if (ilevel==6) { fLevelE=2426; itrans02=2;}
         if (ilevel==7) { fLevelE=2623; itrans02=0;}
         if (ilevel==8) { fLevelE=2700; itrans02=2; fEK=0.003;}
         if (ilevel==9) { fLevelE=2713; itrans02=2; fEK=0.002;}
      }
      else if(strncmp(chn,"Ru104",5)==0 || strncmp(chn,"RU104",5)==0 || strncmp(chn,"ru104",5)==0){
         sprintf(chnuclide,"Ru104");
         fQbb=1.301;
	 fZdbb=46.;
         fEK=0.;
         if (ilevel<0 || ilevel>1){
	    printf(" GENBBsub: illegal Pd104 level: %d",ilevel);
	    ier=1;
            return;
         }
	 if (ilevel==0) { fLevelE=0;   itrans02=0;}
	 if (ilevel==1) { fLevelE=556; itrans02=2;}
      }
      else if(strncmp(chn,"Cd106",5)==0 || strncmp(chn,"CD106",5)==0 || strncmp(chn,"cd106",5)==0){
         sprintf(chnuclide,"Cd106");
         fQbb=2.771;
         fZdbb=-46.;
         fEK=0.024;
	 if (ilevel<0 || ilevel>5){ 
	    printf(" GENBBsub: illegal Pd106 level: %d",ilevel);
	    ier=1;
	    return;
         }
	 if (ilevel==0) { fLevelE=0;    itrans02=0;}
	 if (ilevel==1) { fLevelE=512;  itrans02=2;}
	 if (ilevel==2) { fLevelE=1128; itrans02=2;}
	 if (ilevel==3) { fLevelE=1134; itrans02=0;}
	 if (ilevel==4) { fLevelE=1562; itrans02=2;}
	 if (ilevel==5) { fLevelE=1706; itrans02=0;}
      }
      else if(strncmp(chn,"Cd108",5)==0 || strncmp(chn,"CD108",5)==0 || strncmp(chn,"cd108",5)==0){
         sprintf(chnuclide,"Cd108");
         fQbb=0.269;
         fZdbb=-46.;
         fEK=0.024;
         if (ilevel!=0){ 
	    printf(" GENBBsub: illegal Pd108 level: %d",ilevel);
            ier=1;
            return;
         }
	 fLevelE=0;
         itrans02=0;
      }
      else if(strncmp(chn,"Cd114",5)==0 || strncmp(chn,"CD114",5)==0 || strncmp(chn,"cd114",5)==0){
         sprintf(chnuclide,"Cd114");
         fQbb=0.536;
         fZdbb=50.;
         fEK=0.;
         if (ilevel!=0){ 
            printf(" GENBBsub: illegal Sn114 level: %d",ilevel);
            ier=1;
            return;
         }
	 fLevelE=0;
         itrans02=0;
      }
      else if(strncmp(chn,"Cd116",5)==0 || strncmp(chn,"CD116",5)==0 || strncmp(chn,"cd116",5)==0){
         sprintf(chnuclide,"Cd116");
         fQbb=2.805;
         fZdbb=50.;
         fEK=0.;
         if (ilevel<0 || ilevel>5){ 
            printf(" GENBBsub: illegal Sn116 level: %d",ilevel);
            ier=1;
            return;
         }
	 if (ilevel==0) { fLevelE=0;    itrans02=0;}
	 if (ilevel==1) { fLevelE=1294; itrans02=2;}
         if (ilevel==2) { fLevelE=1757; itrans02=0;}
         if (ilevel==3) { fLevelE=2027; itrans02=0;}
         if (ilevel==4) { fLevelE=2112; itrans02=2;}
         if (ilevel==5) { fLevelE=2225; itrans02=2;}
      }
      else if(strncmp(chn,"Sn112",5)==0 || strncmp(chn,"SN112",5)==0 || strncmp(chn,"sn112",5)==0){
         sprintf(chnuclide,"Sn112");
         fQbb=1.919;
         fZdbb=-48.;
         fEK=0.027;
         if (ilevel<0 || ilevel>6){ 
            printf(" GENBBsub: illegal Cd112 level: %d",ilevel);
            ier=1;
            return;
         }
	 if (ilevel==0) { fLevelE=0;    itrans02=0;}
         if (ilevel==1) { fLevelE=618;  itrans02=2;}
         if (ilevel==2) { fLevelE=1224; itrans02=0;}
         if (ilevel==3) { fLevelE=1312; itrans02=2;}
         if (ilevel==4) { fLevelE=1433; itrans02=0;}
         if (ilevel==5) { fLevelE=1469; itrans02=2;}
         if (ilevel==6) { fLevelE=1871; itrans02=0;}
      }
      else if(strncmp(chn,"Sn122",5)==0 || strncmp(chn,"SN122",5)==0 || strncmp(chn,"sn122",5)==0){
         sprintf(chnuclide,"Sn122");
         fQbb=0.368;
         fZdbb=52.;
         fEK=0.;
   	 if (ilevel!=0){ 
            printf(" GENBBsub: illegal Te122 level: %d",ilevel);
            ier=1;
            return;
         }
	 fLevelE=0;
         itrans02=0;
      }
      else if(strncmp(chn,"Sn124",5)==0 || strncmp(chn,"SN124",5)==0 || strncmp(chn,"sn124",5)==0){
         sprintf(chnuclide,"Sn124");
         fQbb=2.288;
         fZdbb=52.;
 	 fEK=0.;
         if (ilevel< 0 || ilevel>8){ 
	    printf(" GENBBsub: illegal Te124 level: %d",ilevel);
	    ier=1;
            return;
         }
	 if (ilevel==0) { fLevelE=0;    itrans02=0;}
         if (ilevel==1) { fLevelE=603;  itrans02=2;}
         if (ilevel==2) { fLevelE=1326; itrans02=2;}
         if (ilevel==3) { fLevelE=1657; itrans02=0;}
         if (ilevel==4) { fLevelE=1883; itrans02=0;}
         if (ilevel==5) { fLevelE=2039; itrans02=2;}
         if (ilevel==6) { fLevelE=2092; itrans02=2;}
         if (ilevel==7) { fLevelE=2153; itrans02=0;}
         if (ilevel==8) { fLevelE=2182; itrans02=2;}
      }  
      else if(strncmp(chn,"Te120",5)==0 || strncmp(chn,"TE120",5)==0 || strncmp(chn,"te120",5)==0){
         sprintf(chnuclide,"Te120");
         fQbb=1.698;
         fZdbb=-50.;
         fEK=0.029;
   	 if (ilevel<0 || ilevel>1){ 
            printf(" GENBBsub: illegal Sn120 level: %d",ilevel);
	    ier=1;
	    return;
         }
	 if (ilevel==0) { fLevelE=0;    itrans02=0;}
         if (ilevel==1) { fLevelE=1171; itrans02=2;}
      }
      else if(strncmp(chn,"Te128",5)==0 || strncmp(chn,"TE128",5)==0 || strncmp(chn,"te128",5)==0){
         sprintf(chnuclide,"Te128");
         fQbb=0.867;
	 fZdbb=54.;
	 fEK=0.;
	 if (ilevel<0 || ilevel>1){ 
            printf(" GENBBsub: illegal Xe128 level: %d",ilevel);
            ier=1;
  	    return;
 	 }    
	 if (ilevel==0) { fLevelE=0;   itrans02=0;}
	 if (ilevel==1) { fLevelE=443; itrans02=2;}
      }
      else if(strncmp(chn,"Te130",5)==0 || strncmp(chn,"TE130",5)==0 || strncmp(chn,"te130",5)==0){
         sprintf(chnuclide,"Te130");
         fQbb=2.529;
         fZdbb=54.;
         fEK=0.;
         if (ilevel<0 || ilevel>3){ 
  	    printf(" GENBBsub: illegal Xe130 level: %d",ilevel);
            ier=1;
	    return;
          }
	  if (ilevel==0) { fLevelE=0;    itrans02=0;}
	  if (ilevel==1) { fLevelE=536;  itrans02=2;}
          if (ilevel==2) { fLevelE=1122; itrans02=2;}
          if (ilevel==3) { fLevelE=1794; itrans02=0;}
      }
			else if(strncmp(chn,"Xe124",5)==0 || strncmp(chn,"Xe124",5)==0 || strncmp(chn,"Xe124",5)==0){
          sprintf(chnuclide,"Xe124");
          fQbb=3.068;
          fZdbb=-56.;
          fEK=0.0278;
          if (ilevel<0 || ilevel>0){ 
             printf(" GENBBsub: illegal Te124 level: %d",ilevel);
 	     ier=1;
	     return;
          }
          if (ilevel==0) { fLevelE=0;    itrans02=0;}
		  //if (ilevel==1) { fLevelE=819;  itrans02=2;}
		  //if (ilevel==2) { fLevelE=1551; itrans02=2;}
		  //if (ilevel==3) { fLevelE=1579; itrans02=0;}
	  //if (ilevel==4) { fLevelE=2080; itrans02=2;}
	  //if (ilevel==5) { fLevelE=2129; itrans02=2;}
    //if (ilevel==6) { fLevelE=2141; itrans02=0;}
	  //if (ilevel==7) { fLevelE=2223; itrans02=2;}
	  //if (ilevel==8) { fLevelE=2315; itrans02=0;}
	  //if (ilevel==9) { fLevelE=2400; itrans02=2;}
      } 
      else if(strncmp(chn,"Xe136",5)==0 || strncmp(chn,"XE136",5)==0 || strncmp(chn,"xe136",5)==0){
          sprintf(chnuclide,"Xe136");
          fQbb=2.468;
          fZdbb=56.;
          fEK=0.;
          if (ilevel<0 || ilevel>9){ 
             printf(" GENBBsub: illegal Ba136 level: %d",ilevel);
 	     ier=1;
	     return;
          }
          if (ilevel==0) { fLevelE=0;    itrans02=0;}
          if (ilevel==1) { fLevelE=819;  itrans02=2;}
	  if (ilevel==2) { fLevelE=1551; itrans02=2;}
	  if (ilevel==3) { fLevelE=1579; itrans02=0;}
	  if (ilevel==4) { fLevelE=2080; itrans02=2;}
	  if (ilevel==5) { fLevelE=2129; itrans02=2;}
          if (ilevel==6) { fLevelE=2141; itrans02=0;}
	  if (ilevel==7) { fLevelE=2223; itrans02=2;}
	  if (ilevel==8) { fLevelE=2315; itrans02=0;}
	  if (ilevel==9) { fLevelE=2400; itrans02=2;}
      } 
      else if(strncmp(chn,"Ce136",5)==0 || strncmp(chn,"CE136",5)==0 || strncmp(chn,"ce136",5)==0){
          sprintf(chnuclide,"Ce136");
          fQbb=2.419;
          fZdbb=-56.;
          fEK=0.037;
          if (ilevel<0 || ilevel>9){ 
	     printf(" GENBBsub: illegal Ba136 level: %d",ilevel);
	     ier=1;
	     return;
          }
	  if (ilevel==0) { fLevelE=0;    itrans02=0;}
	  if (ilevel==1) { fLevelE=819;  itrans02=2;}
	  if (ilevel==2) { fLevelE=1551; itrans02=2;}
	  if (ilevel==3) { fLevelE=1579; itrans02=0;}
	  if (ilevel==4) { fLevelE=2080; itrans02=2;}
	  if (ilevel==5) { fLevelE=2129; itrans02=2;}
	  if (ilevel==6) { fLevelE=2141; itrans02=0;}
	  if (ilevel==7) { fLevelE=2223; itrans02=2;}
	  if (ilevel==8) { fLevelE=2315; itrans02=0;}
	  if (ilevel==9) { fLevelE=2400; itrans02=2; fEK=0.006;}
      }
      else if(strncmp(chn,"Ce138",5)==0 || strncmp(chn,"CE138",5)==0 || strncmp(chn,"ce138",5)==0){
          sprintf(chnuclide,"Ce138");
          fQbb=0.693;
	  fZdbb=-56.;
	  fEK=0.037;
	  if (ilevel!=0){ 
	     printf(" GENBBsub: illegal Ba138 level: %d",ilevel);
	     ier=1;
	     return;
           }
           fLevelE=0;
           itrans02=0;
      }
      else if(strncmp(chn,"Ce142",5)==0 || strncmp(chn,"CE142",5)==0 || strncmp(chn,"ce142",5)==0){
           sprintf(chnuclide,"Ce142");
           fQbb=1.417;
	   fZdbb=60.;
	   fEK=0.;
	   if (ilevel!=0){ 
	      printf(" GENBBsub: illegal Nd142 level:%d",ilevel);
	      ier=1;
	      return;
           }
	   fLevelE=0;
           itrans02=0;
      }
      else if(strncmp(chn,"Nd148",5)==0 || strncmp(chn,"ND148",5)==0 || strncmp(chn,"nd148",5)==0){
           sprintf(chnuclide,"Nd148");
           fQbb=1.929;
	   fZdbb=62.;
	   fEK=0.;
	   if (ilevel<0 || ilevel>2){ 
	      printf(" GENBBsub: illegal Sm148 level: %d",ilevel);
	      ier=1;
	      return;
	   }  
	   if (ilevel==0) { fLevelE=0;    itrans02=0;}
	   if (ilevel==1) { fLevelE=550;  itrans02=2;}
	   if (ilevel==2) { fLevelE=1455; itrans02=2;}
      } 
      else if(strncmp(chn,"Nd150",5)==0 || strncmp(chn,"ND150",5)==0 || strncmp(chn,"nd150",5)==0){
           sprintf(chnuclide,"Nd150");
           fQbb=3.367;
	   fZdbb=62.;
	   fEK=0.;
	   if (ilevel<0 || ilevel>5){ 
	      printf(" GENBBsub: illegal Sm150 level: %d",ilevel);
	      ier=1;
	      return;
	   }  
	   if (ilevel==0) { fLevelE=0;    itrans02=0;}
	   if (ilevel==1) { fLevelE=334;  itrans02=2;}
	   if (ilevel==2) { fLevelE=740;  itrans02=0;}
	   if (ilevel==3) { fLevelE=1046; itrans02=2;}
	   if (ilevel==4) { fLevelE=1194; itrans02=2;}
	   if (ilevel==5) { fLevelE=1256; itrans02=0;}
      }
      else if(strncmp(chn,"W180",4)==0 ||  strncmp(chn,"w180",4)==0){ 
           sprintf(chnuclide,"W180");
           fQbb=0.144;
	   fZdbb=-72.;
	   fEK=0.065;
	   if (ilevel!=0){ 
	      printf(" GENBBsub: illegal Hf180 level: %d",ilevel);
	      ier=1;
	      return;
	   }   
	   fLevelE=0;
           itrans02=0;
      }
      else if(strncmp(chn,"W186",4)==0 ||  strncmp(chn,"w186",4)==0){
           sprintf(chnuclide,"W186");
            fQbb=0.490;
	    fZdbb=76.;
	    fEK=0.;
	    if (ilevel<0 || ilevel>1){ 
	       printf(" GENBBsub: illegal Os186 level: %d",ilevel);
	       ier=1;
	       return;
	    }  
	    if (ilevel==0) { fLevelE=0;   itrans02=0;}
	    if (ilevel==1) { fLevelE=137; itrans02=2;}
      }
      else if(strncmp(chn,"Bi214",5)==0 || strncmp(chn,"BI214",5)==0 || strncmp(chn,"bi214",5)==0){
            sprintf(chnuclide,"Bi214");
            fQbb=2.180;
	    fZdbb=85.;
	    fEK=0.;
	    if (ilevel!=0){ 
	       printf(" GENBBsub: illegal At214 level: %d",ilevel);
	       ier=1;
	       return;
	    }  
	    fLevelE=0;
	    itrans02=0;
      }
      else if(strncmp(chn,"Pb214",5)==0 || strncmp(chn,"PB214",5)==0 || strncmp(chn,"pb214",5)==0){
            sprintf(chnuclide,"Pb214");
            fQbb=4.289;
	    fZdbb=84.;
	    fEK=0.;
	    if (ilevel!=0){ 
	       printf(" GENBBsub: illegal Po214 level: %d",ilevel);
	       ier=1;
	       return;
            }
	    fLevelE=0;
	    itrans02=0;
      }
      else if(strncmp(chn,"Po218",5)==0 || strncmp(chn,"PO218",5)==0 || strncmp(chn,"po218",5)==0){
            sprintf(chnuclide,"Po218");
            fQbb=3.141;
            fZdbb=86.;
	    fEK=0.;
	    if (ilevel!=0){ 
	       printf(" GENBBsub: illegal Rn218 level: %d",ilevel);
	       ier=1;
	       return;
            }
	    fLevelE=0;
	    itrans02=0; 
      }
      else if(strncmp(chn,"Rn222",5)==0 || strncmp(chn,"RN222",5)==0 || strncmp(chn,"rn222",5)==0){
            sprintf(chnuclide,"Rn222");
            fQbb=2.052;
	    fZdbb=88.;
	    fEK=0.;
	    if (ilevel!=0){ 
	       printf(" GENBBsub: illegal Ra222 level: %d",ilevel);
	       ier=1;
	       return;
	    }  
	    fLevelE=0;
	    itrans02=0;
      }
      else{
            printf(" unknown double beta nuclide\n");
	    ier=1;
	    return;
      }
      
      if (itrans02==0) { sprintf(chdspin,"0+");}
      if (itrans02==2) { sprintf(chdspin,"2+");}
      if (itrans02==0 && strncmp(chnuclide,"Bi214",5)==0) { sprintf(chdspin,"1-");}

      if (fModebb<1 || fModebb>16){ 
	 printf(" unknown bb mode:%d\n",fModebb);
	 ier=1;
	 return;
      }
      else {
	  if (fModebb==1)   sprintf(chmodebb,"0nubb(mn) 0+ -> 0+     {2n}");
	  if (fModebb==2)   sprintf(chmodebb,"0nubb(rc) 0+ -> 0+     {2n}");
          if (fModebb==3)   sprintf(chmodebb,"0nubb(rc) 0+ -> 0+, 2+ {N*}");
	  if (fModebb==4)   sprintf(chmodebb,"2nubb     0+ -> 0+     {2n}");
	  if (fModebb==5)   sprintf(chmodebb,"0nubbM1   0+ -> 0+     {2n}");
	  if (fModebb==6)   sprintf(chmodebb,"0nubbM3   0+ -> 0+     {2n}");
	  if (fModebb==7)   sprintf(chmodebb,"0nubb(rc) 0+ -> 2+     {2n}");
	  if (fModebb==8)   sprintf(chmodebb,"2nubb     0+ -> 2+     {2n}, {N*}");
	  if (fModebb==9)   sprintf(chmodebb,"0nuKb+    0+ -> 0+, 2+");
	  if (fModebb==10)  sprintf(chmodebb,"2nuKb+    0+ -> 0+, 2+");
	  if (fModebb==11)  sprintf(chmodebb,"0nu2K     0+ -> 0+, 2+");
	  if (fModebb==12)  sprintf(chmodebb,"2nu2K     0+ -> 0+, 2+");
	  if (fModebb==13)  sprintf(chmodebb,"0nubbM7   0+ -> 0+     {2n}");
	  if (fModebb==14)  sprintf(chmodebb,"0nubbM2   0+ -> 0+     {2n}");
	  if (fModebb==15)  sprintf(chmodebb,"2nubb     0+ -> 0+  bosonic nu");
	  if (fModebb==16)  sprintf(chmodebb,"2nubb     0+ -> 2+  bosonic nu");
      }   
      El=fLevelE/1000.;

      if (fZdbb>=0.)  e0=fQbb;
      if (fZdbb<0.)   e0=fQbb-4.*emass;

      if (fModebb== 9 || fModebb==10) e0=fQbb-fEK-2.*emass;
      if (fModebb==11 || fModebb==12) e0=fQbb-2.*fEK;

      if (e0<=El){ 
	      printf(" not enough energy for transition to this level: ");
	      printf(" full energy release and Elevel: %f %f",e0,El);
	      ier=1;
	      return;
      }	  
      m=fModebb;

      if (fZdbb>=0 && (m==9 || m==10 || m==11 || m==12)){
	      printf("decay mode and nuclide are inconsistent:%s %s\n ",chmodebb,chnuclide);
	      ier=1;
	      return;
      }
      if ((itrans02==0)&& (m==7 || m==8 || m==16)){
	      printf("decay mode and spin of daughter nucleus level are inconsistent:\n%s and  %s\n",chmodebb,chdspin);
	      ier=1;
	      return;
      }
      if((itrans02==2) && (m!=3&& m!=7&& m!=8&& m!=9&& m!=10&& m!=11&& m!=12&& m!=16)){
	      printf("decay mode and spin of daughter nucleus level are inconsistent:\n%s and  %s\n",chmodebb,chdspin);
	      ier=1;
	      return;

      }
   }
   if (i2bbs==2){
      printf("isotop: %s\n",chn); 
      if (strncmp(chn,"Ac228",5)==0 || strncmp(chn,"AC228",5)==0 || strncmp(chn,"ac228",5)==0){
         sprintf(chnuclide,"Ac228");
         fZdtr=90;
         fEbindeK=0.110;
         printf("%s\n",chnuclide);
      }
      else if (strncmp(chn,"Ar39",4)==0 || strncmp(chn,"AR39",4)==0 || strncmp(chn,"ar39",4)==0){
         sprintf(chnuclide,"Ar39");
         fZdtr=19;
      }
      else if (strncmp(chn,"Ar42",4)==0 || strncmp(chn,"AR42",4)==0 || strncmp(chn,"ar42",4)==0){
         sprintf(chnuclide,"Ar42");
         fZdtr=19;
      }
      else if (strncmp(chn,"As79",4)==0 || strncmp(chn,"AS79",4)==0 || strncmp(chn,"as79",4)==0){
         sprintf(chnuclide,"As79");
         fZdtr=34;
      }
      else if (strncmp(chn,"Bi207",5)==0 || strncmp(chn,"BI207",5)==0 || strncmp(chn,"bi207",5)==0){
         sprintf(chnuclide,"Bi207");
         fZdtr=-82; 
      }
      else if (strncmp(chn,"Bi208",5)==0 || strncmp(chn,"BI208",5)==0 || strncmp(chn,"bi208",5)==0){
         sprintf(chnuclide,"Bi208");

      }
      else if (strncmp(chn,"Bi210",5)==0 || strncmp(chn,"BI210",5)==0 || strncmp(chn,"bi210",5)==0){
         sprintf(chnuclide,"Bi210");
         fZdtr=84; 
      }
      else if (strncmp(chn,"Bi212",5)==0 || strncmp(chn,"BI212",5)==0 || strncmp(chn,"bi212",5)==0){
         sprintf(chnuclide,"Bi212");
         fZdtr=84; 
      }
      else if (strncmp(chn,"Bi214",5)==0 || strncmp(chn,"BI214",5)==0 || strncmp(chn,"bi214",5)==0){
         sprintf(chnuclide,"Bi214");
         fZdtr=84; 
      }
      else if (strncmp(chn,"C14",3)==0 || strncmp(chn,"c14",3)==0){
         sprintf(chnuclide,"C14");
         fZdtr=7; 
      }
      else if (strncmp(chn,"Ca48",4)==0 || strncmp(chn,"CA48",4)==0 || strncmp(chn,"ca48",4)==0){
         sprintf(chnuclide,"Ca48");
         fZdtr=21; 
      }
      else if (strncmp(chn,"Co60",4)==0 || strncmp(chn,"CO60",4)==0 || strncmp(chn,"co60",4)==0){
         sprintf(chnuclide,"Co60");
         fZdtr=28; 
      }
      else if (strncmp(chn,"Cd113",5)==0 || strncmp(chn,"CD113",5)==0 || strncmp(chn,"cd113",5)==0){
         sprintf(chnuclide,"Cd113");
         fZdtr=49; 
      }
      else if (strncmp(chn,"Cs136",5)==0 || strncmp(chn,"CS136",5)==0 || strncmp(chn,"cs136",5)==0){
         sprintf(chnuclide,"Cs136");
         fZdtr=56; 
      }
      else if (strncmp(chn,"Cs137",5)==0 || strncmp(chn,"CS137",5)==0 || strncmp(chn,"cs137",5)==0){
         sprintf(chnuclide,"Cs137");
         fZdtr=56;
      }
      else if (strncmp(chn,"Eu147",5)==0 || strncmp(chn,"EU147",5)==0 || strncmp(chn,"eu147",5)==0){
         sprintf(chnuclide,"Eu147");
         fZdtr=-62;
      }
      else if (strncmp(chn,"Eu152",5)==0 || strncmp(chn,"EU152",5)==0 || strncmp(chn,"eu152",5)==0){
         sprintf(chnuclide,"Eu152");
         fZdtr=64;
      }
      else if (strncmp(chn,"Eu154",5)==0 || strncmp(chn,"EU154",5)==0 || strncmp(chn,"eu154",5)==0){
         sprintf(chnuclide,"Eu154");
         fZdtr=64;
      }
      else if (strncmp(chn,"Gd146",5)==0 || strncmp(chn,"GD146",5)==0 || strncmp(chn,"gd146",5)==0){
         sprintf(chnuclide,"Gd146");
         fZdtr=-63;
      }
      else if (strncmp(chn,"Hf182",5)==0 || strncmp(chn,"HF182",5)==0 || strncmp(chn,"hf182",5)==0){
         sprintf(chnuclide,"Hf182");
         fZdtr=73;
      }
      else if (strncmp(chn,"I126",4)==0 || strncmp(chn,"i126",4)==0){
         sprintf(chnuclide,"I126");
         fZdtr=54;
      }
      else if (strncmp(chn,"I133",4)==0 || strncmp(chn,"i133",4)==0){
         sprintf(chnuclide,"I133");
         fZdtr=54;
      }
      else if (strncmp(chn,"I134",4)==0 || strncmp(chn,"i134",4)==0){
         sprintf(chnuclide,"I134");
         fZdtr=54;
      }
      else if (strncmp(chn,"I135",4)==0 || strncmp(chn,"i135",4)==0){
         sprintf(chnuclide,"I135");
         fZdtr=54;
      }
      else if (strncmp(chn,"K40",3)==0 || strncmp(chn,"k40",3)==0){
         sprintf(chnuclide,"K40");
         fZdtr=20;
      }
      else if (strncmp(chn,"K42",3)==0 || strncmp(chn,"k42",3)==0){
         sprintf(chnuclide,"K42");
         fZdtr=20;
      }
      else if (strncmp(chn,"Kr81",4)==0 || strncmp(chn,"KR81",4)==0 || strncmp(chn,"kr81",4)==0){
         sprintf(chnuclide,"Kr81");
      }
      else if (strncmp(chn,"Kr85",4)==0 || strncmp(chn,"KR85",4)==0 || strncmp(chn,"kr85",4)==0){
         sprintf(chnuclide,"Kr85");
         fZdtr=37;
      }
      else if (strncmp(chn,"Mn54",4)==0 || strncmp(chn,"MN54",4)==0 || strncmp(chn,"mn54",4)==0){
         sprintf(chnuclide,"Mn54");
      }
      else if (strncmp(chn,"Na22",4)==0 || strncmp(chn,"NA22",4)==0 || strncmp(chn,"na22",4)==0){
         sprintf(chnuclide,"Na22");
         fZdtr=-10;
      }
      else if (strncmp(chn,"P32",3)==0 || strncmp(chn,"p32",3)==0){
         sprintf(chnuclide,"P32");
         fZdtr=16;
      }
      else if (strncmp(chn,"Pa234m",6)==0 || strncmp(chn,"PA234m",6)==0 || strncmp(chn,"pa234m",6)==0){
         sprintf(chnuclide,"Pa234m");
         fZdtr=92;
      }
      else if (strncmp(chn,"Pb210",5)==0 || strncmp(chn,"PB210",5)==0 || strncmp(chn,"pb210",5)==0){
         sprintf(chnuclide,"Pb210");
         fZdtr=83;
      }
      else if (strncmp(chn,"Pb211",5)==0 || strncmp(chn,"PB211",5)==0 || strncmp(chn,"pb211",5)==0){
         sprintf(chnuclide,"Pb211");
         fZdtr=83;
      }
      else if (strncmp(chn,"Pb212",5)==0 || strncmp(chn,"PB212",5)==0 || strncmp(chn,"pb212",5)==0){
         sprintf(chnuclide,"Pb212");
         fZdtr=83;
      }
      else if (strncmp(chn,"Pb214",5)==0 || strncmp(chn,"PB214",5)==0 || strncmp(chn,"pb214",5)==0){
         sprintf(chnuclide,"Pb214");
         fZdtr=83;
      }
      else if (strncmp(chn,"Ra228",5)==0 || strncmp(chn,"RA228",5)==0 || strncmp(chn,"ra228",5)==0){
         sprintf(chnuclide,"Ra228");
         fZdtr=89;
      }
      else if (strncmp(chn,"Rb87",4)==0 || strncmp(chn,"RB87",4)==0 || strncmp(chn,"rb87",4)==0){
         sprintf(chnuclide,"Rb87");
         fZdtr=38;
      }
      else if (strncmp(chn,"Rh106",5)==0 || strncmp(chn,"RH106",5)==0 || strncmp(chn,"rh106",5)==0){
         sprintf(chnuclide,"Rh106");
         fZdtr=46;
      }
      else if (strncmp(chn,"Sb125",5)==0 || strncmp(chn,"SB125",5)==0 || strncmp(chn,"sb125",5)==0){
         sprintf(chnuclide,"Sb125");
         fZdtr=52;
      }
      else if (strncmp(chn,"Sb126",5)==0 || strncmp(chn,"SB126",5)==0 || strncmp(chn,"sb126",5)==0){
         sprintf(chnuclide,"Sb126");
         fZdtr=52;
      }
      else if (strncmp(chn,"Sb133",5)==0 || strncmp(chn,"SB133",5)==0 || strncmp(chn,"sb133",5)==0){
         sprintf(chnuclide,"Sb133");
         fZdtr=52;
      }
      else if (strncmp(chn,"Sr90",4)==0 || strncmp(chn,"SR90",4)==0 || strncmp(chn,"sr90",4)==0){
         sprintf(chnuclide,"Sr90");
         fZdtr=39;
      }
      else if (strncmp(chn,"Ta182",5)==0 || strncmp(chn,"TA182",5)==0 || strncmp(chn,"ta182",5)==0){
         sprintf(chnuclide,"Ta182");
         fZdtr=74;
      }
      else if (strncmp(chn,"Te133m",6)==0 || strncmp(chn,"TE133m",6)==0 || strncmp(chn,"te133m",6)==0){
         sprintf(chnuclide,"Te133m");
         fZdtr=53;
      }
      else if (strncmp(chn,"Te133",5)==0 || strncmp(chn,"TE133",5)==0 || strncmp(chn,"te133",5)==0){
         sprintf(chnuclide,"Te133");
         fZdtr=53;
      }
      else if (strncmp(chn,"Te134",5)==0 || strncmp(chn,"TE134",5)==0 || strncmp(chn,"te134",5)==0){
         sprintf(chnuclide,"Te134");
         fZdtr=53;
      }
      else if (strncmp(chn,"Th234",5)==0 || strncmp(chn,"TH234",5)==0 || strncmp(chn,"th234",5)==0){
         sprintf(chnuclide,"Th234");
         fZdtr=91;
      }
      else if (strncmp(chn,"Tl207",5)==0 || strncmp(chn,"TL207",5)==0 || strncmp(chn,"tl207",5)==0){
         sprintf(chnuclide,"Tl207");
         fZdtr=82;
      }
      else if (strncmp(chn,"Tl208",5)==0 || strncmp(chn,"TL208",5)==0 || strncmp(chn,"tl208",5)==0){
         sprintf(chnuclide,"Tl208");
         fZdtr=82;
      }
      else if (strncmp(chn,"Xe129m",6)==0 || strncmp(chn,"XE129m",6)==0 || strncmp(chn,"xe129m",6)==0){
         sprintf(chnuclide,"Xe129m");
         fZdtr=54;
      }
      else if (strncmp(chn,"Xe131m",6)==0 || strncmp(chn,"XE131m",6)==0 || strncmp(chn,"xe131m",6)==0){
         sprintf(chnuclide,"Xe131m");
         fZdtr=54;
      }
      else if (strncmp(chn,"Xe133",5)==0 || strncmp(chn,"XE133",5)==0 || strncmp(chn,"xe133",5)==0){
         sprintf(chnuclide,"Xe133");
         fZdtr=55;
      }
      else if (strncmp(chn,"Xe135",5)==0 || strncmp(chn,"XE135",5)==0 || strncmp(chn,"xe135",5)==0){
         sprintf(chnuclide,"Xe135");
         fZdtr=55;
      }
      else if (strncmp(chn,"Y88",3)==0 || strncmp(chn,"y88",3)==0 ){
         sprintf(chnuclide,"Y88");
         fZdtr=-38;
      }
      else if (strncmp(chn,"Y90",3)==0 || strncmp(chn,"y90",3)==0 ){
         sprintf(chnuclide,"Y90");
         fZdtr=40;
      }
      else if (strncmp(chn,"Zn65",4)==0 || strncmp(chn,"ZN65",4)==0 || strncmp(chn,"zn65",4)==0){
         sprintf(chnuclide,"Zn65");
      }
      else if (strncmp(chn,"Zr96",4)==0 || strncmp(chn,"ZR96",4)==0 || strncmp(chn,"zr96",4)==0){
         sprintf(chnuclide,"Zr96");
         fZdtr=41;
      }
      else if (strncmp(chn,"Art",3)==0 || strncmp(chn,"ART",3)==0 || strncmp(chn,"art",3)==0){
         sprintf(chnuclide,"Artificial");
      
         nartparts=min(10,nartparts);   
         for (int i=0;i<nartparts;++i){
             if (strncmp(chart[i],"Be",2)==0){
                artemin[i]=artQb[i];
                artemax[i]=artQb[i];
	        arttmin[i]=0.;
                arttmax[i]=pi;
                artfmin[i]=0.;
                artfmax[i]=2.*pi;
             }    
             else if (strncmp(chart[i],"Pi",2)==0){
                artemin[i]=artQb[i];
                artemax[i]=artQb[i];
                arttmin[i]=0.;
                arttmax[i]=pi;
                artfmin[i]=0.;
                artfmax[i]=2.*pi;
             }
             else{
                arttmin[i]=arttmin[i]/180.*pi;
	        arttmax[i]=arttmax[i]/180.*pi;
	        artfmin[i]=artfmin[i]/180.*pi;
	        artfmax[i]=artfmax[i]/180.*pi;
             }
         }
      }
      else if (strncmp(chn,"Com",3)==0 || strncmp(chn,"COM",3)==0 || strncmp(chn,"com",3)==0){
         sprintf(chnuclide,"Compton");
         arttmin[1]=arttmin[1]/180.*pi;
	 arttmax[1]=arttmax[1]/180.*pi;
	 artfmin[1]=artfmin[1]/180.*pi;
	 artfmax[1]=artfmax[1]/180.*pi;
      }
      else if (strncmp(chn,"Mol",3)==0 || strncmp(chn,"MOL",3)==0 || strncmp(chn,"mol",3)==0){
         sprintf(chnuclide,"Moller");   
         arttmin[1]=arttmin[1]/180.*pi;
         arttmax[1]=arttmax[1]/180.*pi;
         artfmin[1]=artfmin[1]/180.*pi;
         artfmax[1]=artfmax[1]/180.*pi;
      }
      else if (strncmp(chn,"Pai",3)==0 || strncmp(chn,"PAI",3)==0 || strncmp(chn,"pai",3)==0){
         sprintf(chnuclide,"E+E- external");
         arttmin[1]=arttmin[1]/180.*pi;
         arttmax[1]=arttmax[1]/180.*pi;
         artfmin[1]=artfmin[1]/180.*pi;
         artfmax[1]=artfmax[1]/180.*pi;
      }
      else{
	 printf("unknown background & source nuclide %s",chn);
         ier=1;
         return;
      }   
   }

   if(ievstart<=0) ievstart=1;
   if(ievstart==1) {rnd->SetSeed(0);}
   if(ievstart!=1){ rnd->SetSeed(0); }
   if(iwrfile==0) sprintf(chfile,"no file");
   if (strncmp(chfile,"no file",7)!=0){
      time_t rawtime;
      struct tm * timeinfo;
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );

      file=fopen(chfile,"w");     
      fprintf(file," DECAY0 generated file: %s \n\n",chfile);
      fprintf(file,"  date and hour         : %s ",asctime (timeinfo)); 
      fprintf(file,"  initial random number : %d \n",irndmst); 
      if (i2bbs==1){
    	 fprintf(file,"  event type: %s           %s\n",chnuclide,chmodebb);
	 fprintf(file,"              level, Elevel:  %s, %f MeV\n",chdspin,fLevelE/1000.);
	 if (toallevents>1.){
	    fprintf(file,"range for sum of energies of emitted e-/e+: %f, %f (MeV)",ebb1,ebb2);
	    allevents=(int)round(nevents*toallevents);
	    fprintf(file,"corresponding number of events in full energy range: %d\n",allevents);   
	 }
      }
      if (i2bbs==2){
         fprintf(file,"  event type: %s           \n",chnuclide);
         printf("  event type: %s           \n",chnuclide);
         if (strncmp(chnuclide,"Artificial",10)==0){
            printf("number of parts in artificial event %d",nartparts);
            for (int i=0;i<nartparts;i++){
                if (strncmp(chart[i],"GP",2)==0){
                   if(nartnpg[i]== 1) sprintf(chn,"Gamma");
	           if(nartnpg[i]== 2) sprintf(chn,"Positron");
 	           if(nartnpg[i]== 3) sprintf(chn,"Electron");
	           if(nartnpg[i]== 4) sprintf(chn,"Neutrino");
	           if(nartnpg[i]== 5) sprintf(chn,"Muon+");  
	           if(nartnpg[i]== 6) sprintf(chn,"Muon-");  
 	           if(nartnpg[i]== 7) sprintf(chn,"Pion0");  
 	           if(nartnpg[i]== 8) sprintf(chn,"Pion+");  
	           if(nartnpg[i]== 9) sprintf(chn,"Pion-");  
 	           if(nartnpg[i]==10) sprintf(chn,"Kaon0 long");
                   if(nartnpg[i]==11) sprintf(chn,"Kaon+ "); 
 	           if(nartnpg[i]==12) sprintf(chn,"Kaon-");  
	           if(nartnpg[i]==13) sprintf(chn,"Neutron");
	           if(nartnpg[i]==14) sprintf(chn,"Proton"); 
	           if(nartnpg[i]==15) sprintf(chn,"Antiproton");               
                   if(nartnpg[i]==16) sprintf(chn,"Kaon0 short");              
                   if(nartnpg[i]==17) sprintf(chn,"Eta");    
	           if(nartnpg[i]==18) sprintf(chn,"Lambda"); 
	           if(nartnpg[i]==19) sprintf(chn,"Sigma+"); 
	           if(nartnpg[i]==20) sprintf(chn,"Sigma0"); 
	           if(nartnpg[i]==21) sprintf(chn,"Sigma-"); 
	           if(nartnpg[i]==22) sprintf(chn,"Xi0");    
	           if(nartnpg[i]==23) sprintf(chn,"Xi- ");   
	           if(nartnpg[i]==24) sprintf(chn,"Omega "); 
                   if(nartnpg[i]==25) sprintf(chn,"Antineutron");
                   if(nartnpg[i]==26) sprintf(chn,"Antilambda");
                   if(nartnpg[i]==27) sprintf(chn,"Antisigma-");
                   if(nartnpg[i]==28) sprintf(chn,"Antisigma0"); 
                   if(nartnpg[i]==29) sprintf(chn,"Antisigma+ ");              
                   if(nartnpg[i]==30) sprintf(chn,"Antixi0");
 	           if(nartnpg[i]==31) sprintf(chn,"Antixi+");
	           if(nartnpg[i]==32) sprintf(chn,"Antiomega+");  
                   if(nartnpg[i]==45) sprintf(chn,"Deuteron ");	 
                   if(nartnpg[i]==46) sprintf(chn,"Tritium");
 	           if(nartnpg[i]==47) sprintf(chn,"Alpha "); 
	           if(nartnpg[i]==48) sprintf(chn,"Geantino");
	           if(nartnpg[i]==49) sprintf(chn,"He3");    
	           if(nartnpg[i]==50) sprintf(chn,"Cerenkov");
               }
               if (strncmp(chart[i],"Be",2)==0)sprintf(chn,"Beta"); 
                  if (strncmp(chart[i],"Pi",2)==0)sprintf(chn,"E+E- internal"); 
                     fprintf(file,"%d %s min and max E= %f %f\n",i,chn,artemin[i],artemax[i]);
                     if (strncmp(chart[i],"Be",2)==0){
                        fprintf(file,"Z of daughter nucleus %d\n",fArtZd[i]);
                     }
                     fprintf(file,"min and max teta %f %f and phi [rad] %f %f= \n",arttmin[i],arttmax[i],artfmin[i],artfmax[i]);  
                     fprintf(file,"time of delay and halflife= %f %f\n",arttdelay[i],artthalf[i]);  
                 }
              }     
              if (strncmp(chnuclide,"Compton",7)==0){
                 fprintf(file,"initial gamma: min and max E %f %f",artemin[1],artemax[1]);
                 fprintf(file,"min and max teta %f %f and phi [rad] %f %f ",arttmin[1],arttmax[1],artfmin[1],artfmax[1]);
              }  
              if (strncmp(chnuclide,"Moller",6)==0){
                 fprintf(file,"initial e-: min and max E %f %f",artemin[1],artemax[1]);
                 fprintf(file,"min and max teta %f %f and phi [rad]  %f %f ",arttmin[1],arttmax[1],artfmin[1],artfmax[1]);
                 fprintf(file,"energy threshold for delta ray %f",artQb[1]); 
              }
              if (strncmp(chnuclide,"E+E- external",13)==0){
                 fprintf(file,"initial gamma: min and max E %f %f",artemin[1],artemax[1]);
                 fprintf(file,"min and max teta %f %f and phi [rad] %f %f ",arttmin[1],arttmax[1],artfmin[1],artfmax[1]);
                 fprintf(file,"atomic number of target %d",fArtZd[1]);
              }
           }
	   fprintf(file," Format of data:\n");
	   fprintf(file,"  for each event    - event's number,    time of event's start, \n");
	   fprintf(file,"                      number of emitted, particles;\n");
	   fprintf(file,"  for each particle - GEANT number of particle, (x,y,z) components of momentum,\n");
	   fprintf(file,"                      time shift from previous time\n");
	   fprintf(file," Time - in sec, momentum - in MeV/c\n");
	   fprintf(file," First event and full number of events: \n %d  %d\n\n",ievstart,nevents);
   } 
   fStartbb=0;
   if (istart==-1 && i2bbs==1) {
      bb();
      istart=1;
   }
   istart=1;
  }
  fNbPart=0;

  if (i2bbs==1) {
	  tevst=0.;
	  bb();

	  if (strncmp(chnuclide,"Ca48",4)==0       && fLevelE!=0.)  Ti48low();
	  else if (strncmp(chnuclide,"Ni58",4)==0  && fLevelE!=0.)  Fe58low();
	  else if (strncmp(chnuclide,"Ge76",4)==0  && fLevelE!=0.)  Se76low();
	  else if (strncmp(chnuclide,"Se74",4)==0  && fLevelE!=0.)  Ge74low();
	  else if (strncmp(chnuclide,"Se82",4)==0  && fLevelE!=0.)  Kr82low();
	  else if (strncmp(chnuclide,"Sr84",4)==0  && fLevelE!=0.) {
             fTclev=0.;
             fThlev=4.35e-12;
             nucltransK(0.882,0.014,6.8e-4,0.);
          } 
	  else if (strncmp(chnuclide,"Zr94",4)==0 && fLevelE!=0.) {
             fTclev=0.;
             fThlev=2.9e-12;
             nucltransK(0.871,0.020,9.5e-4,0.);
          }
          else if (strncmp(chnuclide,"Zr96",4)==0  && fLevelE!=0.) Mo96low(); 
	  else if (strncmp(chnuclide,"Mo92",4)==0  && fLevelE!=0.) Zr92low();
	  else if (strncmp(chnuclide,"Mo100",5)==0 && fLevelE!=0.) Ru100low();
          else if (strncmp(chnuclide,"Ru96",4)==0  && fLevelE!=0.) Mo96low(); 
	  else if (strncmp(chnuclide,"Ru104",5)==0 && fLevelE!=0.) {
             fTclev=0.;
             fThlev=9.9e-12;
             nucltransK(0.556,0.024,4.5e-3,0.);
          }
	  else if (strncmp(chnuclide,"Cd106",5)==0 && fLevelE!=0.) Pd106low();
	  else if (strncmp(chnuclide,"Cd116",5)==0 && fLevelE!=0.) Sn116low();
	  else if (strncmp(chnuclide,"Sn112",5)==0 && fLevelE!=0.) Cd112low();
	  else if (strncmp(chnuclide,"Sn124",5)==0 && fLevelE!=0.) Te124low();
          else if (strncmp(chnuclide,"Te120",5)==0 && fLevelE==1171) {
             fTclev=0.;
             fThlev=0.642e-12;
             nucltransK(1.171,0.029,8.0e-4,0.);
          }
          else if (strncmp(chnuclide,"Te128",5)==0 && fLevelE==443) {
             fTclev=0.;
             fThlev=23.8e-12;
             nucltransK(0.443,0.035,8.0e-3,0.);
          }
	  else if (strncmp(chnuclide,"Te130",5)==0 && fLevelE!=0.) Xe130low();
	else if ( (strncmp(chnuclide,"Xe136",5)==0 || strncmp(chnuclide,"Ce136",5)==0) 
                    && fLevelE!=0.) Ba136low();
          else if (strncmp(chnuclide,"Nd148",5)==0 && fLevelE!=0.) Sm148low();
	  else if (strncmp(chnuclide,"Nd150",5)==0 && fLevelE!=0.) Sm150low();
          else if (strncmp(chnuclide,"W186",4)==0  && fLevelE==137) {
             fTclev=0.;
             fThlev=875.e-12;
             // KLM ratios in accordance with BrIcc calculation
             nucltransKLM(0.137,0.074,0.44,0.012,0.64,0.003,0.20,0.);
          }
          else if (strncmp(chnuclide,"Bi214",5)==0){
             fTclev=0.;
             int fNbPart0=fNbPart;
             At214(0.);
             fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
          }
          else if (strncmp(chnuclide,"Pb214",5)==0){          
             fTclev=0.;
             int fNbPart0=fNbPart;
	     Po214(0.);
  	     fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
	  }
          else if (strncmp(chnuclide,"Po218",5)==0){
             fTclev=0.;
             int fNbPart0=fNbPart;              
             Rn218(0.);
             fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
             fNbPart0=fNbPart;
             Po214(0.);
             fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc; 
	  }
          else if (strncmp(chnuclide,"Rn222",5)==0){
             fTclev=0.;
             int fNbPart0=fNbPart;
             Ra222(0.);
             fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
             fNbPart0=fNbPart;
             Rn218(0.);
             fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
             fNbPart0=fNbPart;
             Po214(0.);
             fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
          }
  }
  if (i2bbs==2) {
     if (strncmp(chnuclide,"Ac228",5)==0){
        Ac228(0.);
        tevst=fTdnuc;
     } 
     else if (strncmp(chnuclide,"Ar39",4)==0){
        // Model for Ar39 decay (Nucl. Phys. A 633(1998)1).
        fThnuc=8.488762e9;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        beta1fu(0.565,tcnuc,0.,0.,0.,0.,0.);
     } 
     else if (strncmp(chnuclide,"Ar42",4)==0){
        // Model for Ar42 decay (Nucl. Data Sheets 92(2001)1).
        fThnuc=1.0382166e9;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        beta1fu(0.599,tcnuc,0.,0.,0.,0.,0.);
     }
     else if (strncmp(chnuclide,"As79",4)==0){
        As79(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Bi207",5)==0){
        Bi207(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Bi208",5)==0){
        // Scheme of Bi208 decay (NDS 47(1986)797 + ToI-1978).
        // VIT, 17.12.1995; 10.05.2005
        float tcnuc=0.0;
        fThnuc=1.161288E+13;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        fTclev=0.; fThlev=0.;
        float pdecay=100.*GetRandom();
        if (pdecay<=43.6) fEgamma=0.088;        //EC-K 43.6%
        else if (pdecay>43.6 && pdecay<=83.8) 
                          fEgamma=0.016;       //EC-L 40.2%  
        else fEgamma=0.004;   //EC-M 16.2%
        particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
        fThlev=32.e-12;
        nucltransK(2.615,0.088,8.5e-3,4.3e-4);
     }
     else if (strncmp(chnuclide,"Bi210",5)==0){
        Bi210(0.);
        tevst=fTdnuc;
     } 
     else if (strncmp(chnuclide,"Bi212",5)==0){
        Bi212(0.);
        tevst=fTdnuc;
        int fNbPart0=fNbPart;
        if (npgeant[1]!=47){ //decay of Po212
           double thnuc=0.299e-6;
           float tcnuc=0.;
           fTdnuc=tcnuc-thnuc/log(2.)*log(GetRandom());
           particle(47,8.785,8.785,0.,pi,0.,twopi,0,0);
           fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
        }
     }
     else if (strncmp(chnuclide,"Bi214",5)==0){
        Bi214(0.);
        tevst=fTdnuc;
        int fNbPart0=fNbPart;
        if (npgeant[1]!=47){
           Po214(0.);
           fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
        }
     }
     else if (strncmp(chnuclide,"C14",3)==0){
        // Scheme of C14 beta decay, NPA 523(1991)1 and ToI'1998
        // VIT, 5.11.2006.
        // experimental corrections to the allowed beta shape from
        // F.E.Wietfeldt et al., PRC 52(1995)1028
        // cf(e)=(1+c1/w+c2*w+c3*w**2+c4*w**3), w=e/emass+1 
        fThnuc=1.798734e11;
        float tcnuc=0.; 
	fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        beta1f(0.156,tcnuc,0.,0.,0.19418,0.,0.); 
     }
     else if (strncmp(chnuclide,"Ca48",4)==0){
        // Scheme of Ca48 beta decay. 
        // It is supposed that decay goes to excited 5+ level of Sc48 (E_exc=131 keV)
        // with T1/2=1.1e21 y calculated in M.Aunola et al., Europhys. Lett. 46(1999)577 
        // (transition to g.s. is suppressed by factor of ~1e-10; to 4+ excited level 
        // with E_exc=252 keV - by factor of ~1e-5).
        // VIT, 07.05.1998; update of 13.08.2007.
        fThnuc=3.47e28; 
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        fZdtr=21; 
        beta(0.151,0.,0.);
        fTclev=0.;
	fThlev=0.;
        nucltransK(0.131,0.004,8.1e-3,0.);
        int fNbPart0=fNbPart;
        Sc48(0.);
        fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
     }
     else if (strncmp(chnuclide,"Cd113",5)==0){
        // Scheme of Cd113 beta decay. 
        // Half-life and coefficients in the shape correction factor 
        // are taken from: F.A.Danevich et al., Phys. At. Nuclei 59(1996)1.
        //  Q_beta=0.320 MeV, G.Audi et al., Nucl. Phys. A 729(2003)337.
        // VIT, 31.03.2006.
        fThnuc=2.42987e23;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        float c1= 1.01*7.;
        float c2= 1.48*7.;
        float c3= 0.684;
        beta2f(0.320,0.,0.,3,c1,c2,c3,0.);
     }
     else if (strncmp(chnuclide,"Co60",4)==0){
        Co60(0.); 
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Cs136",5)==0){
        Cs136(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Cs137",5)==0){
        // Model for scheme of Cs137 decay (Nucl. Data Sheets 72(1994)355)
        // (really it is model for (Cs137 + Ba137m)-decay but not the model of
        // Cs137 decay alone). 
        fThnuc=0.9489110E+09;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        fTclev=0.;
        float pbeta=100.*GetRandom();
        if (pbeta<=94.4){
           beta1fu(0.514,0.,0.,0.,0.,0.,0.);
           fThlev=153.12;
           nucltransKL(0.662,0.037,9.0e-2,0.006,1.6e-2,0.);
        }
        else{
           beta1f(1.176,0.,0.,0.,-0.6060315,0.0921520,0.);
        }
     }
     else if (strncmp(chnuclide,"Eu147",5)==0){
        Eu147(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Eu152",5)==0){
        Eu152(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Eu154",5)==0){
        Eu154(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Gd146",5)==0){
        Gd146(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Hf182",5)==0){
        // Scheme of 182Hf decay("Table of Isotopes", 7th ed., 1978) 
        // VIT, 5.03.1996.
        fThnuc=2.84011e14;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());   
        tevst=fTdnuc;
        fTclev=0.;
        beta(0.160,0.,0.);
        fThlev=1.2e-9;
        float p=100.*GetRandom();
        if (p<=91.77){
           nucltransK(0.270,0.067,3.1e-1,0.);
        }
        else if (p<=92.00){
           nucltransK(0.173,0.067,9.5e-2,0.);
           fThlev=0.;
           nucltransK(0.098,0.067,5.0e-0,0.);
        } 
        else{
           nucltransK(0.156,0.067,1.5e-1,0.);
           fThlev=0.;
           nucltransK(0.114,0.067,4.5e-0,0.);
        }

     }
     else if (strncmp(chnuclide,"I126",4)==0){
        I126(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"I133",4)==0){
        I133(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"I134",4)==0){
        I134(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"I135",4)==0){
        I135(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"K40",3)==0){
        K40(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"K42",3)==0){
        K42(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Kr81",4)==0){
        // Scheme of Kr81 decay (NDS 79(1996)447 and ENSDF at NNDC 
        // site on 9.12.2007). 
        // VIT, 9.12.2007.
        fThnuc=7.226493e12;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        fTclev=0.;
        float pdecay,pklm;
        pdecay=100.*GetRandom();
        if (pdecay<=0.30){
           //capture from only K shell is supposed
           particle(1,0.013,0.013,0.,pi,0.,twopi,0,0);
           fThlev=9.7e-12;
           nucltransKLM(0.276,0.013,7.3e-3,0.002,7.8e-4,0.000,2.6e-4,0.); 
        } 
        else{
           pklm=100.*GetRandom();
           if (pklm<=84.73) particle(1,0.013,0.013,0.,pi,0.,twopi,0,0);
           if (pklm>84.73 && pklm<=97.44) particle(1,0.002,0.002,0.,pi,0.,twopi,0,0);
           if (pklm>97.44) particle(1,0.000,0.000,0.,pi,0.,twopi,0,0); 
        }
     }
     else if (strncmp(chnuclide,"Kr85",4)==0){
        // Scheme of Kr85 decay (NDS 62(1991)271 and ENSDF at NNDC 
        // site on 9.12.2007).
        // VIT, 9.12.2007.
        fThnuc=3.394243e8;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        fTclev=0.;
        float pdecay,p;
        pdecay=100.*GetRandom();
        if (pdecay<=0.434){
           beta(0.173,0.,0.);
           fThlev=1.015e-6;
           p=100.*GetRandom();
           if (p<=99.99947){
              nucltransKLM(0.514,0.014,6.3e-3,0.002,7.1e-4,0.000,2.3e-4,0.); 
           }
           else{
              nucltransKLM(0.363,0.014,2.9e-2,0.002,3.9e-3,0.000,1.3e-3,0.);
              fThlev=0.71e-9;
              nucltransKLM(0.151,0.014,4.3e-2,0.002,4.8e-3,0.000,3.4e-4,0.);
           }
        }
        else beta2f(0.687,0.,0.,1,1.,0.,0.,0.);// 1st forbidden unique beta decay
     }
     else if (strncmp(chnuclide,"Mn54",4)==0){
        // Scheme of Mn54 decay ("Table of Isotopes",8th ed.,1996 
        // + NDS 50(1987)255).
        // Accuracy in description of: 
        // decay branches - 0.001%, deexcitation process - 0.001%.
        // VIT, 16.04.1998.
        // VIT, 1.04.2007, updated to NDS 107(2006)1393.
        fThnuc=2.696717e+7;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        fTclev=0.;
        particle(1,0.006,0.006,0.,pi,0.,twopi,0.,0.);//100% EC to Cr54
        fThlev=7.9e-12;
        nucltransK(0.835,0.006,2.4e-4,0.);
     }
     else if (strncmp(chnuclide,"Na22",4)==0){
        // Scheme of Na22 decay ("Table of Isotopes", 7th ed., 1978).
        // Accuracy in description of 
        // decay branches - 0.001%, deexcitation process - 0.001%. 
        // VIT, 26.07.1993, 22.10.1995.
        // VIT, 12.11.2006 (updated to NDS 106(2005)1 
        // and change to beta spectra with experimental corrections).
        // VIT, 26.08.2007 (corrected beta shapes).
        fThnuc=8.2132717e+7;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        fTclev=0.;
        float pdecay=100.*GetRandom();
        if (pdecay<=99.944){
           fThlev=3.63e-12;
           if (pdecay<=9.618) particle(1,0.001,0.001,0.,pi,0.,twopi,0,0);  
           else   beta1f(0.545,0.,0.,1.e-3,0.,0.,0.); 
           nucltransK(1.275,0.001,6.8e-6,2.1e-5);
        }
        else{
           beta2f(1.820,0.,0.,2,3.333333,1.,0.,0.);
        }
     }
     else if (strncmp(chnuclide,"P32",3)==0){
        // Scheme of P32 beta decay, ToI'1998 and ENSDF'2004. 
        // experimental corrections to the allowed beta shape from
        // H.Daniel, RMP 40(1968)659 and M.J.Canty et al., NP 85(1966)317
        // VIT, 5.11.2006.
        fThnuc=1.2323232e6;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        beta1f(1.710,0.,0.,0.,0.003,0.,0.);
     }
     else if (strncmp(chnuclide,"Pa234",5)==0){
        Pa234(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Pb210",5)==0){
        // Scheme of Pb210 decay in accordance with NDS 99(2003)649 and 
        // ENSDF at the NNDC site on 6.08.2007.
        // VIT, 6.08.2007.
        fThnuc=7.0056e8;
        fThlev=0;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        float pbeta=100.*GetRandom();  // beta decay to U234
        fTclev=0.;
        if (pbeta<=84.){
           beta(0.0170,0.,0.);
           nucltransKLM(0.0465,0.016,14.2,0.004,3.36,0.001,1.14,0.);
        }
        else beta(0.0635,0.,0.); 
     }
     else if (strncmp(chnuclide,"Pb211",5)==0){
        Pb211(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Pb212",5)==0){
        Pb212(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Pb214",5)==0){
        Pb214(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Ra228",5)==0){
        Ra228(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Rb87",4)==0){
        // Scheme of Rb87 decay in accordance with NDS 95(2002)543 
        // and ENSDF at NNDC site on 6.08.2007.
        // VIT, 6.08.2007.
        fThnuc=1.518e18;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom()); 
        tevst=fTdnuc;
        beta2f(0.283,0.,0.,2,27.72,90.91,0.,0.);
     }
     else if (strncmp(chnuclide,"Rh106",5)==0){
        Rh106(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Sb125",5)==0){
        Sb125(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Sb126",5)==0){
        Sb126(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Sb133",5)==0){
        Sb133(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Sr90",4)==0){  
        // Scheme of Sr90 decay ("Table of Isotopes", 7th ed., 1978).
        // Slight update to NDS 82(1997)379.
        // VIT, 9.08.1993, 22.10.1995, 26.10.2006.
        // Change from the allowed shape to the 1st forbidden unique 
        // with empirical correction from: 
        // H.H.Hansen, Appl. Rad. Isot. 34(1983)1241
        fThnuc=0.9085184e+9;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        beta1fu(0.546,0.,0.,0.,-0.032,0.,0.);
     }
     else if (strncmp(chnuclide,"Ta182",5)==0){
        Ta182(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Te133m",6)==0){
        Te133m(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Te133",5)==0){
        Te133(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Te134",5)==0){
        Te134(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Th234",5)==0){
        Th234(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Tl207",5)==0){
        // Scheme of Tl207 decay ("Table of Isotopes", 7th ed., 1978).
        // VIT, 14.08.1992, 22.10.1995; 
        // VIT, 30.10.2006 (update to NDS 70(1993)315 and correction 
        // to the beta shape). 
        fThnuc=286.2;
        fThlev=0;
        fTclev=0.;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        float pbeta, p;
        pbeta=100.*GetRandom();
        if (pbeta<= 0.268){
           beta(0.529,0.,0.);
           fThlev=0.115e-12;
           p=100*GetRandom();
           if (p<= 99.29){
              nucltransK(0.898,0.088,2.5e-2,0.);
           }
           else{
              nucltransK(0.328,0.088,3.5e-1,0.);
              fThlev=130.5e-12;
              nucltransK(0.570,0.088,2.2e-2,0.);
           }
        }
        else{
           // change to forbidden spectrum with experimental correction 
           // from J.M.Trischuk et al., NPA 90(1967)33 and H.Daniel, 
           // RMP 40(1968)659
           beta1f(1.427,0.,0.,0.,0.024,0.,0.);
        }
     }
     else if (strncmp(chnuclide,"Tl208",5)==0){
        Tl208(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Xe129m",6)==0){
        // Scheme of Xe129m decay (NDS 77(1996)631 and 
        // ENSDF at NNDC site on 13.11.2007).
        // VIT, 13.11.2007.
        fThnuc=767232.;
        fThlev=0.;
        fTclev=0.;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        nucltransKLM(0.197,0.035,13.94,0.005,5.34,0.001,1.52,0.);
        fThlev=0.97e-9;
        nucltransKLM(0.040,0.035,10.49,0.005,1.43,0.001,0.39,0.);
     }
     else if (strncmp(chnuclide,"Xe131m",6)==0){
        // Scheme of Xe129m decay (NDS 107(2006)2715 and 
        // ENSDF at NNDC site on 13.11.2007).
        // VIT, 13.11.2007.
        fThnuc=1.022976e+6;
        fThlev=0.;
        fTclev=0.;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        nucltransKLM(0.164,0.035,31.60,0.005,14.75,0.001,4.15,0.);
     }
     else if (strncmp(chnuclide,"Xe133",5)==0){
        Xe133(0.);
        tevst=fTdnuc; 
     }
     else if (strncmp(chnuclide,"Xe135",5)==0){
        Xe135(0.);
        tevst=fTdnuc; 
     }
     else if (strncmp(chnuclide,"Y88",3)==0){
        Y88(0.);
        tevst=fTdnuc; 
     }
     else if (strncmp(chnuclide,"Y90",3)==0){
        // Scheme of Y90 decay ("Table of Isotopes", 7th ed., 1978).
        // Accuracy in description of: decay branches       - 0.001%,
        //                           : deexcitation process - 0.001%.
        // VIT, 9.08.1993, 22.10.1995, 26.10.2006; 27.02.2007.
        fThnuc=230400.;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        fTclev=0.;
        fThlev=0.;
        float pbeta, p;
        pbeta=100.*GetRandom(); 
        if (pbeta<=0.0115){
           beta1fu(0.519,0.,0.,0.,0.,0.,0.); 
           fThlev=61.3e-9;
           p=100.*GetRandom();
           if (p<=27.7) {
              pair(0.739);
           }
           else{
              particle(3,1.743,1.743,0.,pi,0.,twopi,fTclev,fThlev);//electron
              particle(1,0.018,0.018,0.,pi,0.,twopi,0,0);//gamma
           }
        }
        else  beta1fu(2.280,0.,0.,0.,-0.0078,0.,0.);
     }
     else if (strncmp(chnuclide,"Zn65",4)==0){
        Zn65(0.);
        tevst=fTdnuc;
     }
     else if (strncmp(chnuclide,"Zr96",4)==0){
        // Scheme of Zr96 beta decay. It is supposed that decay goes 
        // to the first excited level (5+) of Nb96 (E_exc=44 keV) 
        // with T1/2=2.4e20 y in accordance with: 
        // H.Heiskanen et al.,JPG 34(2007)837. 
        // Transition to the ground state (6+) is suppressed by factor 
        // of ~1e-9, to the excited 4+ (E_exc=146 keV) - by 1e-2.
        // VIT, 07.05.1998.
        // VIT, 18.08.2007:update to NDS 68(1993)165 and ENSDF at NNDC site.
        fThnuc=7.574e27;
        float tcnuc=0.;
        fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
        tevst=fTdnuc;
        fTclev=0.;
        fThlev=0.;
        fZdtr=41;
        beta(0.117,0.,0.);
        nucltransK(0.044,0.019,2.4e0,0.);
        int fNbPart0=fNbPart;
        Nb96(0.);
        fPtime[fNbPart0+1]=fPtime[fNbPart0+1]+fTdnuc;
     }
  } 
  if (strncmp(chfile,"no file",7)!=0){
     fprintf(file,"  % 2d % 2.2lf % 2d\n",icurrent+1,tevst,fNbPart);
     for (int j=1;j<fNbPart+1;j++){
	 fprintf(file,"% 3d % 3.4E % 3.4E % 3.4E % 3.4E\n",
	         npgeant[j],fPmoment[0][j],fPmoment[1][j],fPmoment[2][j],fPtime[j]);
     }

  } 
  if(icurrent!=nevents) return;
  printf("toallevents= %f\n",toallevents);
}

void Decay::bb()
{
  // function for sampling the energies and angles of electrons in various 
  // modes of double beta decay without Primakoff-Rosen approximation.
  //   Edlevel   - energy of level of daughter nucleus to which transition occured;
  //   fStartbb  - must be =0 for first call of bb for a given mode;
  //   ebb1,ebb2 - for modes with continuous distribution of sum of
  //               e-/e+ energies (4,5,6,8,10 and 13), left and right
  //               limits for E sum in which events will be generated;
  
  float Edlevel=fLevelE/1000.;

  if (fZdbb>=0) e0=fQbb-Edlevel;
  if (fZdbb<0.) e0=fQbb-Edlevel-4.*emass;

  if (fModebb==9) {
     e0=fQbb-Edlevel-fEK-2.*emass;
     particle(2,e0,e0,0.,pi,0.,twopi,0.,0.);
     particle(1,fEK,fEK,0.,pi,0.,twopi,0.,0.);
     return;
  } 
  else if (fModebb==10) e0=fQbb-Edlevel-fEK-2.*emass;
  else if (fModebb==11) {
     e0=fQbb-Edlevel-2.*fEK;
     particle(1,e0,e0,0.,pi,0.,twopi,0.,0.);
     particle(1,fEK,fEK,0.,pi,0.,twopi,0.,0.);
     particle(1,fEK,fEK,0.,pi,0.,twopi,0.,0.);
     return;
  }
  else if (fModebb==12){
     e0=fQbb-Edlevel-2.*fEK; 
     particle(1,fEK,fEK,0.,pi,0.,twopi,0.,0.);
     particle(1,fEK,fEK,0.,pi,0.,twopi,0.,0.);
     return;
  }
  if (fStartbb==0){ 
     spmax=-1.;
     //calculate the theoretical energy spectrum of first particle 
     //with step of 1 keV and find its maximum
     printf("wait, please: calculation of theoretical spectrum\n");
     if(ebb1<0)  ebb1=0;
     if(ebb2>e0) ebb2=e0;
     int maxi=(int) (e0*1000);
     for (int i=1;i<maxi+1;i++){
 	 e1=i/1000.;
	 spthe1[i]=0.;
 	 if ((fModebb>=1 && fModebb<=3) ||fModebb==7 || fModebb==10)  
            spthe1[i]=fe1_mod(e1);

	 elow=max(1e-4,ebb1-e1+1e-4);
	 ehigh=max(1e-4,ebb2-e1+1e-4);
	 if (((fModebb>=4 && fModebb<=6) ||fModebb==8 || (fModebb>=13 && fModebb<=16))&& e1<e0) {
	    TF1  f1 ("f1",this,0,e0,4,"Decay");
            ROOT::Math::WrappedTF1 wf1(f1);
	    ROOT::Math::GaussLegendreIntegrator ig;
	    ig.SetFunction(wf1);
	    ig.SetRelTolerance(0.001);
	    ig.SetNumberPoints(40);
	    spthe1[i]=  ig.Integral(elow,ehigh);
	 }
	if(spthe1[i]>spmax) spmax=spthe1[i]; 
     }

     toallevents=1.;
     printf("toallevents=%1.2lf\n",toallevents);

   }
   fStartbb=1;
   float re2s, re2f;
   int k,ke2s,ke2f;
   float fe2=0.;
   do{
      if (fModebb==10) e1=ebb1+(ebb2-ebb1)*GetRandom();
      if (fModebb!=10) e1=ebb2*GetRandom();

      k = (int)(round(e1*1000.));
      if (k<1) k=1;
   } while ( (spmax*GetRandom()) > spthe1[k]);   

   if ((fModebb>=1 && fModebb<=3) || fModebb==7) e2=e0-e1; 
   else if ((fModebb>=4 && fModebb<=6) ||fModebb==8 || (fModebb>=13 && fModebb<=16)){
      float zero=0.0;
      re2s=max(zero,ebb1-e1);
      re2f=ebb2-e1;
      float f2max=-1.;
      ke2s=(int)(max(1.0,re2s*1000.));
      ke2f=(int)(round(re2f*1000.));

      for (int ke2=ke2s; ke2<ke2f+1; ke2++){
          e2=ke2/1000.;
	  spthe2[ke2]=fe2_mod(e2);
	  if (spthe2[ke2]>f2max) f2max=spthe2[ke2];
      }
      do{
          e2=re2s+(re2f-re2s)*GetRandom();
	  if ((fModebb>=4&& fModebb<=6) || fModebb==8 || (fModebb>=13 && fModebb<=16) ) 
	  fe2=fe2_mod(e2);
      }while (f2max*GetRandom()>fe2); 
    }
    else if( fModebb==10){
      particle(2,e1,e1,0.,pi,0.,twopi,0.,0.);
      particle(1,fEK,fEK,0.,pi,0.,twopi,0.,0.);
      return;
    }

    p1=sqrt(e1*(e1+2.*emass));
    p2=sqrt(e2*(e2+2.*emass));
    float b1=p1/(e1+emass);
    float b2=p2/(e2+emass);
    float a,b,c;
    a=1.;
    b=-b1*b2;
    c=0.;
    float w1, w2;
    if (fModebb==2) b=b1*b2;
    else if (fModebb==3) {
       w1=e1+emass;
       w2=e2+emass;
       a=3.*(w1*w2+emass*emass)*(p1*p1+p2*p2);
       b=-p1*p2*(pow(w1+w2,2)+4.*(w1*w2+emass*emass));
       c=2.*p1*p1*p2*p2;
    }
    else if (fModebb==7) {
       w1=e1+emass;
       w2=e2+emass;
       a=5.*(w1*w2+emass*emass)*(p1*p1+p2*p2)-p1*p1*p2*p2;
       b=-p1*p2*(10.*(w1*w2+emass*emass)+p1*p1+p2*p2);
       c=3.*p1*p1*p2*p2;
    }
    else if (fModebb==8|| fModebb==16) b=b1*b2/3.;
    else if (fModebb==15) {
       	a=9.*pow(e0-e1-e2,2)+21.*pow(e2-e1,2);
	b=-b1*b2*(9.*pow(e0-e1-e2,2)-7.*pow(e2-e1,2));
    }
    float ctet,stet1,stet2,ctet1,ctet2;
    float romaxt=a+abs(b)+c;
    do{
	phi1=twopi*GetRandom();
	ctet1=1.-2.*GetRandom();
	stet1=sqrt(1.-ctet1*ctet1);
	phi2=twopi*GetRandom();
	ctet2=1.-2.*GetRandom();
	stet2=sqrt(1.-ctet2*ctet2);
	ctet=ctet1*ctet2+stet1*stet2*cos(phi1-phi2);
    } while(romaxt*GetRandom()> a+b*ctet+c*ctet*ctet);

    fNbPart=fNbPart+1;
    if(fZdbb>=0.) npgeant[fNbPart]=3;
    if(fZdbb<0.) npgeant[fNbPart]=2;
    fPmoment[0][fNbPart]=p1*stet1*cos(phi1);
    fPmoment[1][fNbPart]=p1*stet1*sin(phi1);
    fPmoment[2][fNbPart]=p1*ctet1;
    fPtime[fNbPart]=0.;
    fNbPart=fNbPart+1;
    if(fZdbb>=0.) npgeant[fNbPart]=3;
    if(fZdbb<0.) npgeant[fNbPart]=2;
    fPmoment[0][fNbPart]=p2*stet2*cos(phi2);
    fPmoment[1][fNbPart]=p2*stet2*sin(phi2);
    fPmoment[2][fNbPart]=p2*ctet2;
    fPtime[fNbPart]=0.;
    return;
}

//-----------------------------------------------
//-----------------------------------------------
//-----------------------------------------------
void Decay::Ti48low()
{
  // Describes the deexcitation process in Ti48 nucleus  after 2b-decay of Ca48 
  // to ground and excited 0+ and 2+ levels of Ti48 ("Table of Isotopes", 7th ed., 1978).
  fTclev=0.;
  if (fLevelE==2421.){
     fThlev=24.e-15;
     float p=100.*GetRandom();
     if (p<8) nucltransK(2.421,0.005,1.5e-5,5.0e-4);
     else {
	      nucltransK(1.437,0.005,3.1e-5,1.8e-4);
	      fThlev=4.3e-12;
	      nucltransK(0.984,0.005,1.2e-4,0.);
     }
     return;
  }
  if (fLevelE==984.){
     fThlev=4.3e-12;
     nucltransK(0.984,0.005,1.2e-4,0.);
     return;
  } 
  return;
}

//-----------------------------------------------
void Decay::Fe58low()
{
  // describes the deexcitation process in Fe58 nucleus after 2b-decay of Ni58 
  // to ground and excited 0+ and 2+ levels of Fe58 ("Table of Isotopes", 7th ed., 1978).
  fTclev=0.;
  float p;
  if (fLevelE==1675.){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<43.) nucltransK(1.675,0.007,1.0e-4,0.6e-4);
     else{
	nucltransK(0.864,0.007,3.0e-4,0.);
	fThlev=9.0e-12;
	nucltransK(0.811,0.007,5.0e-4,0.);
     }
     return;
   }
   if (fLevelE==811.){
      fThlev=9.0e-12;
      nucltransK(0.811,0.007,5.0e-4,0.);
      return;
   }
   return;
}

//-----------------------------------------------
void Decay::Se76low()
{
   // describes the deexcitation process in Se76 nucleus after 2b-decay of Ge76 
   // to ground and excited 0+ and 2+ levels of Se76
   int npg563=0;
   int npg559=0;
   float p;
   fEbindeK=0.013;
   float cg=1.;
   float cK=2.0e-3;
   fTclev=0.;
   bool next=false;
   if (fLevelE==1216){
      fThlev=3.4e-12;
      p=100.*GetRandom();
      if (p<36.) nucltransK(1.216,0.013,4.3e-4,0.1e-4);
      else{
         nucltransK(0.657,0.013,2.1e-3,0.);
         next=true;
      }
      if(!next) return;
   }
   if (fLevelE==1122){
      fThlev=11.e-12;
      fEgamma=0.563;
      p=(cg+cK)*GetRandom();
      if (p<=cg){ 
         particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
	 npg563=fNbPart;
      }
      else{ 
         particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
	 particle(1,fEbindeK,fEbindeK,0.,pi,0.,twopi,0.,0.);
      }
         next=true;
   }
      
   
   if (fLevelE==559 || next){ 
      fThlev=12.3e-12;
      fEgamma=0.559;
      p=(cg+cK)*GetRandom();
      if (p<cg){
         particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
         npg559=fNbPart;
      }
      else{ 
         particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
         particle(1,fEbindeK,fEbindeK,0.,pi,0.,twopi,0.,0.);
      }
      // Angular correlation between gammas 559 and 563 keV, L.Pandola + VIT

      if (npg559!=0 && npg563!=0){  
         float p559=sqrt(pow(fPmoment[0][npg559],2)+pow(fPmoment[1][npg559],2)+ 
                          pow(fPmoment[2][npg559],2));
         float p563=sqrt(pow(fPmoment[0][npg563],2)+pow(fPmoment[1][npg563],2)+ 
                          pow(fPmoment[2][npg563],2));
         // Coefficients in formula 1+a2*ctet**2+a4*ctet**4 are from: R.D.Evans, 
         //"The Atomic Nucleus", Krieger Publ. Comp., 1985, p. 240, 0(2)2(2)0 cascade.
         float ctet,stet1,stet2,ctet1,ctet2;
         float a2=-3.0;
         float a4=4.0;
         do{
           phi1=twopi*GetRandom();
           ctet1=1.-2.*GetRandom();
           stet1=sqrt(1.-ctet1*ctet1);
           phi2=twopi*GetRandom();
           ctet2=1.-2.*GetRandom();
  	   stet2=sqrt(1.-ctet2*ctet2);
  	   ctet=ctet1*ctet2+stet1*stet2*cos(phi1-phi2);
         }while(GetRandom()*(1.+abs(a2)+abs(a4)) > 1.+a2*ctet*ctet+a4*pow(ctet,4)); 
      
         fPmoment[0][npg559]=p559*stet1*cos(phi1);
         fPmoment[1][npg559]=p559*stet1*sin(phi1);
         fPmoment[2][npg559]=p559*ctet1;
         fPmoment[0][npg563]=p563*stet2*cos(phi2);
         fPmoment[1][npg563]=p563*stet2*sin(phi2);
         fPmoment[2][npg563]=p563*ctet2;

         return;
      }
   }
   return;
}

//-----------------------------------------------
void Decay::Ge74low()
{
  // Subroutine describes the deexcitation process in Ge74 nucleus after 2b-decay 
  // of Se74 to ground and excited 0+ and 2+ levels of Ge74

  fTclev=0.;
  bool next=false;

  if (fLevelE==1204){
     fThlev=6.0e-12;
     float p=100.*GetRandom();
     if (p<34) nucltransK(1.204,0.011,1.9e-4,0.1e-4);
     else {
        nucltransK(0.608,0.011,1.1e-3,0.);
        next=true;     
     }
     if(!next) return;
  }
  if (fLevelE==596 || next){
     fThlev=12.0e-12;
     nucltransK(0.596,0.011,1.1e-3,0.);
     return;
  }
  return;
}

//-----------------------------------------------
void Decay::Kr82low()
{
  //Subroutine describes the deexcitation process in Kr82 nucleus after 2b-decay 
  //of Se82 to ground and excited 0+ and 2+ levels of Kr82
  fTclev=0.;
  bool next=false;

  if (fLevelE==1475){
     fThlev=0.;
     float p=100.*GetRandom();
     if (p<=36.7) nucltransK(1.475,0.014,2.0e-4,0.5e-4); 
     else{
        nucltransK(0.698,0.014,1.3e-3,0.);
        next=true;
    }
    if(!next) return;
  }
  if (fLevelE==776 || next){
     fThlev=5.e-12;
     nucltransK(0.776,0.014,9.3e-4,0.);
     return; 
  }
  return;
}

//-----------------------------------------------
void Decay::Mo96low()
{
  //Subroutine describes the deexcitation process in Mo96 nucleus after 2b- decay 
  //of Zr96 and 2b+/eb+/2e decay of Ru96 to ground and excited 0+ and 2+ levels 
  //of Mo96 ("Table of Isotopes", 7th ed., 1978).  
  //VIT, 20.05.2009: four levels (2096, 2426, 2623 and 2700 keV) are added 
  //and information is updated in accordance with NNDC site on 19.05.2009 
  //and  NDS 109(2008)2501.
  fTclev=0.;
  bool next778=false;
  bool next1148=false;
  bool next1498=false;
  bool next1626=false;

   if (fLevelE==2713){
      fThlev=0.;
      nucltransK(0.272,0.020,5.9e-1,0.);
      fThlev=0.208e-12;
      nucltransK(0.813,0.020,1.3e-3,0.);
      fThlev=1.2e-12;
      nucltransK(0.850,0.020,1.1e-3,0.);
      next778=true;
   }
   if (fLevelE==2700){
      fThlev=0.103e-12;
      bool go=false;
      float p=100.*GetRandom();
      if (p<=3.04){
         nucltransK(0.160,0.020,6.9e-2,0.);
         fThlev=0.069e-12;
         float p=100.*GetRandom();
         if (p<=8.17){
            nucltransK(0.915,0.020,9.7e-4,0.);
            next1626=go=true;

         }
         else if (p<=28.31){
            nucltransK(1.043,0.020,7.3e-4,0.);
            next1498=go=true; 
         }
         else{
            nucltransK(1.762,0.020,2.5e-4,1.7e-4);
            next778=go=true;
         }

      }
      else if (p<=12.12){
         nucltransK(1.074,0.020,6.7e-4,0.);
         next1626=go=true;
      }
      else if (p<=53.40){
         nucltransK(1.202,0.020,5.5e-4,0.);
         next1498=go=true;
      }
      else if (p<=86.67){
         nucltransK(1.922,0.020,2.1e-4,2.5e-4);
         next778=go=true;
      }
      else{
         nucltransK(2.701,0.020,1.1e-4,6.4e-4);
         return;
      }
     if(!go) return;
   }
   if (fLevelE==2623){
      fThlev=0.6e-12;
      nucltransK(1.844,0.020,2.2e-4,2.3e-4);
      next778=true;
   }
   if (fLevelE==2426){
      fThlev=0.19e-12;
      bool go=false;
      float p=100.*GetRandom();
      if (p<=2.50){
         nucltransK(0.448,0.020,6.4e-3,0.);
         fThlev=2.29e-12;
         float p=100.*GetRandom();
         if (p<=0.18){
            nucltransK(0.109,0.020,2.0e-1,0.);
            fThlev=6.4e-12;
            float p=100.*GetRandom();
            if (p<=7.52){
               nucltransK(0.241,0.020,2.4e-2,0.);
               fThlev=1.2e-12;
               nucltransK(0.850,0.020,1.1e-3,0.);
               next778=go=true;
            }
            else if (p<=12.22){
               nucltransK(0.372,0.020,1.2e-2,0.);
               next1498=go=true;
            }
         }
         else if (p<=5.96){
            nucltransK(0.350,0.020,1.2e-2,0.); 
            fThlev=1.2e-12;
            nucltransK(0.850,0.020,1.1e-3,0.);
            next778=go=true;
         }
         else if (p<=9.27){
            nucltransK(0.353,0.020,1.2e-2,0.);
            next1626=go=true;
         }
         else if (p<=30.34){
            nucltransK(0.481,0.020,4.3e-3,0.);
            next1498=go=true;
         }
         else{
            nucltransK(1.200,0.020,5.4e-4,6.1e-6);
            next778=go=true;
         }
     }
     else if (p<=38.16){
         nucltransK(0.800,0.020,1.3e-3,0.);
         next1626=go=true;
     }
     else if (p<=42.07){
         nucltransK(0.928,0.020,9.2e-4,0.);
         next1498=go=true;
     }
     else if (p<=95.22){
         nucltransK(1.648,0.020,2.9e-4,1.2e-4);
         next778=go=true;
     }
     else{
         nucltransK(2.426,0.020,1.4e-4,5.1e-4);
	 return;
     }
     if(!go) return;
      
   }
   if (fLevelE==2096){
      fThlev=0.097e-12;
      bool go=false;
      float p=100.*GetRandom();
      if (p<=3.06){
         nucltransK(0.948,0.020,8.7e-4,0.);
         next1148=go=true;
      }
      else if (p<=98.55){
         nucltransK(1.317,0.020,4.4e-4,2.5e-5);
         next778=go=true;
      }
      else{
         nucltransK(2.096,0.020,1.8e-4,3.5e-4);
	 return;
      }
      if(!go) return;
   }
   if (fLevelE==1626|| next1626){
      bool go=false;
      fThlev=1.4e-12;
      float p=100.*GetRandom();
      if (p<=8.47){
         nucltransK(1.626,0.020,2.8e-4,1.3e-4);
         return;
      }
      else if (p<=98.58){
         nucltransK(0.848,0.020,1.2e-3,0.);
         next778=go=true;
      }
      else{
        nucltransK(0.128,0.020,1.3e-1,0.);
        next1498=go=true; 
      }      
      if(!go) return;
   }
   if (fLevelE==1498 || next1498){
      fThlev=0.78e-12;
      float p=100.*GetRandom();
      if (p<=29.73){
         nucltransK(1.498,0.020,3.3e-4,8.3e-5);
         return;
      }
      else{
         nucltransK(0.720,0.020,1.7e-3,0.);
         next778=true;
      }
      if(!next778) return;
   }
   if (fLevelE==1148 || next1148){
      fThlev=61.e-12;
      nucltransK(0.370,0.020,1.2e-2,0.);
      next778=true;
   }
   if (fLevelE==778 || next778){
      fThlev=3.67e-12;
      nucltransK(0.778,0.020,1.4e-3,0.);
      return;
   }
   return;
}

//-----------------------------------------------
void Decay::Zr92low()
{
  // Subroutine describes the deexcitation process in Zr92 nucleus after 2b-decay
  // of Mo92 to ground and excited 0+ and 2+ levels of Zr92

  fTclev=0.;
  bool next=false;

  if (fLevelE==1383){
     fThlev=0.17e-9;
     nucltransK(0.449,0.018,5.5e-3,0.);
     next=true;
  }
  if (fLevelE==934 || next){
     fThlev=5.0e-12;
     nucltransK(0.934,0.018,8.0e-4,0.);
     return;
  }
  return;
}

//-----------------------------------------------
void Decay::Ru100low()
{
  // Subroutine describes the deexcitation process in Ru100 nucleus after 2b-decay 
  // of Mo100 to ground and excited 0+ and 2+ levels of Ru100
  int  npg591=0;
  int  npg540=0;
  float cg,cK,p;
  bool next540=false;
  bool next1130=false;
  bool next1362=false;
  float ctet,stet1,stet2,ctet1,ctet2;
  fTclev=0.;
  if (fLevelE==1741){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=0.05){
        particle(3,1.741-0.022,1.741-0.022,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.022,0.022,0.,pi,0.,twopi,0,0);
        return;
     }
     else if(p<=59.00){
        nucltransK(1.201,0.022,6.2e-4,0.1e-4);
        next540=true;
     }
     else if(p<=59.03){
        particle(3,0.611-0.022,0.611-0.022,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.022,0.022,0.,pi,0.,twopi,0,0);
        next1130=true;
     }
     else{
        nucltransK(0.379,0.022,1.3e-2,0.);
        next1362=true;
     }
  }
  if (fLevelE==1362 || next1362){
     fThlev=1.2e-12;
     p=100.*GetRandom();
     if (p<=43.){
        nucltransK(1.362,0.022,4.2e-4,0.2e-4);
        return;
     }
     else{
        nucltransK(0.822,0.022,1.7e-3,0.);
        next540=true;
     }
  }
  if (fLevelE==1130 || next1130){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=0.02){
        particle(3,1.130-0.022,1.130-0.022,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.022,0.022,0.,pi,0.,twopi,0,0);
        return;
     }
     else{
        fEgamma=0.591;
	fEbindeK=0.022;
	cg=1.;
	cK=3.3e-3;
        p=GetRandom()*(cg+cK);
        if (p<=cg){
           particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
           npg591=fNbPart;
        }
        else{
           particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,fEbindeK,fEbindeK,0.,pi,0.,twopi,0,0);
        }
     }
    next540=true;
  }
  if (fLevelE==540 || next540){
     fThlev=11.e-12;
     fEgamma=0.540;
     fEbindeK=0.022;
     cg=1.;
     cK=4.4e-3;
     p=GetRandom()*(cg+cK);
     if (p<=cg){ 
        particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
        npg540=fNbPart;
     }
     else{
        particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,fEbindeK,fEbindeK,0.,pi,0.,twopi,0,0);
     }
     if (npg591!=0 && npg540!=0){ 
        float p591=sqrt(pow(fPmoment[0][npg591],2)+pow(fPmoment[1][npg591],2)+
                  pow(fPmoment[2][npg591],2));
	float p540=sqrt(pow(fPmoment[0][npg540],2)+pow(fPmoment[1][npg540],2)+
                  pow(fPmoment[2][npg540],2));
        // Coefficients in formula 1+a2*ctet**2+a4*ctet**4 are from:
        // R.D.Evans, "The Atomic Nucleus", Krieger Publ. Comp., 
        // 1985, p.240, 0(2)2(2)0 cascade.
        float a2=-3.0;
        float a4=4.0;
        do{
           phi1=twopi*GetRandom();
           ctet1=1.-2.*GetRandom();
           stet1=sqrt(1.-ctet1*ctet1);
           phi2=twopi*GetRandom();
           ctet2=1.-2.*GetRandom();
           stet2=sqrt(1.-ctet2*ctet2);
           ctet=ctet1*ctet2+stet1*stet2*cos(phi1-phi2);
        } while(GetRandom()*(1.+abs(a2)+abs(a4)) > 1.+a2*pow(ctet,2)+a4*pow(ctet,4));

        fPmoment[0][npg591]=p591*stet1*cos(phi1);
        fPmoment[1][npg591]=p591*stet1*sin(phi1);
        fPmoment[2][npg591]=p591*ctet1;
        fPmoment[0][npg540]=p540*stet2*cos(phi2);
        fPmoment[1][npg540]=p540*stet2*sin(phi2);
        fPmoment[2][npg540]=p540*ctet2;
        return;
      }
  }
  return;
}

//-----------------------------------------------
void Decay::Pd106low()
{
  // Subroutine describes the deexcitation process in Pd106 nucleus after 2b-decay 
  // of Cd106 to ground and excited 0+ and 2+ levels of Pd106
  bool next1128=false;
  bool next1134=false;
  bool next512=false;
  float p;
  fTclev=0.;
  
  if (fLevelE==1706){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=87.2){
        nucltransK(1.195,0.024,6.9e-4,6.7e-6);
        next512=true;
     }
     else{
        nucltransK(0.578,0.024,4.0e-3,0.);
        next1128=true;
     }
  } 
  if (fLevelE==1562){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=9.0){
        nucltransK(1.562,0.024,3.9e-4,1.1e-4);
        return;
     }
     else if(p<=95.0){
        nucltransK(1.050,0.024,1.0e-3,0.);
        next512=true;
     }
     else if(p<=96.1){
        nucltransK(0.434,0.024,7.7e-3,0.);
        next1128=true;
     }
     else{
        nucltransK(0.428,0.024,9.5e-3,0.);
        next1134=true;
     }
  } 
  if (fLevelE==1134 || next1134){
     fThlev=6.8e-12;
     p=100.*GetRandom();
     if (p<=5.7e-2){
        particle(3,1.110,1.110,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.024,0.024,0.,pi,0.,twopi,0.,0.);
        return;
     }
     else{
        nucltransK(0.622,0.024,3.2e-3,0.);
        next512=true;
     }
  }
  if (fLevelE==1128 || next1128){
     fThlev=3.12e-12;
     p=100.*GetRandom();
     if (p<=35.0){
        nucltransK(1.128,0.024,7.7e-4,0.);
        return;
     }
     else{
        nucltransK(0.616,0.024,3.4e-3,0.);
        next512=true;
     }
  }
  if (fLevelE==512 || next512){
     fThlev=12.1e-12;
     nucltransK(0.512,0.024,5.6e-3,0.);
     return;
  } 
}

//-----------------------------------------------
void Decay::Sn116low()
{
  // Subroutine describes the deexcitation process in Sn116 nucleus after 2b-decay 
  // of Cd116 to ground and excited 0+ and 2+ levels of Sn116
  fTclev=0.;
  bool next1294=false;
  bool next1757=false;
  float p;

  if (fLevelE==2225){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=37.){
        nucltransK(2.225,0.029,2.7e-4,3.4e-4);
        return;
     }
     else{
        nucltransK(0.932,0.029,1.5e-3,0.);
        next1294=true;
     }
  } 
  if (fLevelE==2112){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=54.9){
        nucltransK(2.112,0.029,3.1e-4,2.7e-4);
        return;
     }
     else if(p<=96.9){
        nucltransK(0.819,0.029,2.6e-3,0.);
        next1294=true;
     }
     else{
        nucltransK(0.355,0.029,1.8e-2,0.);
        next1757=true;
     }
  }  
  if (fLevelE==2027){
     fThlev=0.;
     nucltransK(0.733,0.029,2.7e-3,0.);
     next1294=true;
   }
  if (fLevelE==1757 || next1757){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=0.29){
        particle(3,1.757-0.029,1.757-0.029,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.029,0.029,0.,pi,0.,twopi,0.,0.);
        return; 
     }
     else{
        nucltransK(0.463,0.029,9.0e-3,0.);
        next1294=true;
     }
  }
  if (fLevelE==1294 || next1294){
     fThlev=0.36e-12;
     nucltransK(1.294,0.029,7.5e-4,0.5e-4);
     return;
  }

}
//-----------------------------------------------
void Decay::Cd112low()
{
  // Subroutine describes the deexcitation process in Cd112 nucleus after 2b-decay 
  // of Sn112 to ground 0+ and excited 2+ levels of Cd112
  fTclev=0.;
  float p;
  bool next618=false;
  bool next1224=false;
  bool next1312=false;
  bool next1469=false;
  
  if (fLevelE==1871){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=86.91){
        nucltransK(1.253,0.027,7.2e-4,1.5e-5);
        next618=true; 
     }
     else if (p<=89.88){
        nucltransK(0.559,0.027,4.9e-3,0.);
        next1312=true;
     }
     else{
        nucltransK(0.402,0.027,1.3e-2,0.);
        next1469=true;
     }
  }
  if (fLevelE==1469 || next1469){
     fThlev=2.7e-12;
     p=100.*GetRandom();
     if (p<=36.98){
        nucltransK(1.469,0.027,5.8e-4,7.1e-5);
        return;
     } 
     else if(p<=99.14){
        nucltransK(0.851,0.027,1.8e-3,0.);
        next618=true;
     }
     else{
        nucltransK(0.245,0.027,6.4e-2,0.);
        next1224=true;
     }
  } 
  if (fLevelE==1433){
     fThlev=1.9e-9;
     p=100.*GetRandom();
     if (p<=0.66){
        p=100.*GetRandom();
        if (p<=3.8){
           pair(0.411);
        }
        else{
          particle(3,1.406,1.406,0.,pi,0.,twopi,fTclev,fThlev);
          particle(1,0.027,0.027,0.,pi,0.,twopi,0.,0.); 
        }
        return;
     }
     else if (p<=39.36){
        nucltransK(0.816,0.027,1.8e-3,0.);
        next618=true;
     }
     else if (p<=60.61){
        particle(3,0.182,0.182,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.027,0.027,0.,pi,0.,twopi,0.,0.);
        next1224=true;
     }
     else{
        nucltransK(0.121,0.027,7.6e-1,0.);
        next1312=true;
     }
  }
  if (fLevelE==1312 || next1312){
      fThlev=2.0e-12;
      p=100.*GetRandom();
      if (p<=26.59){
         nucltransK(1.312,0.027,6.6e-4,2.6e-5);
         return;
      }
      else{
         nucltransK(0.695,0.027,2.8e-3,0.);
         next618=true;
      }
  } 
  if (fLevelE==1224 || next1224){
     fThlev=4.2e-12;
     p=100.*GetRandom();
     if (p<=0.12){
        p=100.*GetRandom();
        if (p<=0.4){
           pair(0.202);
        }
        else{
           particle(3,1.197,1.197,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.027,0.027,0.,pi,0.,twopi,0.,0.);
        }
      return;
     }
     else{
        nucltransK(0.607,0.027,3.9e-3,0.);
        next618=true;
     }
  }
  if (fLevelE==618 || next618){
     fThlev=6.51e-12;
     nucltransK(0.618,0.027,3.7e-3,0.);
     return;
  }
  return;
}

//-----------------------------------------------
void Decay::Te124low()
{
  // Subroutine describes the deexcitation process in Te124 nucleus after 2b-decay 
  // of Sn124 to ground and excited 0+ and 2+ levels of Te124.
  bool next603=false;  
  bool next1657=false;  
  bool next1326=false;  
  bool next1248=false;  
  
  float p;

  fTclev=0.;
  
  if (fLevelE==2182){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=9.53){
        nucltransK(2.182,0.032,7.0e-4,3.9e-4);
        return;
     }   
     else if (p<=92.38){
        nucltransK(1.580,0.032,7.7e-4,1.0e-4);
        next603=true;
     }
     else{
        nucltransK(0.857,0.032,2.3e-3,0.);
        next1326=true;
     }
  }
  if (fLevelE==2153){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=19.29){
        nucltransK(1.550,0.032,6.8e-4,1.0e-4);
        next603=true;
     }
     else{
        nucltransK(0.828,0.032,2.2e-3,0.);
        next1326=true;
     }
  }
  if (fLevelE==2092){
     fThlev=0.28e-12;
     p=100.*GetRandom();
     if (p<=97.97){
        nucltransK(1.488,0.032,8.3e-4,7.1e-5);
        next603=true;
     }
     else if (p<=98.24){
        nucltransK(0.844,0.032,2.1e-3,0.);
        next1248=true;
     }
     else{
        nucltransK(0.766,0.032,3.0e-3,0.);
        next1326=true;
     }
  }
  if (fLevelE==2039){
     fThlev=0.49e-12;
     p=100.*GetRandom();
     if (p<=34.26){
        nucltransK(2.039,0.032,6.7e-4,3.2e-4);
        return;
     }
     else if (p<=94.57){
        nucltransK(1.437,0.032,8.7e-4,5.5e-5);
        next603=true;
     }
     else if (p<=96.80){
        nucltransK(0.791,0.032,2.5e-3,0.);
        next1248=true;
     }
     else if (p<=99.04){
        nucltransK(0.714,0.032,4.0e-3,0.);
        next1326=true;
     }
     else{
        nucltransK(0.382,0.032,1.8e-2,0.);
        next1657=true; 
     }
  }
  if (fLevelE==1883){
     fThlev=0.76e-12;
     p=100.*GetRandom();
     if (p<=0.31){
        p=100.*GetRandom();
        if (p<=21.89){
           pair(0.861);
        }
        else{
           particle(3,1.851,1.851,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.032,0.032,0.,pi,0.,twopi,0.,0.);
        }
        return;
     }
     else if (p<=99.93){
        nucltransK(0.557,0.032,6.0e-3,0.);
        next1326=true;   
     }
     else{
        particle(3,0.194,0.194,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.032,0.032,0.,pi,0.,twopi,0.,0.);
        next1657=true;
     } 
  }
  if (fLevelE==1657 || next1657){
     fThlev=0.55-12;
     p=100.*GetRandom();
     if (p<=0.02){
        p=100.*GetRandom();
        if (p<=10.68){
           pair(0.636);
        }
        else{
           particle(3,1.626,1.626,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.032,0.032,0.,pi,0.,twopi,0.,0.);
        }
        return;
     }
     else{
        nucltransK(1.055,0.032,1.3e-3,0.);
        next603=true;
     }
  }
  if (fLevelE==1326 || next1326){
     fThlev=1.04e-12;
     p=100.*GetRandom();
     if (p<=13.84){
        nucltransK(1.326,0.032,8.3e-4,2.8e-5);
        return;
     }
     else{
        nucltransK(0.723,0.032,3.1e-3,0.);
        next603=true;
     }
  }
  if (next1248){
     fThlev=1.4e-12;
     nucltransK(0.646,0.032,4.1e-3,0.);
     next603=true;
  }
  if (fLevelE==603 || next603){
     fThlev=6.2e-12;
     nucltransK(0.603,0.032,4.9e-3,0.);
     return;
  }

  return;
}

//-----------------------------------------------
void Decay::Xe130low()
{
  // Subroutine describes the deexcitation process in Xe130 nucleus after 2b-decay 
  // of Te130 to ground 0+ and excited 2+ levels of Xe130
  float p;
  bool next536=false;
  bool next1122=false;
  fTclev=0.;
  
  if (fLevelE==1794){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=1.0){
        p=100.*GetRandom();
        if (p<=12.7){
           pair(0.772);
        }
        else{
           particle(3,1.759,1.759,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.035,0.035,0.,pi,0.,twopi,0.,0.); 
        }
        return;
     }
     else if (p<=86.3){
        nucltransK(1.258,0.035,1.0e-3,1.5e-5);
        next536=true; 
     }
     else{
        nucltransK(0.672,0.035,4.1e-3,0.);
        next1122=true;
     }
  }
  if (fLevelE==1122 || next1122){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=13.3){
        nucltransK(1.122,0.035,1.3e-3,9.0e-7);
        return;
     }
     else{
        nucltransK(0.586,0.035,7.5e-3,0.);
        next536=true;
     }   
  }
  if (fLevelE==536 || next536){
     fThlev=7.0e-12;
     nucltransK(0.536,0.035,7.4e-3,0.);
     return;
  }
  return;
}
//-----------------------------------------------
void Decay::Ba136low()
{
  // Subroutine describes the deexcitation process in Ba136 nucleus after 2b-decay
  // of Xe136 or Ce136 to ground and excited 0+ and 2+ levels of Ba136
  bool next819=false;
  bool next1551=false;
  float p;
  fTclev=0.;

  if (fLevelE==2400){
     fThlev=0.;
     nucltransK(1.581,0.037,1.0e-3,1.1e-4);
     next819=true;
  }
  if (fLevelE==2315){
     fThlev=0.;
     nucltransK(1.497,0.037,8.7e-4,7.7e-5);
     next819=true;
  }
  if (fLevelE==2223){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=4.3){
        nucltransK(2.223,0.037,7.8e-4,4.0e-4);
        return;
     }
     else if (p<=51.9){
        nucltransK(1.404,0.037,9.6e-4,4.8e-5);
        next819=true;
     }
     else{
        nucltransK(0.672,0.037,6.5e-3,0.);
        next1551=true;
     }
  }
  if (fLevelE==2141){
     fThlev=0.;
     nucltransK(1.323,0.037,1.0e-3,2.6e-5);
     next819=true;
  }
  if (fLevelE==2129){
     fThlev=0.051e-12;
     p=100.*GetRandom();
     if (p<=33.3){
        nucltransK(2.129,0.037,7.7e-4,3.6e-4);
        return;
     }
     else{
        nucltransK(1.310,0.037,1.4e-3,2.3e-5);
        next819=true;
     }
  }
  if (fLevelE==2080){
     fThlev=0.6e-12;
     p=100.*GetRandom();
     if (p<=35.4){
        nucltransK(2.080,0.037,7.6e-4,3.3e-4);
        return;
     }
     else if (p<=93.8){
        nucltransK(1.262,0.037,1.3e-3,1.4e-5);
        next819=true;
     }
     else{
        nucltransK(0.529,0.037,1.0e-2,0.);
        next1551=true;
     }
  }
  if (fLevelE==1579){
     fThlev=0.;
     nucltransK(0.761,0.037,3.4e-3,0.);
     next819=true;
  }
  if (fLevelE==1551 || next1551){
     fThlev=1.01e-12;
     p=100.*GetRandom();
     if (p<=52.1){
        nucltransK(1.551,0.037,8.4e-4,9.7e-5);
	return; 
     }
     else{
        nucltransK(0.733,0.037,4.5e-3,0.);
        next819=true;
     }
  }  
  if (fLevelE==819 || next819){
     fThlev=1.93e-12;
     nucltransK(0.819,0.037,2.9e-3,0.);
     return;
  }
  return;
}
//-----------------------------------------------

void Decay::Sm148low()
{
  // Subroutine describes the deexcitation process in Sm148 nucleus after 2b-decay 
  // of Nd148 to ground and excited 0+ and 2+ levels of Sm148  
  bool next550=false;
  float p;
  fTclev=0.;
  if (fLevelE==1455){
     fThlev=0.6e-12;
     p=100.*GetRandom();
     if (p<=42.){
        nucltransK(1.455,0.047,1.1e-3,0.3e-4);
        return;
     }
     else{
        nucltransK(0.904,0.047,2.8e-3,0.);
        next550=true; 
     }
  }
  if (fLevelE==550 || next550){
     fThlev=7.3e-12;
     nucltransK(0.550,0.047,9.0e-3,0.);
     return;
  }  
  return;
}

//-----------------------------------------------
void Decay::Sm150low()
{
  // Subroutine describes the deexcitation process in Sm150 nucleus after 2b-decay 
  // of Nd150 to ground and excited 0+ and 2+ levels of Sm150
  bool next334=false;
  bool next740=false;
  bool next1046=false;
  float p;

  fTclev=0.;

  if (fLevelE==1256){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=93.){
        nucltransK(0.922,0.047,2.6e-3,0.);
        next334=true;
     }
     else{
        nucltransK(0.210,0.047,1.7e-1,0.);
        next1046=true;
     }
  }
  if (fLevelE==1194){
     fThlev=1.3e-12;
     p=100.*GetRandom();
     if (p<=55.9){
        nucltransK(1.194,0.047,1.6e-3,0.1e-4);
        return;
     }
     else if (p<=96.9){
        nucltransK(0.860,0.047,3.2e-3,0.);
        next334=true;
     }
     else if (p<=98.7){
        nucltransK(0.453,0.047,1.5e-2,0.);
        next740=true;
     }
     else{
        nucltransK(0.420,0.047,1.9e-2,0.);
        fThlev=6.6e-12;
        nucltransK(0.439,0.047,1.7e-2,0.);
        next334=true;
     }
  } 
  if (fLevelE==1046 || next1046){
     fThlev=0.7e-12;
     p=100.*GetRandom();
     if (p<=7.0){
        nucltransK(1.046,0.047,2.0e-3,0.);
        return;
     }
     else if (p<=94.3){
        nucltransK(0.712,0.047,7.6e-3,0.);
        next334=true;
     }
     else if (p<=97.0){
        nucltransK(0.306,0.047,4.9e-2,0.);
        next740=true;
     }
     else{
        nucltransK(0.273,0.047,7.0e-2,0.);
        fThlev=6.6e-12;
        nucltransK(0.439,0.047,1.7e-2,0.);
        next334=true;
     }
  }
  if (fLevelE==740 || next740){
     fThlev=20.e-12;
     p=100.*GetRandom();
     if (p<=1.33){
        particle(3,0.740-0.047,0.740-0.047,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.047,0.047,0.,pi,0.,twopi,0.,0.);
        return;
     }
     else{
        nucltransK(0.407,0.047,2.0e-2,0.);
        next334=true;
     }
  }
  if (fLevelE==334 || next334){
     fThlev=48.5e-12;
     nucltransK(0.334,0.047,3.7e-2,0.);
     return;
  }
  return;
}

//-----------------------------------------------
void Decay::At214(float tcnuc)
{
  // Model for scheme of At214 decay (NDS 99(2003)649 and ENSDF at NNDC site on 9.12.2007).
  
  fThnuc=558.e-9;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float palpha=100.*GetRandom();
  if (palpha<=0.32){
     particle(47,8.267,8.267,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     nucltransK(0.563,0.091,8.6e-2,0.);
     return;
  }
  else if (palpha<=0.90){
     particle(47,8.478,8.478,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     nucltransK(0.348,0.091,7.9e-2,0.);
     return;
  }
  else if (palpha<=1.05){
     particle(47,8.505,8.505,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     nucltransK(0.320,0.091,3.9e-1,0.);
     return;
  }
  else{
     particle(47,8.819,8.819,0.,pi,0.,twopi,0,0);
     return;
  }
}

//-----------------------------------------------
void Decay::Po214(float tcnuc)
{
  // Scheme of Po214 decay (Nucl. Data Sheets 55(1988)665).
  // Alpha decay to excited level 1097.7 keV of Pb210 is neglected (6e-5%).
  
  fThnuc=164.3e-6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  float palpha=100.*GetRandom();
  if (palpha<=0.0104){
     particle(47,6.902,6.902,0.,pi,0.,twopi,0,0);
     nucltransK(0.800,0.088,1.1e-2,0.);
     return;
  }
  else{
     particle(47,7.687,7.687,0.,pi,0.,twopi,0,0);
     return;
  }
}

//-----------------------------------------------
void Decay::Rn218(float tcnuc)
{
  // Model for scheme of Rn218 decay (NDS 76(1995)127 and ENSDF at NNDC site on 9.12.2007).
  
  fThnuc=35.e-3;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float palpha=100.*GetRandom();
  if (palpha<=0.127){
     particle(47,6.532,6.532,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     nucltransK(0.609,0.093,2.1e-2,0.);
     return;
  }  
  else{
     particle(47,7.130,7.130,0.,pi,0.,twopi,0,0);
     return;
  }
}
//-----------------------------------------------
void Decay::Ra222(float tcnuc)
{
  // Model for scheme of Ra222 decay (NDS 107(2006)1027 and ENSDF at NNDC site on 9.12.2007)
  
  fThnuc=36.17;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float palpha=100.*GetRandom();
  if (palpha<=0.0042){
     particle(47,5.734,5.734,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     float p=100.*GetRandom();
     if (p<=66.35){
        nucltransK(0.840,0.098,2.9e-2,0.);
        return;
     }
     else{
        nucltransK(0.516,0.098,2.5e-2,0.);
        fThlev=0.;
        nucltransK(0.324,0.098,1.1e-1,0.);
        return;
     }
   }
   else if (palpha<=0.0083){
        particle(47,5.776,5.776,0.,pi,0.,twopi,0,0);
        fThlev=0.;
        float p=100.*GetRandom();
        if (p<=96.75){
           nucltransK(0.473,0.098,1.0e-2,0.);
           fThlev=0.;
           nucltransK(0.324,0.098,1.1e-1,0.);
           return;
        }
        else{
           nucltransK(0.144,0.098,1.9e-1,0.);
           fThlev=0.;
           nucltransK(0.329,0.098,1.1e-1,0.);
           fThlev=0.;
           nucltransK(0.324,0.098,1.1e-1,0.);
           return;
        }
    }
    else if (palpha<=0.0124){
        particle(47,5.917,5.917,0.,pi,0.,twopi,0,0);
        fThlev=0.;
        nucltransK(0.329,0.098,1.1e-1,0.);
        fThlev=0.;
        nucltransK(0.324,0.098,1.1e-1,0.);
        return;   
    }
    else if (palpha<=3.0635){
        particle(47,6.240,6.240,0.,pi,0.,twopi,0,0);
        fThlev=0.;
        nucltransK(0.324,0.098,1.1e-1,0.);
        return;
    }
    else{
        particle(47,6.559,6.559,0.,pi,0.,twopi,0,0);
        return;
    } 
    return;
}
/////////////////////////////////////////////////
//-----------------------------------------------
//--------Section of radioactive isotopes
//-----------------------------------------------
/////////////////////////////////////////////////
void Decay::Ac228(float tcnuc)
{
  // Scheme of 228Ac decay ("Table of Isotopes", 7th ed., 1978)  
  // VIT, 8.08.1992, 22.10.1995.  
  float p;
  bool next58=false,  next129=false, next328=false, next338=false;
  bool next332=false, next969=false, next964=false, next832=false;
  bool next874=false, next651=false, next1017=false,next10910=false;
  bool next988=false, next191=false, next1123=false,next1110=false;
  bool next952=false, next641=false, next1040=false,next1154=false;
  bool next1374=false,next1054=false,next178=false, next1724=false;
   
  fThnuc=22068;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  
  float pbeta=100.*GetRandom();
  if (pbeta<= 0.20){   // 0.20%
     beta(0.110,0.,0.);
     p=100*GetRandom();
     if (p<=31.){
        nucltransK(1.509,fEbindeK,3.2e-3,0.);
        next969=true;
     } 
     else{
        nucltransK(1.004,fEbindeK,3.0e-3,0.);
        next964=true;
     }
  }
   else if (pbeta<=0.40){   // 0.20%
     beta(0.127,0.,0.);
     p=100*GetRandom();
     if (p<=42.){
        nucltransK(1.952,fEbindeK,2.8e-3,1.6e-4);
        next58=true;
     }
     else if (p<=72.){
        nucltransK(1.823,fEbindeK,7.0e-3,1.7e-4);
        next129=true;
     }
     else{
        nucltransK(0.920,fEbindeK,0.7e+0,0.);
        next10910=true;
     }
  }
  else if (pbeta<=0.44){      // 0.04%
     beta(0.172,0.,0.);
     p=100*GetRandom();
     if (p<=33.){
        nucltransK(1.966,fEbindeK,1.2e-3,2.9e-4);
        return;
     } 
     else{
        nucltransK(1.907,fEbindeK,1.8e-3,2.7e-4);
        next58=true;
     }
  }
  else if (pbeta<=0.67){      // 0.23%
     beta(0.192,0.,0.);
     p=100*GetRandom();
     if (p<=31.4){
        nucltransK(1.887,fEbindeK,3.0e-3,1.3e-4);
        next58=true; 
     } 
     else if (p<=42.4){
        nucltransK(1.758,fEbindeK,1.0e-2,2.0e-4);
        next129=true;
     }
     else if (p<=44.3){
        nucltransK(1.549,fEbindeK,1.7e-3,0.6e-4);
        next338=true;
     }
     else if (p<=75.7){
        nucltransK(0.975,fEbindeK,4.0e-2,0.);
        next969=true;
     }
     else if (p<=87.8){       
        nucltransK(0.922,fEbindeK,4.5e-2,0.);
        next964=true; 
     }
     else if (p<=89.0){
        nucltransK(0.853,fEbindeK,5.5e-2,0.);
        next10910=true;
     }
     else if (p<=94.0){
        nucltransK(0.791,fEbindeK,6.5e-2,0.);
        next1154=true;
     }
     else if (p<=96.0){
        nucltransK(0.745,fEbindeK,7.5e-2,0.);
        next178=true;
     }
     else{
        nucltransK(0.220,fEbindeK,2.1e+0,0.);
        next1724=true;   
     }
  } 
  else if (pbeta<=0.77){     // 0.10%
     beta(0.237,0.,0.);
     p=100*GetRandom();
     if (p<=1.7){
        nucltransK(1.900,fEbindeK,3.0e-3,1.3e-4);
        return;
     }
     else if (p<=20.7){
        nucltransK(1.842,fEbindeK,8.0e-3,2.0e-4);
        next58=true;
     }
     else if (p<=21.7){
        nucltransK(1.712,fEbindeK,3.5e-3,0.8e-4);
        next129=true;
     }
     else if (p<=31.7){
        nucltransK(1.504,fEbindeK,1.7e-3,0.5e-4);
        next338=true;
     }
     else if (p<=67.9){
        nucltransK(0.884,fEbindeK,4.3e-3,0.);
        next1017=true; 
     }
     else{
        nucltransK(0.449,fEbindeK,1.6e-2,0.);
        next1054=true;
     }
  }
  else if (pbeta<=0.91){     // 0.14%
     beta(0.244,0.,0.);
     p=100*GetRandom();
     if (p<=17.0){
        nucltransK(1.835,fEbindeK,6.0e-3,2.0e-4);
        next58=true;
     }
     else if (p<=22.6){
        nucltransK(1.706,fEbindeK,1.0e-2,1.0e-4);
        next129=true;
     } 
     else if (p<=38.6){
        nucltransK(0.940,fEbindeK,1.0e-1,0.);
        next952=true;
     }
     else if (p<=47.6){
        nucltransK(0.924,fEbindeK,4.5e-2,0.);
        next969=true;
     }
     else if (p<=74.8){
        nucltransK(0.870,fEbindeK,5.2e-2,0.);
        next964=true;
     }
     else if (p<=76.8){
        nucltransK(0.739,fEbindeK,3.3e-1,0.);
        next1154=true;
     }
     else if (p<=77.8){
        nucltransK(0.693,fEbindeK,1.2e-1,0.);
        next178=true;  
     }
     else{
        nucltransK(0.461,fEbindeK,2.8e-1,0.);
        next1374=true;
     }
  }
  else if (pbeta<=1.02){     // 0.11%
     beta(0.377,0.,0.);
     p=100*GetRandom();
     if (p<=62.){
        nucltransK(1.702,fEbindeK,1.4e-3,1.1e-4);
        next58=true;
     }
     else{
        nucltransK(1.574,fEbindeK,1.7e-3,0.7e-4);
        next129=true;
     }
  }
  else if (pbeta<=1.31){    // 0.29%
     beta(0.393,0.,0.);
     p=100*GetRandom();
     if (p<=32.){
        nucltransK(1.686,fEbindeK,1.4e-3,1.6e-4);
        next58=true;
     }
     else if (p<=94.){
        nucltransK(1.557,fEbindeK,1.7e-3,0.6e-4);
        next129=true;
     }
     else{
        nucltransK(1.348,fEbindeK,2.0e-3,0.2e-4);
        next338=true;
     }
  }
  else if (pbeta<=2.81){    // 1.50%
     beta(0.413,0.,0.);
     next1724=true;
  }
  
  else if (pbeta<=4.91){    // 2.10%
     beta(0.448,0.,0.);
    
     p=100*GetRandom();
     if (p<=72.3){
        nucltransK(1.631,fEbindeK,7.4e-3,1.2e-4);
        next58=true;
     }
     else if (p<=98.0){
        nucltransK(1.502,fEbindeK,1.7e-3,0.5e-4);
        next129=true;
     }
     else{
        nucltransK(0.666,fEbindeK,7.5e-3,0.);
        next964=true;
     }
  }
  else if (pbeta<=6.41){    // 1.50%
     beta(0.454,0.,0.);
     p=100*GetRandom();
     if (p<=22.){
        nucltransK(1.625,fEbindeK,1.5e-3,0.8e-4);
        next58=true;
     }
     else if (p<=90.){
        nucltransK(1.496,fEbindeK,1.7e-3,0.5e-4);
        next129=true;
     }
     else if (p<=96.){
        nucltransK(1.287,fEbindeK,2.2e-3,0.2e-4);
        next338=true;
     }
     else{
        nucltransK(1.165,fEbindeK,2.6e-3,0.1e-4);
        next332=true;
     }
  }
  else if (pbeta<=11.11){   // 4.70%
     beta(0.491,0.,0.);
    
     p=100*GetRandom();
     if (p<=75.0){
        nucltransK(1.588,fEbindeK,4.7e-3,0.5e-4);
        next58=true;
     }
     else if (p<=95.0){
        nucltransK(1.459,fEbindeK,5.0e-3,0.3e-4);
        next129=true;
     }
     else if (p<=96.0){
        nucltransK(0.677,fEbindeK,2.2e-2,0.);
        next969=true;  
     }
     else if (p<=96.4){
        nucltransK(0.624,fEbindeK,8.0e-3,0.);
        next964=true;
     }
     else if (p<=97.4){
        nucltransK(0.555,fEbindeK,1.6e-1,0.);
        next10910=true;
     }
     else if (p<=99.5){
        nucltransK(0.523,fEbindeK,1.2e-2,0.);
        next1123=true;
     }
     else{
        nucltransK(0.420,fEbindeK,1.8e-2,0.);
        next1040=true;
     }
  }
  else if (pbeta<=11.91){  // 0.80%
     beta(0.494,0.,0.);
    
     p=100*GetRandom();
     if (p<=1.0){
        nucltransK(1.315,fEbindeK,1.8e-2,0.6e-4);
        next328=true;
     }
     else if (p<=51.7){
        nucltransK(1.247,fEbindeK,2.1e-2,0.4e-4);
        next338=true;
     }
     else if (p<=61.7){
        nucltransK(0.675,fEbindeK,7.0e-3,0.);
        next969=true;
     }
     else if (p<=70.7){
        nucltransK(0.620,fEbindeK,8.0e-3,0.);
        next964=true;  
     }  
     else{
        nucltransK(0.210,fEbindeK,7.9e-2,0.);
        next1374=true;
     }
  } 
  else if (pbeta<=13.11){  // 1.20%
     beta(0.499,0.,0.);
    
     p=100*GetRandom();
     if (p<=38.){
        nucltransK(1.638,fEbindeK,4.0e-3,0.5e-4);
        return;
     }
     else if (p<=96.){
        nucltransK(1.581,fEbindeK,1.1e-2,1.1e-4);
        next58=true;
     }
     else{
        nucltransK(0.516,fEbindeK,1.2e-2,0.);
        next1123=true;
     }
  }
  else if (pbeta<=13.18){  // 0.07%
     beta(0.590,0.,0.);
    
     p=100*GetRandom();
     if (p<=50.){
        nucltransK(1.169,fEbindeK,2.7e-3,0.2e-4);
        next191=true;
     }
     else if (p<=86.){
        nucltransK(0.378,fEbindeK,2.2e-2,0.);
        next1110=true;
     }
     else{
        nucltransK(0.373,fEbindeK,2.2e-2,0.);
        next988=true;
     }
  }
  else if (pbeta<=13.38){  // 0.20%
     beta(0.598,0.,0.);
    
     p=100*GetRandom();
     if (p<=12.){
        nucltransK(1.480,fEbindeK,1.8e-3,0.5e-4);
        next58=true;
     }
     else if (p<=18.){
        nucltransK(1.143,fEbindeK,2.6e-2,0.3e-4);
        next338=true;
     }
     else if (p<=32.){
        nucltransK(1.020,fEbindeK,1.0e-2,0.);
        next332=true;
     }
     else{
        nucltransK(0.571,fEbindeK,1.5e-1,0.);
        next641=true;
     }
  }
  else if (pbeta<=21.38){ // 8.00%
     beta(0.605,0.,0.);
    
     p=100*GetRandom();
     if (p<=0.4){
        nucltransK(1.136,fEbindeK,2.8e-3,0.2e-4);
        next338=true;
     }
     else if (p<=30.7){
        nucltransK(0.563,fEbindeK,5.0e-2,0.);
        next969=true;
     }
     else if (p<=46.7){
        nucltransK(0.509,fEbindeK,6.0e-2,0.);
        next964=true;
     }
     else if (p<=51.3){
        nucltransK(0.441,fEbindeK,3.0e-1,0.);
        next10910=true;
     }
     else if (p<=52.3){
        nucltransK(0.378,fEbindeK,4.5e-1,0.);
        next1154=true; 
     }
     else if (p<=53.0){
        nucltransK(0.357,fEbindeK,1.7e+0,0.);
        next988=true;
     }
     else{
        nucltransK(0.100,0.020,4.0e+0,0.);
        next1374=true; 
     }
  }
  else if (pbeta<=21.58){  // 0.20%
     beta(0.648,0.,0.);
    
     p=100*GetRandom();
     if (p<=22.){
        nucltransK(0.399,fEbindeK,2.0e-2,0.);
        next10910=true;  
     }
     else{
        nucltransK(0.314,fEbindeK,0.6e+0,0.);
        next988=true;
     }
  }
  else if (pbeta<=23.48){  //1.90%
     beta(0.687,0.,0.);
     next1054=true;
  }
  else if (pbeta<=24.98){  //1.5%
     beta(0.705,0.,0.);
     next1374=true;
  }
  else if (pbeta<=25.18){  // 0.20%
     beta(0.793,0.,0.);
    
     p=100*GetRandom();
     if (p<=27.0){
        nucltransK(1.017,fEbindeK,3.5e-3,0.);
        next328=true;
     }
     else if (p<=64.5){
        nucltransK(0.948,fEbindeK,3.7e-3,0.);
        next338=true;
     }
     else if (p<=95.0){
        nucltransK(0.825,fEbindeK,5.0e-3,0.);
        next332=true;
     }
     else{
        nucltransK(0.169,fEbindeK,1.4e-1,0.);
        next988=true;
     }
  }
  else if (pbeta<=25.98){  // 0.80%
     beta(0.910,0.,0.);
     next1040=true;
  }
  else if (pbeta<=29.58){  // 3.60%
     beta(0.968,0.,0.);
     next1110=true;
  }
  else if (pbeta<=34.58){   // 5.00%
     beta(0.983,0.,0.);
     next1154=true;
  }
  else if (pbeta<=39.98){   // 5.40%
      beta(1.014,0.,0.);
      next1123=true;
   }
  else if (pbeta<=40.08){   // 0.10%
     beta(1.077,0.,0.); 
     p=100*GetRandom();
     if (p<=74.){
        nucltransK(1.002,fEbindeK,3.5e-3,0.);
        next58=true;
     }    
     else if (p<=87.){
         nucltransK(0.664,fEbindeK,7.0e-3,0.);
         next338=true;
     }
     else{
         nucltransK(0.541,fEbindeK,1.1e-2,0.);
         next332=true;
     }
  }
  else if (pbeta<=42.08){   // 2.00%
     beta(1.115,0.,0.);
     next964=true;
  }
  else if (pbeta<=42.28){   // 0.20%
     beta(1.121,0.,0.);
     next1017=true;
  }
  else if (pbeta<=42.38){   // 0.10%
     beta(1.158,0.,0.);
     next651=true;
  }
  else if (pbeta<=74.33){   // 31.95%
     beta(1.168,0.,0.);
     next969=true;
  }
  else if (pbeta<=74.63){   // 0.30%
     beta(1.169,0.,0.);
     next641=true;
  }
  else if (pbeta<=74.86){   // 0.23%
     beta(1.193,0.,0.);
    
     p=100*GetRandom();
     if (p<=42.5){
        nucltransK(0.944,fEbindeK,1.1e-2,0.);
        return;
     }
     else if(p<=57.5){
        nucltransK(0.888,fEbindeK,7.5e-1,0.);
        next58=true;
     }
     else{
        nucltransK(0.616,fEbindeK,8.5e-3,0.);
        next328=true;
     }
  }
  else if (pbeta<=75.00){   // 0.14%
     beta(1.262,0.,0.);
     next874=true;
  }
  else if (pbeta<=75.20){   // 0.20%
     beta(1.618,0.,0.);
     next332=true;
  }
   else if (pbeta<=88.20){   // 13.00%
      beta(1.741,0.,0.);
      next338=true;
   }
  else if (pbeta<=89.00){   //0.80%
     beta(1.950,0.,0.);      
     next129=true;
  }
  else{ 
     beta(2.079,0.,0.);
     next58=true;
  }
 
  if (next1054){
    
     p=100*GetRandom();
     if (p<=3.0){
        nucltransK(1.054,fEbindeK,3.2e-2,0.);
        next338=true; 
     }
     else if (p<=15.0){
        nucltransK(0.498,fEbindeK,4.2e-2,0.);
        next952=true;
     }
     else if (p<=16.5){
        nucltransK(0.481,fEbindeK,2.5e-1,0.);
        next641=true;
     }
     else if (p<=58.0){
        nucltransK(0.328,fEbindeK,4.4e-2,0.);
        next1123=true;
     }
     else if (p<=82.0){
        nucltransK(0.282,fEbindeK,1.3e+0,0.);
        next1110=true;
     }
     else{
        nucltransK(0.224,fEbindeK,1.5e+0,0.);
        next1040=true;   
     }
  }
  if (next1724){
    
     p=100*GetRandom();
     if (p<=2.0){
        nucltransK(1.724,fEbindeK,3.7e-3,0.7e-4);
        return;
     } 
     else if (p<=15.0){
        nucltransK(1.666,fEbindeK,1.0e-2,1.2e-4);
        next58=true;
     }
     else if (p<=15.8){
        nucltransK(1.537,fEbindeK,4.0e-3,0.4e-4);
        next129=true;
     }
     else if (p<=89.0){
        nucltransK(0.755,fEbindeK,6.9e-2,0.);
        next969=true;
     }
     else{
        nucltransK(0.702,fEbindeK,9.5e-2,0.);
        next964=true;
     }
  }
  if (next1374){
    
     p=100*GetRandom();
     if (p<=1.0){
        nucltransK(1.374,fEbindeK,1.4e-2,0.7e-4);
        next58=true;
     }
     else if (p<=3.0){
        nucltransK(1.245,fEbindeK,2.0e-2,0.4e-4);
        next129=true;
     }
     else if (p<=64.8){
        nucltransK(0.463,fEbindeK,4.7e-2,0.);
        next969=true;
     }
     else if (p<=92.2){
        nucltransK(0.410,fEbindeK,8.3e-2,0.);
        next964=true;
     }
     else if (p<=98.2){
        nucltransK(0.341,fEbindeK,1.2e-1,0.);
        next10910=true;
     }
     else if (p<=99.0){
        nucltransK(0.308,fEbindeK,3.5e-2,0.);
        next1123=true;
     }
     else if (p<=99.6){
        nucltransK(0.264,fEbindeK,5.0e-2,0.);
        next1110=true;
     }
     else{
        nucltransK(0.258,fEbindeK,5.0e-2,0.);
        next988=true;
     }
  }
  if (next1040){
    
     p=100*GetRandom();
     if (p<=10.){
        nucltransK(1.040,fEbindeK,3.5e-3,0.);
        next129=true;
     } 
     else if (p<=68.){
        nucltransK(0.830,fEbindeK,1.8e-2,0.);
        next338=true;
     }
     else if (p<=79.){
        nucltransK(0.707,fEbindeK,1.0e-1,0.);
        next332=true;
     }
     else if (p<=97.){
        nucltransK(0.204,fEbindeK,9.0e-2,0.);
        next964=true;
     }
     else{
        nucltransK(0.136,fEbindeK,1.7e+0,0.);
        next10910=true;
     }
  }
  if (next178){
    
     nucltransK(0.178,fEbindeK,6.0e+1,0.);
     next964=true;
  }
  
  if (next1110){
     p=100*GetRandom();
     if (p<=11.0){
        nucltransK(1.110,fEbindeK,2.9e-3,0.1e-4);
        next58=true; 
     }   
     else if (p<=38.8){
        nucltransK(0.840,fEbindeK,1.4e-2,0.);
        next328=true;
     }
     else if (p<=84.5){
        nucltransK(0.772,fEbindeK,1.5e-2,0.);
        next338=true;
     }
     else if (p<=93.8){
        nucltransK(0.200,fEbindeK,9.5e-2,0.);
        next969=true;
     }
     else if (p<=98.9){
        nucltransK(0.146,fEbindeK,1.2e+0,0.);
        next964=true;
     }
     else{
        nucltransK(0.078,0.020,2.2e-1,0.);
        next10910=true;   
     }
  }
  if (next1154){
     p=100*GetRandom();
     if (p<=14.9){
        nucltransK(1.154,fEbindeK,7.5e-3,0.1e-4);
        return; 
     }
     else if (p<=25.8){
        nucltransK(1.096,fEbindeK,2.8e-2,0.1e-4);
        next58=true;
     }
     else if (p<=40.7){
        nucltransK(0.967,fEbindeK,2.0e-2,0.);
        next129=true;
     }
     else if (p<=61.4){
        nucltransK(0.322,fEbindeK,5.2e-1,0.);
        next832=true;
     }
     else if (p<=84.1){
        nucltransK(0.279,fEbindeK,1.3e+0,0.);
        next874=true;
     }
     else if (p<=93.0){
        nucltransK(0.185,fEbindeK,5.4e+1,0.);
        next969=true;
     }
     else if (p<=97.0){
        nucltransK(0.174,fEbindeK,1.4e+0,0.);
        next651=true; 
     }
     else{
        nucltransK(0.138,fEbindeK,4.9e+0,0.);
        next1017=true;       
     }
  }
  if (next1123){
     p=100*GetRandom();
     if (p<=1.4){
        nucltransK(1.123,fEbindeK,6.5e-2,0.1e-4);
        return;
     }
     else if (p<=5.4){
        nucltransK(1.065,fEbindeK,3.2e-3,0.);
        next58=true;
     }
     else if (p<=5.6){
        nucltransK(0.936,fEbindeK,0.1e+0,0.);
        next129=true;
     }
     else if (p<=73.8){
        nucltransK(0.795,fEbindeK,1.9e-2,0.);
        next328=true;
     }
     else if (p<=85.8){
        nucltransK(0.727,fEbindeK,1.0e-2,0.);
        next338=true;
     }
     else if (p<=98.1){
        nucltransK(0.154,fEbindeK,0.2e+0,0.);
        next969=true;      
     }
     else{
        nucltransK(0.100,0.020,3.3e-1,0.);
        next964=true;
     }
  }
  
  
  if (next988){
     p=100*GetRandom();
     if (p<=67.){
        nucltransK(0.988,fEbindeK,3.5e-3,0.);
        next129=true;
     }
     else{
        nucltransK(0.796,fEbindeK,5.2e-3,0.);
        next191=true;  
     }
  }
  if (next1017){
     p=100*GetRandom();
     if (p<=19.){
        nucltransK(1.017,fEbindeK,2.4e-2,0.);
        return;
     } 
     else if (p<=66.){
        nucltransK(0.958,fEbindeK,3.8e-3,0.);
        next58=true;
     }
     else if (p<=82.){
        nucltransK(0.688,fEbindeK,2.1e-2,0.);
        next328=true;
     }
     else{
        nucltransK(0.620,fEbindeK,1.4e-1,0.);
        next338=true;
     }
  }
  if (next651){
     p=100*GetRandom();
     if (p<=35.){
        nucltransK(0.651,fEbindeK,7.5e-3,0.);
        next328=true;
     }
     else{
        nucltransK(0.583,fEbindeK,9.5e-3,0.);
        next338=true;
     }
  }
  if (next952){
     p=100*GetRandom();
     if (p<=91.){
        nucltransK(0.894,fEbindeK,4.2e-3,0.);
        next58=true;
     }
     else if (p<=94.){
        nucltransK(0.624,fEbindeK,1.2e-1,0.);
        next328=true;
     }
     else{
        nucltransK(0.556,fEbindeK,3.5e-2,0.);
        next338=true;
     }
  }
  if (next641){
     p=100*GetRandom();
     if (p<=22.){
        nucltransK(0.641,fEbindeK,1.1e-1,0.);
        next328=true;
     }
     else{
        nucltransK(0.572,fEbindeK,1.3e-1,0.);
        next338=true; 
     }
  }
  if (next964){
     p=100*GetRandom();
     if (p<=73.){
        nucltransK(0.964,fEbindeK,9.2e-3,0.);
        next58=true;
     }
     else{
        nucltransK(0.836,fEbindeK,1.4e-2,0.);
        next129=true; 
     }
  }
  if (next969){
    
     p=100*GetRandom();
     if (p<=37.0){
        nucltransK(0.969,fEbindeK,1.0e-2,0.);
        return;
     }
     else if (p<=98.8){
        nucltransK(0.911,fEbindeK,1.2e-2,0.);
        next58=true;
     }
     else{
        nucltransK(0.782,fEbindeK,6.8e-2,0.);
        next129=true; 
     }
  }
  if (next10910){
     p=100*GetRandom();
     if (p<=21.){
        nucltransK(1.033,fEbindeK,9.5e-3,0.);
        next58=true;
     }
     else{
        nucltransK(0.904,fEbindeK,2.6e-2,0.);
        next129=true;
     }
  }
  if (next832){
     p=100*GetRandom();
     if (p<=4.){
        particle(3,0.721,0.721,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,fEbindeK,fEbindeK,0.,pi,0.,twopi,0.,0.);
        return;
     }
     else if (p<=27.){
        nucltransK(0.771,fEbindeK,1.7e-2,0.);
        next58=true;
     }
     else{
        nucltransK(0.504,fEbindeK,1.2e-2,0.);
        next328=true;
     }
  }
  if (next874){
     p=100*GetRandom();
     if (p<=15.){
        nucltransK(0.874,fEbindeK,1.3e-2,0.);
        return;
     } 
     else if (p<=26.){
        nucltransK(0.816,fEbindeK,0.5e+0,0.);
        next58=true;
     }
     else if (p<=60.){
        nucltransK(0.546,fEbindeK,1.1e-2,0.);
        next328=true;
     }
     else{
        nucltransK(0.479,fEbindeK,1.4e-2,0.);
        next338=true;
     }
  } 
  if (next332){
     p=100*GetRandom();
     if (p<=90.){
        nucltransK(0.332,fEbindeK,4.7e-1,0.);
        next129=true;
     }
     else{
        nucltransK(0.141,fEbindeK,0.9e+0,0.);
        next191=true;
     }
  }
  if (next338){
     p=100*GetRandom();
     if (p<=73.){
        nucltransK(0.338,fEbindeK,1.0e-2,0.);
        next58=true;
     }
     else{
        nucltransK(0.209,fEbindeK,7.9e-2,0.);
        next129=true;
     }
  }
  if (next328){
     p=100*GetRandom();
     if (p<=47.){
        nucltransK(0.328,fEbindeK,4.4e-2,0.);
        return; 
     } 
     else{
        nucltransK(0.270,fEbindeK,3.4e-2,0.);
        next58=true;
     } 
  }
  if (next191){
     nucltransK(0.191,fEbindeK,4.2e-1,0.);
     next129=true;
  }
  if (next129){
     fThlev=0.16e-9;
     nucltransK(0.129,0.020,2.7e+0,0.);
     next58=true;    
  }
  if (next58){
     fThlev=0.40e-9;
     nucltransK(0.058,0.020,1.2e+2,0.);
     return;    
  }
}

/////////////////////////////////////////////////
void Decay::As79(float tcnuc)
{
  //  Model for scheme of As79+Se79m decay 
  // ("Table of Isotopes", 8th ed.,1998 and Nucl. Data Sheets 96(2002)1).
  bool next960=false;
  bool next365=false;
  bool next528=false;
  bool next572=false;
  fThnuc=540.6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  float p;
  float pbeta=100.*GetRandom();
  if (pbeta<= 0.24){
     beta(1.192,0.,0.);
     p=100*GetRandom();
     if (p<=53.70){
        nucltransK(0.993,fEbindeK,4.0e-4,0.);
        next960=true;
     }
     else{
       nucltransK(0.724,fEbindeK,8.5e-4,0.);
       next365=true;
     }
  }
  else if (pbeta<= 0.67){
     beta(1.201,0.,0.);
     p=100*GetRandom();
     if (p<=68.97){
        nucltransK(0.715,fEbindeK,8.5e-4,0.);
        next365=true;
     }
     else {
        nucltransK(0.552,fEbindeK,1.5e-3,0.);
        next528=true;
     }
  }
  else if (pbeta<= 2.43){
      beta(1.306,0.,0.);
      p=100*GetRandom();
      if (p<=79.37){
         nucltransK(0.879,fEbindeK,5.5e-4,0.);
         next960=true;
      }
      else if (p<=94.45){
         nucltransK(0.447,fEbindeK,2.5e-3,0.);
         next528=true;
      }
      else{
         nucltransK(0.402,fEbindeK,3.2e-3,0.);
         next572=true; 
      }
  }
  else if (pbeta<= 2.69){
     beta(1.709,0.,0.);
     next572=true;
  }
  else if (pbeta<=3.79){
     beta(1.753,0.,0.);
     next528=true;
  }
  else if(pbeta<=5.24){
     beta(1.916,0.,0.);
     next365=true;
  }  
  else{
     beta(2.185,0.,0.);
     next960=true;
  }

  if (next572){
     fThlev=1.6e-11;
     p=100*GetRandom();
     if (p<=95.33){
        nucltransK(0.476,fEbindeK,1.0e-3,0.);
        next960=true;
     }
     else{
        nucltransK(0.207,fEbindeK,1.8e-2,0.);
        next365=true;
     }
  }
  if (next528){
     fThlev=3.1e-12;
     nucltransK(0.432,fEbindeK,2.7e-3,0.);
     next960=true;
  }
  if (next365){
     fThlev=9.4e-11;
     nucltransK(0.365,fEbindeK,2.0e-3,0.);
     return;
  }
  if (next960){
     fThlev=235.2;
     p=100*GetRandom();
     if (p<=99.944){
        nucltransK(0.096,fEbindeK,9.39,0.); //IT to Se79 (gs)
        return;
     }
     else{
        fZdtr=35;
        beta(0.247,fTclev,fThlev);// beta to Br79 (gs)
        return; 
     }
  }
}  
/////////////////////////////////////////////////
void Decay::Bi207(float tcnuc)
{
  /* Scheme of Bi207+Pb207m decay (Nucl. Data Sheets 70(1993)315)
     with NNDC corrections of 10.10.2000 and 20.04.2005.
     To describe the deexcitation processes in atomic shell 
     of Pb207, the information of PC Nuclear Data File retrieval 
     program and data base (last updated on 16-Aug-1994) was used.
     Gammas, beta+, e+e- pairs, K, L and M conversion electrons, 
     K, L and M X-rays and K and L Auger electrons are emitted.
  */ 
  fThnuc=1.0382166E+09;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fEbindeK=0.088;
  fEbindeL=0.015;
  fEbindeM=0.003;
  bool next57=false, next898=false;
  float p, cg=1.;
  int npe570=0,	npe1064=0, npg570=0, npg1064=0, np570=0, np1064=0;  
  double cK=0., cL=0., cM=0., cp=0., a2=0., a4=0.;
  float pdecay=100.*GetRandom();

  if (pdecay<=7.027){
     bool calc=false;
     if (GetRandom()<=0.70)  PbAtShell(88); 
     else PbAtShell(15);

     fThlev=0.;
     p=100*GetRandom();
     if (p<=98.13){
        fEgamma=1.770;
	cK=3.51e-3; cL=0.49e-3; cM=0.13e-3; cp=2.5e-4;
        calc=true;
        next57=true;
     }
     else{
        fEgamma=1.442;
	cK=2.7e-3; cL=0.4e-3; cM=0.1e-3; cp=0.2e-4;
        calc=true;
        next898=true;
     } 
     if (calc){
        p=GetRandom()*(cg+cK+cL+cM+cp);
        if (p<=cg) particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
        else if (p<=cg+cK){
           particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
           PbAtShell(88);
        }
        else if (p<=cg+cK+cL){
           particle(3,fEgamma-fEbindeL,fEgamma-fEbindeL,0.,pi,0.,twopi,fTclev,fThlev);
           PbAtShell(15);
        }
        else if (p<=cg+cK+cL+cM){
           particle(3,fEgamma-fEbindeM,fEgamma-fEbindeM,0.,pi,0.,twopi,fTclev,fThlev);
           PbAtShell(3);
        }
        else pair(fEgamma-1.022);
     }
  } 
  else if (pdecay<=90.992){
     p=GetRandom();
     if (p<=0.733) PbAtShell(88);
     else if (p<=0.931) PbAtShell(15);
     else PbAtShell(3);
  
     fThlev=0.80;
     fEgamma=1.064;
     cK=9.42e-2, cL=2.47e-2, cM=0.73e-2;
     p=GetRandom()*(cg+cK+cL+cM);
     if (p<=cg){
        particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
        npg1064=fNbPart;
     }
     else if (p<=cg+cK){
        particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
        npe1064=fNbPart;
        PbAtShell(88);
     }
     else if (p<=cg+cK+cL){
        particle(3,fEgamma-fEbindeL,fEgamma-fEbindeL,0.,pi,0.,twopi,fTclev,fThlev);
        npe1064=fNbPart;
        PbAtShell(15);
     }
     else {
        particle(3,fEgamma-fEbindeM,fEgamma-fEbindeM,0.,pi,0.,twopi,fTclev,fThlev);
        npe1064=fNbPart;
        PbAtShell(3);
     }
     next57=true;
  }
  else if (pdecay<=99.988){
     p=GetRandom();
     if (p<=0.7965) PbAtShell(88);
     else if (p<=0.9466)  PbAtShell(15);
     else PbAtShell(3);
     next57=true;
  }  
  else {
     beta(0.807,0.,0.);
     next57=true;
  }
  if (next898){
     p=100*GetRandom();
     if (p<=99.245){
        fEgamma=0.898;
	cK=2.01e-2; cL=0.34e-2;
        p=GetRandom()*(cg+cK+cL);
        if (p<=cg) particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
        else if (p<=cg+cK){
           particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
           PbAtShell(88);
        }
        else {
           particle(3,fEgamma-fEbindeL,fEgamma-fEbindeL,0.,pi,0.,twopi,fTclev,fThlev);
           PbAtShell(15); 
        }
        return; 
     }
     else{
        fEgamma=0.328;
	cK=0.2850; cL=0.0486; cM=0.0151;
        p=GetRandom()*(cg+cK+cL+cM);
        if (p<=cg) particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
        else if (p<=cg+cK){
           particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
           PbAtShell(88);
        }
        else if (p<=cg+cK+cL){
           particle(3,fEgamma-fEbindeL,fEgamma-fEbindeL,0.,pi,0.,twopi,fTclev,fThlev);
           PbAtShell(15);
        }
        else{
           particle(3,fEgamma-fEbindeM,fEgamma-fEbindeM,0.,pi,0.,twopi,fTclev,fThlev);
           PbAtShell(3);
        }
        next57=true;
     } 
  }
  if (next57){
     fThlev=130.5e-12;
     fEgamma=0.570;
     cK=1.55e-2; cL=0.45e-2; cM=0.15e-2;
     p=GetRandom()*(cg+cK+cL+cM);
     if (p<=cg){
        particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
        npg570=fNbPart;
     }
     else if (p<=cg+cK){
        particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
        npe570=fNbPart;
        PbAtShell(88);
     }
     else if (p<=cg+cK+cL){
        particle(3,fEgamma-fEbindeL,fEgamma-fEbindeL,0.,pi,0.,twopi,fTclev,fThlev);
        npe570=fNbPart;
        PbAtShell(15);
     }
     else {
        particle(3,fEgamma-fEbindeM,fEgamma-fEbindeM,0.,pi,0.,twopi,fTclev,fThlev);
        npe570=fNbPart;
        PbAtShell(3);
     }
  }
     if (npg1064!=0 && npg570!=0){
        a2=0.231; a4=-0.023;
        np1064=npg1064;
        np570=npg570; 
     }
     else if (npe1064!=0 && npg570!=0){
        a2=0.223; a4=-0.020;
        np1064=npe1064;
        np570=npg570;
     }
     else if (npg1064!=0 && npe570!=0){
        a2=0.275; a4=-0.012;
        np1064=npg1064;
        np570=npe570;
     }
     else if (npe1064!=0 && npe570!=0){
        a2=0.271; a4=-0.010;
        np1064=npe1064;
        np570=npe570;
     }
     else return;

     float ctet,stet1,stet2,ctet1,ctet2; 
     float p1064=sqrt(pow(fPmoment[0][np1064],2)+pow(fPmoment[1][np1064],2)
                     +pow(fPmoment[2][np1064],2));
     float p570=sqrt(pow(fPmoment[0][np570],2)+pow(fPmoment[1][np570],2)
                    +pow(fPmoment[2][np570],2));
     do{
        phi1=twopi*GetRandom();
        ctet1=1.-2.*GetRandom();
	stet1=sqrt(1.-ctet1*ctet1);
	phi2=twopi*GetRandom();
	ctet2=1.-2.*GetRandom();
	stet2=sqrt(1.-ctet2*ctet2);
	ctet=ctet1*ctet2+stet1*stet2*cos(phi1-phi2);
	p2=(3.*ctet*ctet-1.)/2.;
	p4=(35.*pow(ctet,4)-30.*ctet*ctet+3.)/8.;
     }while(GetRandom()*(1.+abs(a2)+abs(a4))>1.+a2*p2+a4*p4);

     fPmoment[0][np1064]=p1064*stet1*cos(phi1);
     fPmoment[1][np1064]=p1064*stet1*sin(phi1);
     fPmoment[2][np1064]=p1064*ctet1;
     fPmoment[0][np570]=p570*stet2*cos(phi2);
     fPmoment[1][np570]=p570*stet2*sin(phi2);
     fPmoment[2][np570]=p570*ctet2;
     return;
  
}
/////////////////////////////////////////////////
void Decay::Bi210(float tcnuc)
{
  // Scheme of Bi210 decay ("Table of Isotopes", 7th ed., 1978).
  // Three-figured labels correspond to energies of 
  // 206Tl excited levels in keV.
  // VIT, 14.08.1992, 22.10.1995; 30.10.2006.
  // Update to NDS 99(2003)949 and empirical correction 
  // to the beta shape, VIT, 28.10.2006. 
  fThnuc=433036.8;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float pdecay=100.*GetRandom();
  if (pdecay<=1.32e-4){
     float palfa=100.*GetRandom();
     if (palfa<=60.){
        particle(47,4.656,4.656,0.,pi,0.,twopi,0,0);
        fThlev=4.e-12;
        nucltransK(0.304,fEbindeK,3.9e-1,0.);
        return; 
     }
     else{
        particle(47,4.694,4.694,0.,pi,0.,twopi,0,0);
        fThlev=3.e-9;
        nucltransK(0.265,fEbindeK,1.6e-1,0.);
        return; 
     }
  }
  else beta1f(1.162,0.,0.,28.466,0.578,-0.658,0.);    
}
/////////////////////////////////////////////////
void Decay::Bi212(float tcnuc)
{
  // Scheme of Bi212 decay ("Table of Isotopes", 7th ed., 1978).
  // Two-, three- and four-figured labels correspond to energies of
  // 208Tl or 212Po excited levels in keV. Beta-alfa decays to 208Pb
  // are not considered (p=0.014%). 
  bool next400=false, next328=false, next727=false;
  fThnuc=3636.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  float pdecay,p, palfa, pbeta;
  pdecay=100*GetRandom();
  if (pdecay<=36.){  // 36% alfa to 208Tl
     palfa=100.*GetRandom();
     if (palfa<=1.10){  // 1.10% 
        particle(47,5.607,5.607,0.,pi,0.,twopi,0,0);
        p=100.*GetRandom();
        if (p<= 5.){
           nucltransK(0.493,0.086,2.8e-2,0.);
           return;
        }
        else if (p<=99.){
           nucltransK(0.453,0.086,0.18,0.);
           next400=true;
        }
        else{
           nucltransK(0.164,0.086,0.75,0.);
           next328=true;
        }
     }
     else if (palfa<=1.25){  // 0.15%
        particle(47,5.626,5.626,0.,pi,0.,twopi,0,0);
        p=100.*GetRandom();
        if (p<=68.){
           nucltransK(0.474,0.086,0.14,0.);
           return;
        }
        else if (p<=87.){
           nucltransK(0.434,0.086,0.14,0.);
           next400=true;
        }
        else{ 
           nucltransK(0.145,0.086,2.8,0.);
           next328=true;
        }
     } 
     else if (palfa<=2.92){  // 1.67%
        particle(47,5.769,5.769,0.,pi,0.,twopi,0,0);
        next328=true; 
     }
     else if (palfa<=72.80){ // 69.88%
        particle(47,6.051,6.051,0.,pi,0.,twopi,0,0);
        next400=true;
     }
     else{                  // 27.20% 
        particle(47,6.090,6.090,0.,pi,0.,twopi,0,0);
        return;
     }     

     if (next328){
        p=100.*GetRandom();
        if (p<=29.){
           nucltransK(0.328,0.086,0.33,0.);
           return;
        }
        else{
           nucltransK(0.288,0.086,0.53,0.);
           next400=true;
        }  
     }
     if (next400){
        fThlev=6.e-12;
        nucltransK(0.040,0.015,22.55,0.);
        return;
     }
  }
  else{                //64% beta to 212Po
     pbeta=64.*GetRandom();
     if (pbeta<=0.660){ // 0.660%
        beta(0.440,0.,0.);
        p=100*GetRandom();
        if (p<=17.){
           nucltransK(1.806,0.093,2.6e-2,1.7e-4);
           return;
        } 
        else{
           nucltransK(1.079,0.093,2.0e-2,0.);
           next727=true;
        }
     }
     else if (pbeta<=0.687){//0.027%
        beta(0.445,0.,0.);
        p=100*GetRandom();
        if (p<=35.){
           particle(3,1.708,1.708,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.093,0.093,0.,pi,0.,twopi,0,0);
           return;
        }
        else{
           nucltransK(1.074,0.093,7.0e-3,0.);
           next727=true;
        }
     }
     else if (pbeta<=0.937){ //0.250%
        beta(0.566,0.,0.);
        p=100*GetRandom();
        if (p<=28.){
           nucltransK(1.680,0.093,2.8e-3,1.0e-4);
           return;
        }
        else{
           nucltransK(0.952,0.093,4.5e-2,0.);
           next727=true;
        }
     }
     else if (pbeta<=2.837){ //1.900%
        beta(0.625,0.,0.);
        p=100*GetRandom();
        if (p<=80.){
           nucltransK(1.621,0.093,7.0e-3,1.2e-4);
           return;
        }
        else{
           nucltransK(0.893,0.093,4.5e-2,0.);
           next727=true;
        }
     }
     else if (pbeta<=4.337){ //1.500%
        beta(0.733,0.,0.);
        p=100*GetRandom();
        if (p<=22.){
           nucltransK(1.513,0.093,3.5e-3,0.7e-4);
           return;
        }
        else{
           nucltransK(0.786,0.093,4.1e-2,0.);
           next727=true;
        }
     } 
     else if (pbeta<=8.737){ //4.400%
        beta(1.519,0.,0.);
        next727=true;
     }
     else{                   // 55.263%
        beta(2.246,0.,0.); 
        return;
     }

     if (next727){
        nucltransK(0.727,0.093,1.7e-2,0.);
        return; 
     }
  }
}
//-----------------------------------------------
void Decay::Bi214(float tcnuc)
{
  // Scheme of Bi214 decay ("Table of Isotopes", 7th ed., 1978).
  // Two-, three- and four-figured labels correspond to energies of
  // 210Tl or 214Po excited levels in keV. Beta-alfa decays to 210Pb
  // are not considered (p=2.8e-3%).
  // VIT, 13.08.1992, 22.10.1995.
  // VIT, 31.05.2005, updated to NDS 99(2003)649. 
  // Not well known alpha decays to 210Tl levels with E>253 keV 
  // are omitted.
  bool next630=false, next609=false, next1015=false, next1275=false;
  bool next1378=false, next1415=false, next1543=false, next1661=false;
  bool next1713=false, next1730=false, next1743=false, next1764=false;
  bool next2448=false, next2266=false, next2209=false, next2204=false;
  bool next2119=false, next2088=false, next2011=false, next2017=false;
  bool next1995=false, next1847=false;
  fThnuc=1194.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  float p, pdecay, palfa;
   
  pdecay=100.*GetRandom();
  if (pdecay<=0.021){       // 0.021% alfa to 210Tl
     palfa=100.*GetRandom();
     if (palfa<=5.86){   // 5.86%
        particle(47,5.273,5.273,0.,pi,0.,twopi,0,0);
        nucltransK(0.191,0.086,1.3,0.);
        next630=true;
     }
     else if (palfa<=60.36){ // 54.50%
        particle(47,5.452,5.452,0.,pi,0.,twopi,0,0);
        next630=true;
     }
     else{                   // 39.64%
        particle(47,5.516,5.516,0.,pi,0.,twopi,0,0);
        return; 
     }
     if (next630){
        nucltransK(0.063,0.015,6.48,0.);
        return; 
     }
  }
  else{                    // 99.979% beta to 214Po
     float pbeta=100.*GetRandom();
     if (pbeta<=0.001){
        beta(0.088,0.,0.);
        nucltransK(3.184,0.093,4.0e-4,8.0e-4);
        return;
     }
     else if (pbeta<= 0.002){ // 0.001%
        beta(0.112,0.,0.);
        p=100.*GetRandom();
        if (p<=41.){
           nucltransK(3.160,0.093,4.0e-4,8.0e-4);
           return;
        }
        else{
           nucltransK(2.551,0.093,6.0e-4,4.6e-4);
           next609=true;
        }
     }
     else if (pbeta<= 0.003){ // 0.001%
        beta(0.129,0.,0.);
        nucltransK(3.143,0.093,4.0e-4,8.0e-4);
        return;
     }
     else if (pbeta<= 0.008){ // 0.005%
        beta(0.190,0.,0.);
        nucltransK(3.082,0.093,4.2e-4,8.0e-4);
        return;
     }
     else if (pbeta<= 0.012){ // 0.004%
        beta(0.218,0.,0.);
        p=100.*GetRandom();
        if (p<=60.0){
           nucltransK(3.054,0.093,1.0e-3,8.0e-4);
           return;
        }
        else if (p<=82.9){
           nucltransK(2.445,0.093,1.5e-3,5.2e-4);
           next609=true;
        }
        else{
           nucltransK(1.637,0.093,3.5e-3,0.5e-4);
           next1415=true;
        }
     }
     else if (pbeta<= 0.017){ // 0.005%
        beta(0.258,0.,0.);
        p=100.*GetRandom();
        if (p<=0.9){
           nucltransK(2.405,0.093,3.0e-3,4.1e-4);
           next609=true;
        }
        else if (p<=27.8){
           nucltransK(1.636,0.093,7.0e-3,1.2e-4);
           next1378=true;
        } 
        else if (p<=41.7){
           nucltransK(1.598,0.093,8.0e-3,1.1e-4);
           next1415=true;
        }
        else if (p<=61.8){
           nucltransK(1.471,0.093,9.0e-3,0.9e-4);
           next1543=true;
        }
        else{
           nucltransK(1.285,0.093,1.2e-2,0.5e-4);
           next1730=true;
        }
     }
     else if (pbeta<= 0.024){ // 0.007%
        beta(0.269,0.,0.);
        nucltransK(1.156,0.093,1.8e-2,3.0e-6);
        next1847=true;
     }
     else if (pbeta<= 0.034){ // 0.010%
        beta(0.272,0.,0.);
        p=100.*GetRandom();
        if (p<=84.6){
           nucltransK(3.000,0.093,1.4e-3,6.5e-4);
           return;
        }
        else{
           nucltransK(2.391,0.093,2.8e-3,4.0e-4);
           next609=true;
        }
     }
     else if (pbeta<= 0.043){ // 0.009%
        beta(0.286,0.,0.);
        p=100.*GetRandom();
        if (p<=83.0){
           nucltransK(2.377,0.093,2.8e-3,4.0e-4);
           next609=true;
        }
        else{
           nucltransK(1.711,0.093,6.0e-3,1.3e-4);
           next1275=true;
        }
     }
     else if (pbeta<= 0.060){ // 0.017%
        beta(0.293,0.,0.);  
        p=100.*GetRandom();
        if (p<=83.6){
           nucltransK(2.979,0.093,1.4e-3,6.5e-4);
           return;
        }
        else{
           nucltransK(2.369,0.093,2.8e-3,4.0e-4);
           next609=true;
        }
     }
     else if (pbeta<= 0.102){ // 0.042%
        beta(0.331,0.,0.); 
        p=100.*GetRandom();
        if (p<=52.1){
           nucltransK(2.331,0.093,3.0e-3,3.8e-4);
           next609=true;
        }
        else if (p<=71.7){
           nucltransK(1.666,0.093,7.0e-3,1.1e-4);
           next1275=true;
        }
        else{
           nucltransK(1.279,0.093,1.2e-2,1.5e-5);
           next1661=true;
        }
     }
     else if (pbeta<= 0.104){ // 0.002%
        beta(0.337,0.,0.);
        p=100.*GetRandom();
        if (p<=21.3){
           nucltransK(2.935,0.093,1.5e-3,6.3e-4);
           return;
        }
        else{
           nucltransK(2.325,0.093,3.0e-3,3.7e-4);
           next609=true;
        }
     }
     else if (pbeta<= 0.105){ // 0.001%
        beta(0.343,0.,0.);
        p=100.*GetRandom();
        if (p<=73.3){
           nucltransK(2.929,0.093,1.5e-3,6.2e-4);
           return;
        }
        else{
           nucltransK(2.319,0.093,3.0e-3,3.7e-4);
           next609=true;
        }
     }
     else if (pbeta<= 0.128){ // 0.023%
        beta(0.350,0.,0.);
        p=100.*GetRandom();
        if (p<=60.9){
           nucltransK(2.922,0.093,4.5e-4,8.0e-4);
           return;
        }
        else{
           nucltransK(2.312,0.093,3.0e-3,3.7e-4);
           next609=true;
        }
     }
     else if (pbeta<= 0.129){ // 0.001%
        beta(0.352,0.,0.);
        nucltransK(2.310,0.093,3.0e-3,3.7e-4);
        next609=true;
     }
     else if (pbeta<= 0.134){ // 0.005%
        beta(0.375,0.,0.);
        nucltransK(2.287,0.093,3.0e-3,3.6e-4);
        next609=true;
     }
     else if (pbeta<= 0.145){ // 0.011%
        beta(0.378,0.,0.);
        p=100.*GetRandom();
        if (p<=26.1){
           nucltransK(2.894,0.093,1.7e-3,6.1e-4);
           return; 
        }
        else if (p<=48.3){
           nucltransK(2.284,0.093,3.0e-3,3.6e-4);
           next609=true;
        }
        else if (p<=78.3){
           nucltransK(1.516,0.093,9.0e-3,6.7e-5);
           next1378=true;
        }
        else{
           nucltransK(0.626,0.093,7.5e-2,0.);
           next2266=true;
        }
     }
     else if (pbeta<= 0.156){ // 0.011%
        beta(0.392,0.,0.);
        p=100.*GetRandom();
        if (p<=87.6){
           nucltransK(2.880,0.093,1.7e-3,6.0e-4);
           return; 
        }
        else{
           nucltransK(2.270,0.093,3.0e-3,3.5e-4);
           next609=true;
        }
     }
     else if (pbeta<= 0.170){ // 0.014%
        beta(0.402,0.,0.); 
        p=100.*GetRandom();
        if (p<=63.5){
           nucltransK(2.260,0.093,3.0e-3,3.5e-4);
           next609=true;
        }
        else{
           nucltransK(1.595,0.093,8.0e-3,9.1e-5);
           next1275=true;
        }
     }
     else if (pbeta<= 0.184){ // 0.014%
        beta(0.411,0.,0.);
        p=100.*GetRandom();
        if (p<=2.7){
           nucltransK(2.861,0.093,1.8e-3,6.0e-4);
           return;
        }
        else if (p<=41.5){
           nucltransK(2.252,0.093,3.5e-3,3.4e-4);
           next609=true;
        }
        else{
           nucltransK(1.014,0.093,2.2e-2,0.);
           next1847=true;
        }
     }
     else if (pbeta<= 0.186){ // 0.002%
        beta(0.445,0.,0.);
        nucltransK(2.827,0.093,1.8e-3,5.8e-4);
        return;
     }
     else if (pbeta<= 0.222){ // 0.036%
        beta(0.486,0.,0.);
        p=100.*GetRandom();
        if (p<=15.4){
           nucltransK(2.786,0.093,1.8e-3,5.7e-4);
           return;
        }
        else if (p<=24.4){
           nucltransK(2.177,0.093,3.5e-3,3.1e-4);
           next609=true;
        }
        else if (p<=63.6){
           nucltransK(1.021,0.093,1.2e-2,0.);
           next1764=true;
        }
        else{
           nucltransK(0.939,0.093,1.4e-2,0.);
           next1847=true;
        }
     }
     else if (pbeta<= 0.258){ // 0.036%
        beta(0.502,0.,0.);
        p=100.*GetRandom();
        if (p<=54.6){
           nucltransK(2.770,0.093,1.8e-3,5.6e-4);
           return;
        }
        else if (p<=58.5){
           nucltransK(2.160,0.093,3.5e-3,3.1e-4);
           next609=true;
        }
        else{
           nucltransK(1.393,0.093,1.0e-2,3.6e-5);
           next1378=true;
        }
     }
     else if (pbeta<= 0.819){ // 0.561%
        beta(0.543,0.,0.);
        p=100.*GetRandom();
        if (p<=1.2){
           nucltransK(2.120,0.093,3.5e-3,2.9e-4);
           next609=true; 
        }
        else if (p<=6.0){
           nucltransK(1.067,0.093,2.1e-2,1.4e-7);
           next1661=true;
        } 
        else if (p<=71.5){
           nucltransK(0.964,0.093,1.4e-2,0.);
           next1764=true;
        }
        else if (p<=79.4){
           nucltransK(0.734,0.093,2.7e-2,0.);
           next1995=true;
        }
        else if (p<=82.6){
           nucltransK(0.525,0.093,6.0e-2,0.);
           next2204=true;
        }
        else if (p<=85.7){
           nucltransK(0.520,0.093,7.0e-2,0.);
           next2209=true;
        }
        else{
           nucltransK(0.281,0.093,3.3e-1,0.);
           next2448=true;
        }
     }
     else if (pbeta<= 1.097){ // 0.278%
        beta(0.553,0.,0.);
        p=100.*GetRandom();
        if (p<=0.7){
           nucltransK(2.719,0.093,2.0e-3,5.4e-4);
           return; 
        }
        else if (p<=33.6){
           nucltransK(2.110,0.093,3.5e-3,2.9e-4);
           next609=true;
        }
        else if (p<=41.8){
           nucltransK(1.341,0.093,1.0e-2,2.6e-5);
           next1378=true;
        }
        else if (p<=83.6){
           nucltransK(1.304,0.093,1.0e-2,1.9e-5);
           next1415=true;
        }
        else if (p<=90.7){
           nucltransK(0.976,0.093,2.6e-2,0.);
           next1743=true;
        }
        else if (p<=97.0){
           nucltransK(0.709,0.093,2.9e-2,0.);
           next2011=true;
        }
        else{
           nucltransK(0.600,0.093,9.0e-2,0.);
           next2119=true;
        }
     }
     else if (pbeta<= 1.151){ // 0.054%
        beta(0.573,0.,0.);
        p=100.*GetRandom();
        if (p<=5.3){
           nucltransK(2.699,0.093,2.0e-3,5.3e-4);
           return;
        }
        else{
           nucltransK(2.090,0.093,4.0e-3,2.8e-4);
           next609=true;
        }
     }
     else if (pbeta<= 1.204){ // 0.053%
        beta(0.574,0.,0.);
        p=100.*GetRandom();
        if (p<=29.5){
           nucltransK(1.156,0.093,1.6e-2,2.9e-6);
           next1543=true;
        }
        else if (p<=44.8){
           nucltransK(1.038,0.093,2.0e-2,0.);
           next1661=true;
        }
        else if (p<=63.3){
           nucltransK(0.935,0.093,2.8e-2,0.);
           next1764=true;
        }
        else if (p<=76.0){
           nucltransK(0.688,0.093,6.5e-2,0.);
           next2011=true;
        }
        else{
           nucltransK(0.494,0.093,7.0e-2,0.);
           next2204=true;
        }
     }
     else if (pbeta<= 1.467){ // 0.263%
        beta(0.577,0.,0.);
        p=100.*GetRandom();
        if (p<=11.9){
           nucltransK(2.695,0.093,2.0e-3,5.3e-4);
           return;
        }
        else if (p<=15.4){
           nucltransK(2.085,0.093,4.0e-3,2.8e-4);
           next609=true;
        }
        else if (p<=17.4){
           nucltransK(1.420,0.093,9.5e-3,4.2e-5);
           next1275=true;
        }
        else if (p<=48.7){
           nucltransK(1.317,0.093,1.2e-2,2.1e-5);
           next1378=true;
        }
        else if (p<=57.3){
           nucltransK(1.033,0.093,2.2e-2,0.);
           next1661=true; 
        }
        else if (p<=59.6){
           nucltransK(0.952,0.093,2.8e-2,0.);
           next1743=true;
        }
        else if (p<=72.3){
           nucltransK(0.930,0.093,1.5e-2,0.);
           next1764=true;
        }
        else if (p<=82.3){
           nucltransK(0.847,0.093,1.9e-2,0.);
           next1847=true; 
        }
        else if (p<=88.4){
           nucltransK(0.700,0.093,3.0e-2,0.);
           next1995=true; 
        }
        else if (p<=90.7){
           nucltransK(0.677,0.093,7.5e-2,0.);
           next2017=true;
        }
        else{
           nucltransK(0.486,0.093,8.0e-2,0.);
           next2209=true;
        }
     }
     else if (pbeta<= 1.594){ // 0.127%
        beta(0.610,0.,0.);
        p=100.*GetRandom();
        if (p<=0.2){
           nucltransK(2.662,0.093,2.0e-3,5.1e-4);
           return;
        }
        else if (p<=55.0){
           nucltransK(2.053,0.093,4.0e-3,2.6e-4);
           next609=true;
        }
        else if (p<=63.7){
           nucltransK(1.284,0.093,1.3e-2,1.6e-5);
           next1378=true;
        }
        else if (p<=95.6){
           nucltransK(1.119,0.093,1.9e-2,1.2e-6);
           next1543=true;
        }
        else{
           nucltransK(0.950,0.093,2.8e-2,0.);
           next1713=true;
        }
     }
     else if (pbeta<= 1.615){ // 0.021%
        beta(0.641,0.,0.);
        p=100.*GetRandom();
        if (p<=3.8){
           nucltransK(2.631,0.093,2.2e-3,5.0e-4);
           return;
        }
        else{
           nucltransK(2.022,0.093,4.0e-3,2.5e-4);
           next609=true;
        }
     }
     else if (pbeta<= 1.699){ // 0.084%
        beta(0.667,0.,0.);
        p=100.*GetRandom();
        if (p<=0.5){
           nucltransK(2.605,0.093,2.2e-3,4.9e-4);
           return;
        } 
        else if (p<=7.1){
           nucltransK(1.995,0.093,4.5e-3,2.4e-4);
           next609=true;
        }
        else if (p<=21.7){
           nucltransK(1.330,0.093,1.1e-2,2.4e-5);
           next1275=true;
        }
        else if (p<=44.2){
           nucltransK(0.943,0.093,2.8e-2,0.);
           next1661=true;
        }
        else if (p<=56.1){
           nucltransK(0.840,0.093,4.0e-2,0.);
           next1764=true;
        }
        else{
           nucltransK(0.396,0.093,1.3e-1,0.);
           next2209=true;
        }
     }
     else if (pbeta<= 1.752){ // 0.053%
        beta(0.727,0.,0.);
        p=100.*GetRandom();
        if (p<=77.4){
           nucltransK(1.936,0.093,3.0e-3,2.1e-4);
           next609=true; 
        }
        else{
           nucltransK(1.167,0.093,1.6e-2,3.6e-6);
           next1378=true;
        }
     }
     else if (pbeta<= 1.892){ // 0.140%
        beta(0.764,0.,0.); 
        p=100.*GetRandom();
        if (p<=49.8){
           nucltransK(1.899,0.093,5.0e-3,2.0e-4);
           next609=true;
        }
        else if (p<=84.8){
           nucltransK(1.130,0.093,1.8e-2,1.6e-6);
           next1378=true;
        }
        else if (p<=93.5){
           nucltransK(0.965,0.093,2.6e-2,0.);
           next1543=true;
        }
        else{
           nucltransK(0.497,0.093,7.0e-2,0.);
           next2011=true; 
        }
     }
     else if (pbeta<= 2.085){ // 0.193%
        beta(0.767,0.,0.);        
        p=100.*GetRandom();
        if (p<=3.0){
           nucltransK(2.505,0.093,2.4e-3,4.5e-4);
           return;
        }
        else if (p<=86.0){
           nucltransK(1.896,0.093,5.0e-3,2.0e-4);
           next609=true;
        }
        else if (p<=93.8){
           nucltransK(1.231,0.093,1.4e-2,9.2e-6);
           next1275=true;
        }
        else{
           nucltransK(0.962,0.093,2.6e-2,0.);
           next1543=true;
        }
     }
     else if (pbeta<= 3.417){ // 1.332%
        beta(0.790,0.,0.); 
        p=100.*GetRandom();
        if (p<=0.1){
           nucltransK(2.483,0.093,2.5e-3,4.4e-4);
           return;
        }
        else if (p<=16.7){
           nucltransK(1.873,0.093,5.0e-3,1.9e-4);
           next609=true;
        }
        else if (p<=51.0){
           nucltransKLM(1.208,0.093,1.6e-3,0.017,2.5e-4,0.004,8.0e-5,6.4e-5);
           next1275=true; 
        }
        else if (p<=56.9){
           nucltransK(1.105,0.093,1.0e-2,7.8e-7);
           next1378=true;
        }
        else if (p<=58.3){
           nucltransKLM(0.940,0.093,1.4e-2,0.017,2.5e-3,0.004,5.0e-4,0.);
           next1543=true;
        }
        else if (p<=70.7){
           nucltransK(0.821,0.093,3.7e-2,0.);
           next1661=true;
        }
        else if (p<=80.9){
           nucltransKLM(0.753,0.093,2.4e-2,0.017,4.4e-3,0.004,6.0e-4,0.);
           next1730=true;
        }
        else if (p<=81.4){
           nucltransKLM(0.635,0.093,3.6e-2,0.017,7.0e-3,0.004,7.0e-3,0.);
           next1847=true;
        }
        else if (p<=83.5){
           nucltransKLM(0.488,0.093,8.9e-3,0.017,1.5e-3,0.004,4.5e-4,0.);
           next1995=true;
        }
        else if (p<=84.8){
           nucltransK(0.394,0.093,1.3e-1,0.);
           next2088=true;
        }
        else{
           nucltransK(0.274,0.093,3.6e-1,0.);
           next2209=true;
        }
     }
     else if (pbeta<= 6.232){ // 2.815%
        beta(0.824,0.,0.);
        next2448=true;
     } 
     else if (pbeta<= 6.312){ // 0.080%
        beta(0.849,0.,0.);
        p=100.*GetRandom();
        if (p<=7.2){
           nucltransK(2.423,0.093,6.5e-4,9.2e-4);
           return;
        }
        else if (p<=24.5){
           nucltransK(1.814,0.093,1.2e-3,5.1e-4);
           next609=true;
        }
        else if (p<=65.4){
           nucltransK(1.046,0.093,2.8e-3,1.4e-6);
           next1378=true;
        }
        else if (p<=674.8){
           nucltransK(0.693,0.093,5.7e-3,0.);
           next1730=true;
        }
        else{
           nucltransK(0.659,0.093,3.5e-2,0.);
           next1764=true;
        }
     }
     else if (pbeta<= 6.872){ // 0.560%
        beta(0.979,0.,0.);
        p=100.*GetRandom();
        if (p<=54.6){
           nucltransK(2.293,0.093,3.0e-3,3.6e-4);
           return;
        }
        else if (p<=93.2){
           nucltransK(1.684,0.093,7.0e-3,1.2e-4);
           next609=true;
        }
        else if (p<=97.9){
           nucltransK(0.916,0.093,1.8e-2,0.);
           next1378=true;
        }
        else{
           nucltransK(0.878,0.093,2.0e-2,0.);
           next1415=true;
        }
     }
     else if (pbeta<= 7.073){ // 0.201%
        beta(1.006,0.,0.);
        next2266=true;
     }
     else if (pbeta<=12.802){ // 5.729%
        beta(1.068,0.,0.);
        next2204=true;
     }
        
     else if (pbeta<=13.635){ // 0.833%
        beta(1.079,0.,0.);
        p=100.*GetRandom();
        if (p<=4.1){
           nucltransK(2.193,0.093,2.0e-3,4.1e-4);
           return;
        }
        else if (p<=87.2){
           nucltransK(1.583,0.093,5.4e-3,8.7e-5);
           next609=true;
        }
        else if (p<=87.8){
           nucltransKLM(0.918,0.093,2.6e-3,0.017,4.1e-4,0.004,1.4e-4,0.);
           next1275=true;
        }
        else if (p<=92.5){
           nucltransKLM(0.815,0.093,1.9e-2,0.017,3.6e-3,0.004,1.4e-3,0.);
           next1378=true;
        }
        else{
           nucltransKLM(0.649,0.093,3.4e-2,0.017,6.0e-3,0.004,3.0e-3,0.);
           next1543=true; 
        }
     }
     else if (pbeta<=14.056){ // 0.421%
        beta(1.122,0.,0.);
        p=100.*GetRandom();
        if (p<=3.2){
           nucltransK(2.148,0.093,3.5e-3,3.0e-4);
           return; 
        } 
        else if (p<=89.0){
           nucltransK(1.539,0.093,8.0e-3,1.0e-4);
           next609=true;
        }
        else if (p<=93.1){
           nucltransK(0.873,0.093,3.5e-2,0.);
           next1275=true;
        }
        else{
           nucltransK(0.770,0.093,4.5e-2,0.);
           next1378=true; 
        }
     }
     else if (pbeta<=18.323){ // 4.267%
        beta(1.151,0.,0.);
        next2119=true;
     }
     else if (pbeta<=18.419){ // 0.096%
        beta(1.181,0.,0.);
        next2088=true;
     }
     else if (pbeta<=20.623){ // 2.204%
        beta(1.253,0.,0.);
        next2017=true;
     }
     else if (pbeta<=21.995){ // 1.372%
        beta(1.259,0.,0.);
        next2011=true;
     }
     else if (pbeta<=23.137){ // 1.142%
        beta(1.275,0.,0.);
        next1995=true;
     }
     else if (pbeta<=24.730){ // 1.593%
        beta(1.380,0.,0.);
        p=100.*GetRandom();
        if (p<=5.0){
           nucltransK(1.890,0.093,2.8e-3,2.6e-4);
           return;
        } 
        else if (p<=96.2){
           nucltransKLM(1.281,0.093,9.5e-3,0.017,1.6e-3,0.004,5.5e-4,1.6e-5);
           next609=true;  
        }
        else{
           nucltransKLM(0.616,0.093,5.6e-3,0.017,9.0e-4,0.004,3.0e-4,0.);
           next1275=true;
        }
     }
     else if (pbeta<=32.923){ // 8.193%
        beta(1.423,0.,0.);
        next1847=true;
     }
     else if (pbeta<=49.996){ // 17.073%
        beta(1.505,0.,0.);
        next1764=true;
     }
     else if (pbeta<=50.109){ // 0.113%
        beta(1.527,0.,0.);
        next1743=true;
     }
     else if (pbeta<=67.963){ // 17.854%
        beta(1.540,0.,0.);
        next1730=true;
     } 
     else if (pbeta<=68.113){ // 0.150%
        beta(1.540,0.,0.);
        next1713=true;
     }
     else if (pbeta<=68.834){ // 0.721%
        beta(1.609,0.,0.);
        next1661=true;
     }
     else if (pbeta<=71.789){ // 2.955%
        beta(1.727,0.,0.);
        next1543=true;
     }
     else if (pbeta<=72.600){ // 0.811%
        beta(1.855,0.,0.);
        next1415=true;
     }
     else if (pbeta<=80.042){ // 7.442%
        beta(1.892,0.,0.);
        next1378=true; 
     }
     else if (pbeta<=81.745){ // 1.703%
        beta(2.661,0.,0.); 
        next609=true;
     }
     else{       // 18.255%
        beta(3.270,0.,0.);
        return;
     }
     
     if (next2266){
        p=100.*GetRandom();
        if (p<=9.0){
           nucltransK(2.266,0.093,1.8e-3,4.4e-4);
           return; 
        }
        else if (p<=31.9){
           nucltransK(1.657,0.093,1.3e-3,3.9e-4);
           next609=true; 
        }
        else if (p<=36.9){
           nucltransK(0.991,0.093,1.3e-2,0.);
           next1275=true;
        } 
        else if (p<=54.8){
           nucltransK(0.723,0.093,2.8e-2,0.);
           next1543=true;
        }
        else if (p<=90.6){
           nucltransK(0.537,0.093,6.0e-2,0.);
           next1730=true;
        }
        else{
           nucltransK(0.502,0.093,7.0e-2,0.);
           next1764=true; 
        }
     }
     if (next2448){
        p=100.*GetRandom();
        if (p<=54.3){
           nucltransK(2.448,0.093,6.0e-4,9.4e-4);
           return;
        }
        else if (p<=66.8){
           nucltransK(1.838,0.093,1.2e-3,5.3e-4);
           next609=true; 
        }
        else if (p<=68.6){
           nucltransKLM(1.173,0.093,4.4e-3,0.017,8.3e-4,0.004,2.8e-4,6.2e-6);
           next1275=true; 
        }
        else if (p<=78.1){
           nucltransKLM(1.070,0.093,2.0e-3,0.017,3.1e-4,0.004,1.0e-4,5.5e-6);
           next1378=true;
        }
        else if (p<=80.8){
           nucltransKLM(1.032,0.093,2.1e-3,0.017,3.3e-4,0.004,1.1e-4,0.); 
           next1415=true;
        }
        else if (p<=83.7){
           nucltransKLM(0.904,0.093,2.7e-3,0.017,4.2e-4,0.004,1.4e-4,0.);
           next1543=true;
        }
        else if (p<=94.4){
           nucltransKLM(0.786,0.093,3.5e-3,0.017,5.5e-4,0.004,1.9e-4,0.);
           next1661=true;
        }
        else if (p<=96.0){
           nucltransKLM(0.705,0.093,4.3e-3,0.017,6.9e-4,0.004,2.2e-4,0.);
           next1743=true; 
        }
        else if (p<=98.8){
           nucltransKLM(0.683,0.093,4.6e-3,0.017,7.3e-4,0.004,2.4e-4,0.);
           next1764=true; 
        }
        else{
           nucltransKLM(0.453,0.093,8.0e-2,0.017,1.7e-2,0.004,5.6e-3,0.);
           next1995=true;  
        }
     }
     if (next2209){
        p=100.*GetRandom();
        if (p<=82.1){
           nucltransK(1.599,0.093,8.0e-3,9.2e-5);
           next609=true;
        }
        else{
           nucltransK(0.934,0.093,2.8e-2,0.);
           next1275=true;
        } 
     }
     if (next2204){
        p=100.*GetRandom();
        if (p<=87.9){
           nucltransK(2.204,0.093,3.5e-3,3.2e-4);
           return;
        }
        else if (p<=92.2){
           nucltransK(1.595,0.093,5.5e-3,9.1e-5);
           next609=true;
        }
        else if (p<=94.2){
           nucltransKLM(0.826,0.093,2.9e-2,0.017,5.0e-3,0.004,2.0e-3,0.);
           next1378=true;
        }
        else if (p<=94.5){
           nucltransKLM(0.789,0.093,3.3e-2,0.017,5.7e-3,0.004,1.9e-3,0.);
           next1415=true;
        }
        else if (p<=95.3){
           nucltransKLM(0.661,0.093,3.3e-2,0.017,6.0e-3,0.004,2.0e-3,0.);
           next1543=true;
        } 
        else if (p<=96.9){
           nucltransKLM(0.543,0.093,5.0e-2,0.017,1.1e-2,0.004,9.0e-3,0.);
           next1713=true;
        }
        else if (p<=99.0){
           nucltransKLM(0.474,0.093,7.0e-2,0.017,1.5e-2,0.004,4.9e-3,0.);
           next1730=true;
        }
        else{
           nucltransKLM(0.461,0.093,1.4e-1,0.017,2.3e-2,0.004,7.4e-3,0.); 
           next1743=true;
        }
     }
     if (next2119){ 
        p=100.*GetRandom();
        if (p<=26.7){
           nucltransK(2.119,0.093,3.5e-3,2.9e-4);
           return;
        }
        else if (p<=76.4){
           nucltransK(1.509,0.093,6.3e-3,6.5e-5);
           next609=true;
        }
        else if (p<=77.3){
           nucltransK(0.741,0.093,3.1e-2,0.);
           next1378=true;
        }
        else if (p<=89.0){
           nucltransKLM(0.703,0.093,4.5e-2,0.017,7.6e-3,0.004,2.6e-3,0.);
           next1415=true;
        }
        else{
           nucltransKLM(0.389,0.093,2.1e-1,0.017,3.7e-2,0.004,1.2e-2,0.);
           next1730=true;
        }
     }
     if (next2088){
        p=100.*GetRandom();
        if (p<=40.5){
           nucltransK(1.479,0.093,9.0e-3,5.7e-5);
           next609=true;
        }
        else{
           nucltransK(0.711,0.093,6.0e-2,0.);
           next1378=true;
        }
     }
     if (next2017){
        p=100.*GetRandom();
        if (p<=0.27){
           p=100.*GetRandom();
           if (p<=95.08){
              particle(3,1.923,1.923,0.,pi,0.,twopi,fTclev,fThlev);
              particle(1,0.093,0.093,0.,pi,0.,twopi,0,0);
           }
           else pair(0.995);   
           return;
        }
        else if (p<=98.01){
           nucltransKLM(1.408,0.093,3.1e-3,0.017,5.7e-4,0.004,1.9e-4,5.7e-5);
           next609=true;
        } 
        else if (p<=99.41){
           nucltransKLM(0.640,0.093,1.4e-2,0.017,3.7e-3,0.004,1.2e-3,0.);
           next1378=true;
        }
        else if (p<=99.75){
           nucltransKLM(0.356,0.093,4.6e-2,0.017,2.4e-2,0.004,8.1e-3,0.);
           next1661=true;
        }
        else{
           nucltransKLM(0.253,0.093,6.9e-1,0.017,1.2e-1,0.004,3.8e-2,0.);
           next1764=true; 
        }
     }
     if (next2011){
        p=100.*GetRandom();
        if (p<=3.4){
           nucltransK(2.011,0.093,2.2e-3,3.2e-4);
           return;
        }
        else if (p<=94.7){
           nucltransKLM(1.402,0.093,4.4e-3,0.017,7.7e-4,0.004,2.3e-4,3.8e-5);
           next609=true;
        }
        else if (p<=98.7){
           nucltransK(0.633,0.093,4.5e-2,0.);
           next1378=true;
        }
        else{
           nucltransK(0.595,0.093,5.3e-2,0.);
           next1415=true;
        }
     }
     if (next1995){
        p=100.*GetRandom();
        if (p<=60.3){
           nucltransKLM(1.385,0.093,1.3e-3,0.017,2.0e-4,0.004,6.4e-5,1.8e-4);
           next609=true;
        }
        else if (p<=90.8){
           nucltransKLM(0.720,0.093,1.1e-2,0.017,2.7e-3,0.004,9.3e-4,0.);
           next1275=true;
        }
        else if (p<=93.5){
           nucltransK(0.617,0.093,6.8e-3,0.);
           next1378=true; 
        }
        else{
           nucltransKLM(0.333,0.093,2.0e-2,0.017,3.5e-3,0.004,1.1e-3,0.);
           next1661=true;
        }
     }
     if (next1764){
        p=100.*GetRandom();
        if (p<=87.61){
           nucltransK(1.764,0.093,6.0e-3,1.5e-4);
           return;
        }
        else if (p<=97.01){
           nucltransKLM(1.155,0.093,1.2e-2,0.017,2.0e-3,0.004,6.3e-4,2.8e-6);
           next609=true;
        }
        else if (p<=99.06){
           nucltransKLM(0.387,0.093,1.3e-1,0.017,2.8e-2,0.004,8.9e-3,0.);
           next1378=true;
        }
        else if (p<=99.97){
           nucltransKLM(0.349,0.093,2.9e-1,0.017,5.0e-2,0.004,1.6e-2,0.);
           next1415=true;
        }
        else{
           nucltransKLM(0.221,0.093,5.6e-1,0.017,1.6e-1,0.004,5.4e-2,0.);
           next1543=true;
        }
     }
     if (next1743){
        nucltransKLM(1.134,0.093,4.7e-3,0.017,9.0e-4,0.004,2.9e-4,2.8e-6);
        next609=true;
     }
     if (next1730){
        p=100.*GetRandom();
        if (p<=15.66){
           nucltransK(1.730,0.093,2.7e-3,1.9e-4);
           return;
        }
        else if (p<=97.92){
           nucltransKLM(1.120,0.093,1.3e-2,0.017,2.2e-3,0.004,6.7e-4,1.2e-6);
           next609=true;
        }
        else if (p<=99.55){
           nucltransKLM(0.455,0.093,1.0e-2,0.017,1.7e-3,0.004,5.3e-4,0.);
           next1275=true;
        }
        else{
           nucltransKLM(0.352,0.093,1.6e-1,0.017,3.7e-2,0.004,1.2e-2,0.);
           next1378=true; 
        }
     }
     if (next1713){
        p=100.*GetRandom();
        if (p<=65.36){
           nucltransK(1.104,0.093,1.2e-2,7.6e-7);
           next609=true;
        }
        else{
           nucltransK(0.698,0.093,6.0e-2,0.);
           next1015=true;
        }
     }
     if (next1661){
        p=100.*GetRandom();
        if (p<=78.23){
           nucltransK(1.661,0.093,3.0e-3,1.5e-4);
           return;
        }
        else{
           nucltransK(1.052,0.093,1.5e-2,0.);
           next609=true;
        }
     }
     if (next1543){
        p=100.*GetRandom();
        if (p<=6.00){
           nucltransK(1.543,0.093,2.7e-3,1.1e-4);
           return;
        }
        else if (p<=99.26){
           nucltransKLM(0.934,0.093,2.0e-2,0.017,3.5e-3,0.004,1.2e-3,0.);
           next609=true; 
        }
        else if (p<=99.38){
           nucltransKLM(0.528,0.093,2.0e-2,0.017,6.4e-3,0.004,2.1e-3,0.);
           next1015=true;
        }
        else{
           nucltransKLM(0.269,0.093,3.3e-2,0.017,5.8e-3,0.004,1.8e-3,0.);
           next1275=true;
        }
     }
     if (next1415){
        fThlev=99.e-12;
        p=100.*GetRandom();
        if (p<=28.05){
           p=100.*GetRandom();
           if (p<=99.63){
              particle(3,1.322,1.322,0.,pi,0.,twopi,fTclev,fThlev);
              particle(1,0.093,0.093,0.,pi,0.,twopi,0.,0.); 
           }
           else{
              pair(0.393);
           }
           return;
        }
        else{
           nucltransKLM(0.806,0.093,8.7e-3,0.017,2.0e-3,0.004,6.6e-4,0.);
           next609=true;
        }
     }
     if (next1847){
        p=100.*GetRandom();
        if (p<=25.7){
           nucltransK(1.847,0.093,2.5e-3,2.4e-4);
           return;
        }
        else if (p<=97.0){
           nucltransKLM(1.238,0.093,1.0e-2,0.017,1.8e-3,0.004,5.5e-4,1.0e-5);
           next609=true;
        }
        else if (p<=97.3){
           nucltransKLM(0.832,0.093,8.2e-3,0.017,1.8e-3,0.004,6.4e-4,0.);
           next1015=true; 
        }
        else if (p<=98.2){
           nucltransKLM(0.573,0.093,6.4e-3,0.017,1.1e-3,0.004,3.4e-4,0.);
           next1275=true;
        }
        else{
           nucltransKLM(0.470,0.093,8.0e-2,0.017,1.6e-2,0.004,5.1e-3,0.);
           next1378=true;
        }
     }
     if (next1378){
        p=100.*GetRandom();
        if (p<=44.47){
           nucltransKLM(1.378,0.093,3.3e-3,0.017,5.9e-4,0.004,2.0e-4,4.8e-5);
           return;
        }
        else{
           nucltransKLM(0.768,0.093,1.3e-2,0.017,2.7e-3,0.004,9.0e-4,0.);
           next609=true;
        }
     } 
     if (next1275){
        fThlev=0.;
        nucltransKLM(0.665,0.093,4.8e-3,0.017,7.7e-4,0.004,2.5e-4,0.);
        next609=true;
     }
     if (next1015){
        fThlev=0.;
        nucltransKLM(0.406,0.093,3.5e-2,0.017,1.5e-2,0.004,5.1e-3,0.);
        next609=true;
     }
     if (next609){
        fThlev=0.;
        nucltransKLM(0.609,0.093,1.5e-2,0.017,4.2e-3,0.004,1.4e-3,0.);
        return;
     }
  }
}
//-----------------------------------------------

void Decay::Sc48(float tcnuc)
{
  // Scheme of Sc48 decay("Table of Isotopes", 8 ed.,1996 + NDS 68(1993)1).
  // Three-figured labels correspond to energies of 48Ti excited
  // levels in keV. 
  // VIT, 7.05.1998; 13.08.2007 update to NDS 107(2006)1747.
  fThnuc=1.57212e5;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  bool next=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  fZdtr=22;
  if (pbeta<=10.02){
     beta(0.483,0.,0.);
     fThlev=1.4e-12;
     p=100.*GetRandom();      
     if (p<=24.14){
        nucltransK(1.213,0.005,8.8e-5,0.1e-4);
        next=true;
     }
     else{
        nucltransK(0.175,0.005,4.5e-3,0.);
        fThlev=221.e-15;
        nucltransK(1.038,0.005,1.1e-4,0.);
        next=true;
     }
  }
  else{
     beta(0.659,0.,0.);
     fThlev=221.e-15;
     nucltransK(1.038,0.005,1.1e-4,0.);
     next=true;   
  } 
  if (next){
     fThlev=0.762e-12;
     nucltransK(1.312,0.005,9.7e-5,0.1e-4);
     fThlev=4.04e-12;
     nucltransK(0.984,0.005,1.3e-4,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Co60(float tcnuc)
{
  // Scheme of Co60 decay ("Table of Isotopes", 7th ed., 1978).
  // Four-figured labels correspond to energies of 60Ni excited levels in keV.
  // VIT, 3.08.1992, 22.10.1995.
  // Updated to NDS 100(2003)347, VIT, 16.10.2006;
  // angular correlation of 1173 and 1333 keV gammas, 
  // L.Pandola + VIT, 18.10.2006;
  // 2nd forbidden unique shape for beta decay to 1333 keV level, 
  // VIT, 27.10.2006.
  fThnuc=0.166344192e+09;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float pbeta, p;
  int npg1173=0, npg1333=0;
  bool next1333=false;
  pbeta=100.*GetRandom();
  if (pbeta<=99.880){
     beta(0.318,0.,0.);
     fThlev=0.3e-12;
     p=100.*GetRandom();
     if (p<=0.000002){
        nucltransK(2.506,0.008,8.6e-5,0.);
        return;
     } 
     else if (p<=99.992449){
        fEgamma=1.173;
        fEbindeK=0.008;
        float cg=1., cK=1.7e-4, cp=6.2e-6;
	p=GetRandom()*(cg+cK+cp);
        if (p<=cg){
           particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
           npg1173=fNbPart;
        }
        else if (p<=cg+cK){
           particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,fEbindeK,fEbindeK,0.,pi,0.,twopi,0.,0.); 
        }
        else{
           pair(fEgamma-1.022);
        }
        next1333=true;
     }
     else{
        nucltransK(0.347,0.008,5.5e-3,0.);
        fThlev=0.;
        p=GetRandom()*100;
        if (p<=13.64){
           nucltransK(2.159,0.008,4.9e-5,3.9e-4);
           return;
        }
        else{
           nucltransK(0.826,0.008,3.3e-4,0.);
           next1333=true;
        }
     }
  }
  else{
     beta2f(1.491,0.,0.,2,3.333333,1.,0.,0.);
     next1333=true;
  }
  if (next1333){
     fThlev=0.9e-12;
     fEgamma=1.333;
     fEbindeK=0.008;
     float cg=1., cK=1.3e-4, cp=3.4e-5;
     p=GetRandom()*(cg+cK+cp);
     if (p<=cg){
        particle(1,fEgamma,fEgamma,0.,pi,0.,twopi,fTclev,fThlev);
        npg1333=fNbPart;
     } 
     else if (p<=cg+cK){
        particle(3,fEgamma-fEbindeK,fEgamma-fEbindeK,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,fEbindeK,fEbindeK,0.,pi,0.,twopi,0.,0.); 
     }
     else{
           pair(fEgamma-1.022);
        }
  }
  if (npg1333!=0 && npg1173!=0){
     float p1333=sqrt(pow(fPmoment[0][npg1333],2)+pow(fPmoment[1][npg1333],2)                     +pow(fPmoment[2][npg1333],2));
     float p1173=sqrt(pow(fPmoment[0][npg1173],2)+pow(fPmoment[1][npg1173],2)
                     +pow(fPmoment[2][npg1173],2)); 
     // Coefficients in formula 1+a2*ctet**2+a4*ctet**4 are from:
     // R.D.Evans, "The Atomic Nucleus", Krieger Publ. Comp., 1985, 
     // p. 240 (4(2)2(2)0 cascade).
     // They correspond to coefficients in formula 
     // 1+a2*p2+a4*p4, a2=0.1020, a4=0.0091 
     // in K.Siegbahn, "Alpha-, Beta- and Gamma-Ray Spectroscopy", 
     // North-Holland Publ. Comp., 1968, p. 1033.
     float ctet,stet1,stet2,ctet1,ctet2;
     float a2=1./8.;
     float a4=1./24.;
     do{
           phi1=twopi*GetRandom();
           ctet1=1.-2.*GetRandom();
           stet1=sqrt(1.-ctet1*ctet1);
           phi2=twopi*GetRandom();
           ctet2=1.-2.*GetRandom();
  	   stet2=sqrt(1.-ctet2*ctet2);
  	   ctet=ctet1*ctet2+stet1*stet2*cos(phi1-phi2);
     } while(GetRandom()*(1.+abs(a2)+abs(a4)) > 1.+a2*ctet*ctet+a4*pow(ctet,4)); 
     fPmoment[0][npg1333]=p1333*stet1*cos(phi1);
     fPmoment[1][npg1333]=p1333*stet1*sin(phi1);
     fPmoment[2][npg1333]=p1333*ctet1;
     fPmoment[0][npg1173]=p1173*stet2*cos(phi2);
     fPmoment[1][npg1173]=p1173*stet2*sin(phi2);
     fPmoment[2][npg1173]=p1173*ctet2;
  }
}
//-----------------------------------------------

void Decay::Cs136(float tcnuc)
{
  // Model for scheme of Cs136 decay ("Table of Isotopes", 8th ed., 
  // 1996 and Nucl. Data Sheets 95(2002)837).
  // VIT, 11.09.2002.
  fThnuc=1.137024e+6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float pbeta, p;
  bool next1867=false, next2054=false, next2140=false, next2207=false;
  bool next8190=false, next2031=false;  
  pbeta=100.*GetRandom();
  if (pbeta<=2.025){       // 2.025%
     beta(0.174,0.,0.);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=48.02){
        nucltransK(0.507,0.037,1.1e-2,0.);
        next1867=true;
     }
     else if (p<=73.36){
        nucltransK(0.320,0.037,3.5e-2,0.);
        next2054=true;
     }
     else if (p<=77.36){
        nucltransK(0.234,0.037,2.0e-2,0.);
        next2140=true;
     }
     else{
        nucltransK(0.167,0.037,2.5e-1,0.);
        next2207=true;
     }
  }
  else if (pbeta<= 2.233){ // 0.208%
     beta(0.191,0.,0.);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=47.12){
        nucltransK(1.538,0.037,6.0e-4,0.2e-4);
        next8190=true;
     }
     else if (p<=85.24){
        nucltransK(0.490,0.037,1.2e-2,0.);
        next1867=true;
     }
     else{
        nucltransK(0.303,0.037,4.5e-2,0.);
        next2054=true;
     }
  }
  else if (pbeta<=72.009){ // 69.776%
     beta(0.341,0.,0.); 
     next2207=true;
  }
  else if (pbeta<=82.431){ // 10.422% 
     beta(0.408,0.,0.);
     next2140=true;
  }
  else if (pbeta<=87.097){ // 4.666%
     beta(0.494,0.,0.);
     next2054=true;
  }
  else{                    // 12.903%
     beta(0.681,0.,0.);
     next1867=true;
  }
  
  if (next2207){
     fThlev=3.11e-9;
     p=100.*GetRandom();
     if (p<=61.84){
        nucltransK(0.341,0.037,3.0e-2,0.);
        next1867=true;
     }
     else if (p<=76.78){
        nucltransK(0.176,0.037,5.0e-2,0.);
        next2031=true;
     }
     else if (p<=88.48){
        nucltransK(0.153,0.037,4.3e-1,0.);
        next2054=true;
     }
     else{
        nucltransK(0.067,0.037,6.9e-1,0.);
        next2140=true;
     }
  }
  if (next2140){
     fThlev=1.6e-9;
     p=100.*GetRandom();
     if (p<=0.27){
        nucltransK(1.321,0.037,1.7e-3,0.1e-4);
        next8190=true; 
     }
     else if (p<=60.13){
        nucltransK(0.273,0.037,1.6e-2,0.);
        next1867=true;
     }
     else if (p<=62.88){
        nucltransK(0.109,0.037,1.47,0.);
        next2031=true;
     }
     else{
        nucltransK(0.086,0.037,3.5e-1,0.);
        next2054=true;
     }
  }
  if (next2054){
     fThlev=0.; 
     p=100.*GetRandom();
     if (p<=97.90){
        nucltransK(1.235,0.037,1.0e-3,0.1e-4);
        next8190=true;
     }
     else{
        nucltransK(0.187,0.037,1.9e-1,0.);
        next1867=true;
     }
  }
  if (next2031){
     fThlev=0.3084;
     nucltransK(0.164,0.037,2.26,0.);
     next1867=true; 
  }
  if (next1867){
     fThlev=0.;
     nucltransK(1.048,0.037,1.5e-3,0.);
     next8190=true;
  }
  if (next8190){
     fThlev=1.930e-12;
     nucltransK(0.819,0.037,2.7e-3,0.);
     return;
  }
}
//-----------------------------------------------
void Decay::Eu147(float tcnuc)
{
  // Scheme of 147Eu decay ("Table of Isotopes", 7th ed., 1978).
  // VIT, 10.03.1996.
  fThnuc=1.9008e6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float pdecay,p;
  bool next7990=false, next1970=false, next1210=false;
  pdecay=100.*GetRandom();
  if (pdecay<=0.002){
     particle(47,2.908,2.908,0.,pi,0.,twopi,0,0);
     return; 
  }
  else if (pdecay<= 0.202){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     nucltransK(1.449,0.047,5.0e-4,0.8e-4);
     next1970=true; 
  }
  else if (pdecay<= 0.502){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=31.){
        nucltransK(1.427,0.047,5.0e-4,0.7e-4);
        next1210=true;
     }
     else{
        nucltransK(0.750,0.047,1.8e-3,0.);
        next7990=true;
     }
  }
  else if (pdecay<= 0.552){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     nucltransK(1.350,0.047,5.5e-4,0.6e-4);
     next1210=true;
  }
  else if (pdecay<= 1.552){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=24.){
        nucltransK(1.332,0.047,5.5e-4,0.5e-4);
        next1210=true;
     }
     else{
        nucltransK(1.256,0.047,6.0e-4,0.4e-4);
        next1970=true;
     }
  }
  else if (pdecay<= 2.052){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=22.8){
        nucltransK(1.318,0.047,5.5e-4,0.5e-4);
       return;
     } 
     else if (p<=68.3){
       nucltransK(1.197,0.047,6.5e-4,0.3e-4);
       next1210=true;
     }
     else{
       nucltransK(1.120,0.047,7.5e-4,0.2e-4);
       next1970=true;
     }
  }
  else if (pdecay<= 2.252){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=15.){
        nucltransK(1.107,0.047,7.5e-4,0.1e-4);
        next1970=true;  
     }
     else{
        nucltransK(0.505,0.047,4.0e-3,0.);
        next7990=true;
     }
  }
  else if (pdecay<=11.352){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=62.5){
        nucltransK(1.077,0.047,4.3e-3,0.1e-4);
        return;
     }
     else if (p<=98.5){
        nucltransK(0.956,0.047,5.5e-3,0.);
        next1210=true;
     }
     else{
        nucltransK(0.880,0.047,9.5e-3,0.);
        next1970=true;
     }
  }
  else if (pdecay<=11.582){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=40.){
        nucltransK(1.063,0.047,1.1e-3,0.1e-4);
        return;
     }
     else{
        nucltransK(0.942,0.047,1.3e-3,0.);
        next1210=true;
     }
  }
  else if (pdecay<=17.102){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=57.){
        nucltransK(0.933,0.047,1.3e-3,0.);
        next1210=true;
     } 
     else{
        nucltransK(0.857,0.047,1.6e-3,0.);
        next1970=true; 
     }
  }
  else if (pdecay<=17.182){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=45.){
        nucltransK(0.886,0.047,1.5e-3,0.);
        next1210=true;
     }
     else{
        nucltransK(0.809,0.047,9.0e-3,0.);
        next1970=true;
     }
  }
  else if (pdecay<=35.382){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     next7990=true;
  }
  else if (pdecay<=57.592){
     p=100.*GetRandom();
     if (p<=99.4) particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     else beta(0.511,0.,0.); 
     next1970=true;
  }
  else if (pdecay<=76.792){
     p=100.*GetRandom();
     if (p<=99.5) particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     else beta(0.587,0.,0.);
     next1210=true;
  }
  else{
     p=100.*GetRandom();
     if (p<=99.3) particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     else beta(0.708,0.,0.);
  }
  if (next7990){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=24.2){
         nucltransK(0.799,0.047,4.5e-3,0.);
         return;
     }
     else if (p<=71.7){
         nucltransK(0.678,0.047,1.4e-2,0.);
         next1210=true;
     }
     else{
         nucltransK(0.601,0.047,1.8e-2,0.);
         next1970=true;
     }
  }
  if (next1970){
     fThlev=1.3e-9;
     p=100.*GetRandom();
     if (p<=97.2){
        nucltransK(0.197,0.047,2.2e-1,0.);
        return; 
     }
     else {
        nucltransK(0.076,0.047,3.5e-0,0.);
        next1210=true;
     }
  }
  if (next1210){
     fThlev=0.77e-9;
     nucltransK(0.121,0.047,1.1e-0,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Eu152(float tcnuc)
{
  // Scheme of 152Eu decay("Table of Isotopes", 8th ed., 1996 
  // and NDS 79(1996)1).
  // VIT, 16.11.1996.
  // VIT, approximate electron capture from K, L or M shells, 27.12.2006.
  // VIT, correction to the 1fnu shape for Qbeta=1.474, 13.11.2007.
  fThnuc=4.273413e8;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  bool next1530=false, next1372=false, next1293=false, next1234=false;
  bool next1222=false, next1086=false, next1041=false, next1023=false;
  bool next963=false,  next810=false,  next707=false,  next685=false;
  bool next366=false,  next122=false;
  bool next1434=false, next1318=false, next1123=false, next1109=false;
  bool next1048=false, next931=false, next755=false, next615=false, next344=false;
  float pdecay, pbeta, pKLM, pec, p;
  pdecay=100.*GetRandom();
  if (pdecay<=72.08){  //  EC to 152Sm
     //approximate electron capture from K (82%), L (14%) or M (4%) shell
     pKLM=100.*GetRandom();
     if (pKLM<=82.) particle(1,0.047,0.047,0.,pi,0.,twopi,0,0);
     if (pKLM>82.&& pKLM<=96.) particle(1,0.008,0.008,0.,pi,0.,twopi,0,0);
     if (pKLM>96.)  particle(1,0.001,0.001,0.,pi,0.,twopi,0,0);
     pec=100.*GetRandom();
     if (pec<=0.071){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=15.21){
           nucltransK(1.769,0.047,8.0e-4,0.9e-4);
           return;
        }        
        else if (p<=25.62){
           nucltransK(1.647,0.047,4.0e-4,0.5e-4);
           next122=true;
        }
        else if (p<=65.65){
           nucltransK(0.959,0.047,2.5e-3,0.);
           next810=true;
        }
        else if (p<=88.07){
           nucltransK(0.806,0.047,1.4e-3,0.);
           next963=true;
        }
        else if (p<=99.28){
           nucltransK(0.536,0.047,1.0e-2,0.);
           next1234=true; 
        }
        else{
           nucltransK(0.239,0.047,1.0e-1,0.);
           next1530=true;
        }
     }
     else if (pec<= 0.118){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=0.43){
           nucltransK(1.635,0.047,3.5e-4,0.5e-4);
           next122=true;
        }
        else if (p<=14.38){
           nucltransK(1.390,0.047,5.0e-4,0.2e-4);
           next366=true;
        }
        else if (p<=81.12){
           nucltransK(0.671,0.047,1.1e-2,0.);
           next1086=true;
        }
        else if (p<=85.46){
           nucltransK(0.523,0.047,3.5e-3,0.);
           next1234=true;
        }
        else{
           nucltransK(0.386,0.047,7.0e-3,0.);
           next1372=true;
        }
     }
     else if (pec<= 0.180){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=11.26){
           nucltransK(1.608,0.047,4.0e-4,0.8e-4);
           next122=true;
        } 
        else if (p<=66.17){
           nucltransK(1.364,0.047,5.3e-4,0.3e-4);
           next366=true;           
        }
        else if (p<=79.35){
           nucltransK(0.644,0.047,2.2e-3,0.);
           next1086=true;
        }
        else{
           nucltransK(0.496,0.047,4.0e-3,0.);
           next1234=true;
        } 
     }
     else if (pec<= 1.416){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=30.00){
           nucltransK(0.769,0.047,8.0e-3,0.);
           next122=true;
        } 
        else if (p<=31.85){
           nucltransK(0.769,0.047,8.0e-3,0.);
           next810=true;
        }
        else if (p<=34.02){
           nucltransK(0.769,0.047,8.0e-3,0.);
           next963=true;
        }
        else if (p<=86.93){
           nucltransK(0.769,0.047,8.0e-3,0.);
           next1086=true;
        }
        else if (p<=99.47){
           nucltransK(0.769,0.047,8.0e-3,0.);
           next1234=true;
        }
        else{
           nucltransK(0.769,0.047,8.0e-3,0.);
           next1293=true;
        }
     }
     else if (pec<= 1.445){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=71.28){
           nucltransK(0.906,0.047,2.8e-3,0.);
           next707=true;
        } 
        else if (p<=94.09){
           nucltransK(0.572,0.047,2.8e-3,0.);
           next1041=true;
        }
        else{
           nucltransK(0.391,0.047,7.0e-3,0.);
           next1222=true;
        }
     }
     else if (pec<= 4.292){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=24.04){
           nucltransK(1.458,0.047,4.5e-4,0.4e-4);
           next122=true; 
        } 
        else if (p<=92.32){ 
           nucltransK(1.213,0.047,6.2e-4,0.1e-4);
           next366=true;
        }
        else if (p<=96.69){
           nucltransK(0.769,0.047,1.5e-3,0.);
           next810=true;
        }
        else if (p<=97.13){
           nucltransK(0.616,0.047,7.0e-3,0.);
           next963=true;
        }
        else if (p<=98.06){
           nucltransK(0.557,0.047,3.0e-3,0.);
           next1023=true;
        }
        else if (p<=98.26){
           nucltransK(0.538,0.047,9.5e-3,0.);
           next1041=true;
        }
        else if (p<=99.79){
           nucltransK(0.494,0.047,4.0e-3,0.);
           next1086=true;
        }
        else{
           nucltransK(0.208,0.047,3.5e-2,0.);
           next1372=true;
        }
     }
     else if (pec<=38.592){
        next1530=true;
     }
     else if (pec<=39.883){
        next1372=true;
     }
     else if (pec<=40.744){
        next1293=true;
     }
     else if (pec<=64.629){
        next1234=true;
     }
     else if (pec<=94.069){
        next1086=true;
     }
     else if (pec<=94.152){
        next1041=true;
     }
     else if (pec<=94.471){
        next1023=true;
     }
     else if (pec<=96.179){
        next810=true;
     }
     else if (pec<=97.359){
        next366=true;
     }
     else{
       next122=true;
     } 
     if (next1530){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=84.56){
           nucltransK(1.408,0.047,5.0e-4,0.3e-4);
           next122=true;
        }
        else if (p<=84.80){
           nucltransK(0.719,0.047,2.0e-3,0.);
           next810=true;
        }
        else if (p<=85.32){
           nucltransK(0.566,0.047,1.4e-2,0.);
           next963=true;           
        }
        else if (p<=86.97){
           nucltransK(0.489,0.047,1.4e-2,0.);
           next1041=true;
        }
        else if (p<=98.22){
           nucltransK(0.444,0.047,5.7e-3,0.);
           next1086=true;
        }
        else{
           nucltransK(0.296,0.047,1.5e-2,0.);
           next1234=true;
        }
     }
     if (next1372){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=19.62){
           nucltransK(1.250,0.047,1.6e-3,0.1e-4);
           next122=true;
        }
        else if (p<=88.70){
           nucltransK(1.005,0.047,2.6e-3,0.);
           next366=true;
        }
        else if (p<=90.70){
           nucltransK(0.665,0.047,6.3e-3,0.);
           next707=true;
        }
        else if (p<=90.81){
           nucltransK(0.561,0.047,3.5e-3,0.);
           next810=true;
        }
        else if (p<=98.82){
           nucltransK(0.331,0.047,1.2e-2,0.);
           next1041=true;
        }
        else{
           nucltransK(0.286,0.047,6.6e-2,0.);
           next1086=true;
        }
     }
     if (next1293){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=16.39){
           nucltransK(1.293,0.047,1.3e-3,0.1e-4);
           return;
        }
        else if (p<=22.12){
           nucltransK(1.171,0.047,1.6e-3,0.);
           next122=true;
        }
        else if (p<=64.58){
           nucltransK(0.926,0.047,3.0e-3,0.);
           next366=true;
        }
        else if (p<=68.91){
           nucltransK(0.482,0.047,2.0e-2,0.);
           next810=true;
        }
        else if (p<=88.66){
           nucltransK(0.329,0.047,1.2e-2,0.);
           next963=true;
        }
        else if (p<=89.98){
           nucltransK(0.270,0.047,7.9e-2,0.);
           next1023=true;
        }
        else{
           nucltransK(0.252,0.047,2.3e-2,0.);
           next1041=true; 
        }
     }
     if (next1234){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=76.33){
           nucltransK(1.112,0.047,2.0e-3,0.2e-4);
           next122=true; 
        }
        else if (p<=99.76){
           nucltransK(0.867,0.047,3.5e-3,0.);
           next366=true;
        }
        else if (p<=99.78){
           nucltransK(0.423,0.047,2.7e-2,0.);
           next810=true;
        }
        else{
           nucltransK(0.148,0.047,5.8e-1,0.);
           next1086=true;
        }
     }
     if (next1222){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=80.65){
           nucltransK(0.855,0.047,1.6e-3,0.);
           next366=true; 
        }
        else{
           nucltransK(0.515,0.047,3.7e-3,0.);
           next707=true;
        }
     }
     if (next1086){
        fThlev=0.85e-12;
        p=100.*GetRandom();
        if (p<=40.36){
           nucltransK(1.086,0.047,2.1e-3,0.);
           return;
        }
        else if (p<=98.77){
           nucltransK(0.964,0.047,2.7e-3,0.);
           next122=true;
        }
        else if (p<=99.86){
           nucltransK(0.719,0.047,5.2e-3,0.);
           next366=true;
        }
        else{
           nucltransK(0.275,0.047,1.0e-1,0.);
           next810=true;
        }
     }
     if (next1041){
        fThlev=0;
        p=100.*GetRandom();
        if (p<=71.23){
           nucltransK(0.919,0.047,1.2e-3,0.);
           next122=true;
        }
        else{
           nucltransK(0.675,0.047,2.3e-3,0.);
           next366=true;
        }
     }
     if (next1023){
        fThlev=6.7e-12;
        p=100.*GetRandom();
        if (p<=35.73){
           nucltransK(0.901,0.047,3.1e-3,0.);
           next122=true;
        }
        else if (p<=91.56){
           nucltransK(0.656,0.047,5.7e-2,0.);
           next366=true;
        } 
        else if (p<=92.40){
           nucltransK(0.316,0.047,4.8e-2,0.);
           next707=true;
        }
        else{
           nucltransK(0.213,0.047,1.7e-1,0.);
           next810=true;
        }
     }
     if (next963){
        fThlev=28.2e-15;
        p=100.*GetRandom();
        if (p<=45.143){
           nucltransK(0.963,0.047,1.1e-3,0.);
           return;
        }
        else if (p<=99.994){
           nucltransK(0.842,0.047,1.4e-3,0.);
           next122=true;
        }
        else if (p<=99.995){
           nucltransK(0.279,0.047,1.8e-2,0.);
           next685=true;
        }
        else{
           nucltransK(0.153,0.047,8.7e-2,0.);
           next810=true;
        }
     } 
     if (next810){
        fThlev=7.2e-12;
        p=100.*GetRandom();
        if (p<=21.66){
           nucltransK(0.810,0.047,4.0e-3,0.);
           return;
        }
        else if (p<=78.06){
           nucltransK(0.689,0.047,4.3e-2,0.);
           next122=true;
        } 
        else if (p<=99.21){
           nucltransK(0.444,0.047,1.8e-2,0.);
           next366=true;
        }
        else{
           nucltransK(0.126,0.047,1.0e-0,0.);
           next685=true;
        }
     }
     if (next707){
        fThlev=10.1e-12;
        nucltransK(0.340,0.047,3.8e-4,0.);
        next366=true; 
     }
     if (next685){
        fThlev=6.2e-12;
        p=100.*GetRandom();
        if (p<= 1.43){
           particle(3,0.638,0.638,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.047,0.047,0.,pi,0.,twopi,0,0);
           return;
        }
        else{
           nucltransK(0.563,0.047,9.5e-3,0.);
           next122=true;
        }
     }
     if (next366){
        fThlev=57.7e-12;
        nucltransK(0.245,0.047,1.1e-1,0.);
        next122=true;
     }
     if (next122){
        fThlev=1.428e-9;
        nucltransK(0.122,0.047,1.2e-0,0.);
        return;
     }
  }    
  else{               //   b- to 152Gd
      pbeta=100.*GetRandom();
      if (pbeta<=0.071){
         beta(0.126,0.,0.);   
         fThlev=0.;
         p=100.*GetRandom();
         if (p<=84.03){
            nucltransK(1.348,0.050,1.5e-3,0.2e-4);
            next344=true; 
         } 
         else{ 
            nucltransK(0.937,0.050,3.2e-3,0.);
            next755=true;
         }
      }
      else if (pbeta<= 6.421){
         beta(0.175,0.,0.); 
         fThlev=0.;
         p=100.*GetRandom();
         if (p<=89.21){
            nucltransK(1.299,0.050,7.5e-4,0.3e-4);
            next344=true;
         } 
         else if (p<=94.47){
            nucltransK(0.713,0.050,2.2e-3,0.);
            next931=true;
         }
         else if (p<=96.83){
            nucltransK(0.534,0.050,4.5e-3,0.);
            next1109=true;
         } 
         else if (p<=99.75){
            nucltransK(0.520,0.050,1.8e-2,0.);
            next1123=true; 
         }
         else if (p<=99.76){
            nucltransK(0.325,0.050,1.3e-2,0.);
            next1318=true;
         }
         else{
            nucltransK(0.209,0.050,4.0e-2,0.);
            next1434=true;
         }
      }
      else if (pbeta<= 6.778){
         beta(0.212,0.,0.);
         fThlev=0.;
         p=100.*GetRandom();
         if (p<=7.64){
            nucltransK(1.606,0.050,1.2e-3,0.5e-4);
            return;
         }
         else if (p<=40.86){
            nucltransK(1.261,0.050,2.7e-3,0.5e-4);
            next344=true;
         }
         else if (p<=72.09){
            nucltransK(0.990,0.050,3.0e-3,0.);
            next615=true;
         }
         else if (p<=91.36){
            nucltransK(0.675,0.050,7.6e-3,0.);
            next931=true;
         }
         else if (p<=95.35){
            nucltransK(0.558,0.050,1.1e-2,0.);
            next1048=true;
         }
         else if (p<=99.99){
            nucltransK(0.496,0.050,9.7e-2,0.);
            next1109=true; 
         }
         else{
            nucltransK(0.482,0.050,5.0e-3,0.);
            next1123=true;
         }
      }
      else if (pbeta<=15.339){
         beta(0.384,0.,0.);  
         next1434=true;
      }
      else if (pbeta<=15.981){
         beta(0.500,0.,0.); 
         next1318=true;
      }
      else if (pbeta<=16.063){
         beta(0.536,0.,0.);
         fThlev=0.;
         p=100.*GetRandom();
         if (p<=58.41){
            nucltransK(0.527,0.050,9.4e-2,0.);
            next755=true; 
         }
         else if (p<=98.13){
            nucltransK(0.352,0.050,3.8e-2,0.);
            next931=true;
         }
         else{
            nucltransK(0.172,0.050,3.7e-1,0.);
            next1109=true;
         }
      }
      else if (pbeta<=65.291){
         beta(0.695,0.,0.);
         next1123=true;
      } 
      else if (pbeta<=66.325){
         beta(0.709,0.,0.);
         next1109=true;
      }
      else if (pbeta<=67.395){
         beta(0.887,0.,0.);
         next931=true;
      }
      else if (pbeta<=70.748){
         beta(1.063,0.,0.);
         next755=true;
      }
      else{
         beta1f(1.474,0.,0.,0.,-0.4026,0.0928,0.);
         next344=true;
      }
      if (next1434){
         fThlev=0.;
         p=100.*GetRandom();
         if (p<=71.07){
            nucltransK(1.090,0.050,2.3e-3,0.);
            next344=true;
         }
         else if (p<=90.54){
            nucltransK(0.679,0.050,6.9e-3,0.);
            next755=true;
         }
         else if (p<=96.87){
            nucltransK(0.503,0.050,1.4e-2,0.);
            next931=true;
         }
         else if (p<=99.99){
            nucltransK(0.325,0.050,6.3e-2,0.);
            next1109=true;
         }
         else{
            nucltransK(0.115,0.050,1.4e-0,0.);
            next1318=true; 
         }
      }
      if (next1318){
         fThlev=0.;
         p=100.*GetRandom();
         if (p<=6.87){
            nucltransK(1.318,0.050,1.6e-3,0.1e-4);
            return;
         }
         else if (p<=62.27){
            nucltransK(0.974,0.050,5.6e-3,0.);
            next344=true;
         }
         else if (p<=76.12){
            nucltransK(0.703,0.050,6.0e-3,0.);
            next615=true;
         } 
         else if (p<=77.39){
            nucltransK(0.563,0.050,1.0e-2,0.);
            next755=true;
         }
         else if (p<=84.76){
            nucltransK(0.388,0.050,4.5e-1,0.);
            next931=true;
         }
         else if (p<=92.35){
            nucltransK(0.271,0.050,8.4e-2,0.);
            next1048=true;
         }
         else if (p<=92.79){
            nucltransK(0.209,0.050,5.0e-1,0.);
            next1109=true; 
         }
         else{
            nucltransK(0.195,0.050,4.9e-2,0.);
            next1123=true;
         }
      }
      if (next1123){
         fThlev=7.3e-12;
         p=100.*GetRandom();
         if (p<=93.11){
            nucltransK(0.779,0.050,1.9e-3,0.);
            next344=true;
         }
         else if (p<=99.26){
            nucltransK(0.368,0.050,9.7e-3,0.);
            next755=true;
         }
         else{
            nucltransK(0.193,0.050,5.0e-2,0.);
            next931=true;
         } 
      } 
      if (next1109){
         fThlev=7.3e-12;
         p=100.*GetRandom();
         if (p<=50.00){
            nucltransK(1.109,0.050,2.2e-3,0.);
            return; 
         }
         else if (p<=97.50){
            nucltransK(0.765,0.050,5.2e-3,0.);
            next344=true;
         }
         else{
            nucltransK(0.494,0.050,1.5e-2,0.);
            next615=true;
         }
      }
      if (next1048){
         fThlev=0.;
         p=100.*GetRandom();
         if (p<=0.88){
            particle(3,0.998,0.998,0.,pi,0.,twopi,fTclev,fThlev);
            particle(1,0.050,0.050,0.,pi,0.,twopi,0,0);
            return;
         }
         else if (p<=65.79){
            nucltransK(0.703,0.050,6.0e-3,0.);
            next344=true;
         }
         else if (p<=83.77){
            particle(3,0.383,0.383,0.,pi,0.,twopi,fTclev,fThlev);
            particle(1,0.050,0.050,0.,pi,0.,twopi,0,0);
            next615=true;
         }
         else{
            nucltransK(0.117,0.050,1.4e-0,0.);
            next931=true;
         }
      }
      if (next931){
         fThlev=7.3e-12;
         p=100.*GetRandom();
         if (p<=12.40){
            nucltransK(0.931,0.050,3.2e-3,0.);
            return;
         }
         else if (p<=90.91){
            nucltransK(0.586,0.050,2.4e-2,0.);
            next344=true;
         }
         else if (p<=99.55){
            nucltransK(0.315,0.050,5.2e-2,0.);
            next615=true;
         }
         else{
            nucltransK(0.175,0.050,3.5e-1,0.);
            next755=true;
         }
      }
      if (next755){
         fThlev=7.3e-12;
         nucltransK(0.411,0.050,2.4e-2,0.);
         next344=true;
      }
      if (next615){
         fThlev=37.e-12;
         p=100.*GetRandom();
         if (p<=11.35){
            particle(3,0.565,0.565,0.,pi,0.,twopi,fTclev,fThlev);
            particle(1,0.050,0.050,0.,pi,0.,twopi,0,0);
            return;
         }
         else{
            nucltransK(0.271,0.050,8.3e-2,0.);
            next344=true;
         }
      }
      if (next344){
         fThlev=31.9e-12;
         nucltransK(0.344,0.050,4.0e-2,0.);
         return;
      }
  }
}
//-----------------------------------------------

void Decay::Eu154(float tcnuc)
{
  // Scheme of 154Eu decay ("Table of Isotopes", 8th ed., 1996 
  // and NDS 69(1993)507).
  // VIT, 13.11.1996.
  // VIT, correction to the 1fnu shape for Qbeta=1.846, 13.11.2007.
  fThnuc=2.711670e8;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  bool next1775=false, next1661=false, next1646=false, next1560=false;
  bool next1531=false, next1510=false, next1418=false, next1414=false;
  bool next1404=false, next1398=false, next1296=false, next1294=false;
  bool next1277=false, next1264=false, next1252=false, next1241=false; 
  bool next1182=false, next1136=false, next1128=false, next1048=false;
  bool next996=false,   next815=false,  next718=false,  next681=false;
  bool next371=false, next123=false;
  float pdecay,pbeta,p;
  pdecay=100.*GetRandom();
  if (pdecay<=0.020){     // 0.020% EC to 154Sm
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     if (pdecay<=0.005){
        fThlev=172.e-12;
        nucltransK(0.185,0.047,2.7e-1,0.);
        fThlev=3.02e-9;
        nucltransK(0.082,0.047,4.9+0,0.); 
        return;
     }
     else{
        fThlev=3.02e-9;
        nucltransK(0.082,0.047,4.9+0,0.);
        return;
     }
  }  
  else{                   // 99.980% b- to 154Gd
     pbeta=100.*GetRandom();
     if (pbeta<=0.019){
        beta(0.130,0.,0.);  
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=4.0){
           nucltransK(1.839,0.050,1.0e-3,1.0e-4);
           return;
        }
        else if (p<=6.8){
           nucltransK(1.717,0.050,1.1e-3,0.7e-4);
           next123=true;
        }
        else if (p<=38.1){
           nucltransK(1.023,0.050,3.8e-3,0.);
           next815=true;
        }
        else if (p<=87.8){
           nucltransK(0.790,0.050,5.0e-3,0.);
           next1048=true;
        }
        else if (p<=94.0){
           nucltransK(0.193,0.050,2.8e-1,0.);
           next1646=true;
        } 
        else{
           nucltransK(0.063,0.050,1.4e+1,0.);
           next1775=true;
        }
     }
     else if (pbeta<= 0.087){
        beta(0.172,0.,0.);
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=2.4){
           nucltransK(1.674,0.050,5.0e-4,0.9e-4);
           next123=true;
        } 
        else if (p<=4.0){
           nucltransK(1.426,0.050,6.0e-4,0.4e-4);
           next371=true;
        }
        else if (p<=15.3){
           nucltransK(0.982,0.050,1.2e-3,0.);
           next815=true;
        }
        else if (p<=58.7){
           nucltransK(0.801,0.050,1.8e-3,0.);
           next996=true;
        }
        else if (p<=77.0){
           nucltransK(0.669,0.050,2.5e-3,0.);
           next1128=true;
        }
        else if (p<=82.0){
           nucltransK(0.556,0.050,1.1e-2,0.);
           next1241=true;
        }
        else if (p<=91.6){
           nucltransK(0.533,0.050,4.0e-3,0.);
           next1264=true;
        }
        else if (p<=99.9){
           nucltransK(0.393,0.050,3.0e-2,0.);
           next1404=true;
        }
        else{
           nucltransK(0.266,0.050,2.0e-2,0.);
           next1531=true;
        }
     }
     else if (pbeta<=28.587){
        beta(0.249,0.,0.);
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=6.272){
           nucltransK(1.596,0.050,5.0e-4,0.7e-4);
           next123=true; 
        }
        else if (p<=9.404){
           nucltransK(0.904,0.050,1.4e-3,0.);
           next815=true;
        } 
        else if (p<=79.950){
           nucltransK(0.723,0.050,2.2e-3,0.);
           next996=true;
        }
        else if (p<=97.361){
           nucltransK(0.592,0.050,3.3e-3,0.);
           next1128=true;
        }
        else if (p<=98.150){
           nucltransK(0.478,0.050,1.6e-2,0.);
           next1241=true;
        }
        else if (p<=98.362){
           nucltransK(0.468,0.050,3.5e-2,0.);
           next1252=true;
        }
        else if (p<=98.594){
           nucltransK(0.322,0.050,9.0e-2,0.);
           next1398=true;
        }
        else if (p<=98.629){
           nucltransK(0.305,0.050,9.5e-2,0.);
           next1414=true;
        }
        else if (p<=98.665){
           nucltransK(0.301,0.050,1.6e-2,0.);
           next1418=true;
        }
        else if (p<=98.674){
           nucltransK(0.209,0.050,3.0e-1,0.);
           next1510=true;
        }
        else if (p<=99.513){
           nucltransK(0.188,0.050,5.5e-2,0.);
           next1531=true;
        }
        else if (p<=99.986){
           nucltransK(0.160,0.050,4.5e-1,0.);
           next1560=true;
        }
        else{
           nucltransK(0.058,0.050,1.2e-0,0.);
           next1661=true;
        }
     }
     else if (pbeta<=29.417){
        beta(0.308,0.,0.);
        next1661=true;
     }
     else if (pbeta<=29.564){
        beta(0.323,0.,0.);
        next1646=true;
     }
     else if (pbeta<=31.174){
        beta(0.352,0.,0.);
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=43.76){
           nucltransK(1.494,0.050,5.6e-4,0.5e-4);
           next123=true;
        }
        else if (p<=97.79){
           nucltransK(1.246,0.050,7.6e-4,0.2e-4);
           next371=true;
        }
        else if (p<=97.80){
           nucltransK(0.801,0.050,1.8e-3,0.);
           next815=true;
        }
        else if (p<=98.37){
           nucltransK(0.621,0.050,2.8e-3,0.);
           next996=true;
        }
        else if (p<=99.00){
           nucltransK(0.569,0.050,3.5e-3,0.);
           next1048=true;
        }
        else if (p<=99.43){
           nucltransK(0.488,0.050,5.0e-3,0.);
           next1128=true;
        }
        else if (p<=99.73){
           nucltransK(0.481,0.050,5.0e-3,0.);
           next1136=true;
        }
        else if (p<=99.85){
           nucltransK(0.375,0.050,4.0e-2,0.);
           next1241=true;
        }
        else{ 
           nucltransK(0.219,0.050,2.0e-1,0.);
           next1398=true;
        }
     }
     else if (pbeta<=31.271){
        beta(0.409,0.,0.);
        next1560=true;
     } 
     else if (pbeta<=31.571){
        beta(0.438,0.,0.);
        next1531=true;
     }
     else if (pbeta<=31.679){
        beta(0.551,0.,0.);
        next1418=true;
     }
     else if (pbeta<=67.879){
        beta(0.571,0.,0.);
        next1398=true;
     } 
     else if (pbeta<=68.599){
        beta(0.705,0.,0.);
        next1264=true;
     }
     else if (pbeta<=68.909){
        beta(0.717,0.,0.);
        next1252=true;
     }
     else if (pbeta<=85.610){
        beta(0.841,0.,0.);
        next1128=true;
     }
     else if (pbeta<=89.110){
        beta(0.973,0.,0.);
        next996=true;
     }
     else if (pbeta<=89.810){
        beta(1.154,0.,0.);
        next815=true;
     }
     else if (pbeta<=90.000){
        beta(1.598,0.,0.);
        next371=true;
     }
     else{
        beta1f(1.846,0.,0.,0.,-0.2280,0.04465,0.);
        next123=true;
     }
    
     if (next1775){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=19.6){
           nucltransK(1.776,0.050,1.0e-3,0.8e-4);
           return;
        }
        else if (p<=62.3){
           nucltransK(1.652,0.050,1.6e-3,0.6e-4);
           next123=true;
        }
        else if (p<=79.4){
           nucltransK(1.405,0.050,1.4e-3,0.2e-4);
           next371=true;
        }
        else if (p<=88.4){
           nucltransK(1.095,0.050,2.3e-3,0.);
           next681=true; 
        }
        else if (p<=92.2){
           nucltransK(0.960,0.050,3.0e-3,0.);
           next815=true;
        }
        else if (p<=98.7){
           nucltransK(0.728,0.050,6.0e-3,0.);
           next1048=true;
        }
        else{
           nucltransK(0.648,0.050,1.4e-2,0.);
           next1128=true;
        }
     }
     if (next1661){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=6.03){
           nucltransK(1.538,0.050,1.2e-3,0.4e-4);
           next123=true;
        }
        else if (p<=8.84){
           nucltransK(1.290,0.050,1.8e-3,0.1e-4);
           next371=true;
        }
        else if (p<=75.84){
           nucltransK(0.845,0.050,4.0e-3,0.);
           next815=true;
        }
        else if (p<=79.11){
           nucltransK(0.665,0.050,7.0e-3,0.);
           next996=true;
        }
        else if (p<=89.76){
           nucltransK(0.613,0.050,1.2e-2,0.);
           next1048=true;
        }
        else if (p<=94.78){
           nucltransK(0.533,0.050,1.2e-2,0.);
           next1128=true;
        } 
        else if (p<=95.18){
           nucltransK(0.420,0.050,1.4e-1,0.);
           next1241=true;
        }
        else if (p<=98.45){
           nucltransK(0.397,0.050,3.5e-2,0.);
           next1264=true;
        }
        else{
           nucltransK(0.130,0.050,1.0e-0,0.);
           next1531=true; 
        }
     }
     if (next1646){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=1.26){
           nucltransK(1.523,0.050,1.2e-3,0.4e-4);
           next123=true;
        }
        else if (p<=2.57){
           nucltransK(1.275,0.050,1.8e-3,0.1e-4);
           next371=true;
        }
        else if (p<=3.97){
           nucltransK(0.928,0.050,3.5e-3,0.);
           next718=true;
        }
        else if (p<=6.77){
           nucltransK(0.830,0.050,5.7e-3,0.);
           next815=true;
        }
        else if (p<=51.79){
           nucltransK(0.650,0.050,7.3e-3,0.);
           next996=true;
        }
        else if (p<=56.52){
           nucltransK(0.598,0.050,1.4e-2,0.);
           next1048=true;
        }
        else if (p<=82.63){
           nucltransK(0.518,0.050,1.3e-2,0.);
           next1128=true;
        }
        else if (p<=84.79){
           nucltransK(0.394,0.050,8.5e-3,0.);
           next1252=true;
        }
        else if (p<=89.70){
           nucltransK(0.382,0.050,3.4e-2,0.);
           next1264=true;
        } 
        else if (p<=91.41){
           nucltransK(0.368,0.050,9.0e-3,0.);
           next1277=true;
        }
        else if (p<=95.37){
           nucltransK(0.352,0.050,3.5e-2,0.);
           next1294=true;
        }
        else if (p<=99.02){
           nucltransK(0.242,0.050,3.5e-2,0.);
           next1404=true;
        }
        else{
           nucltransK(0.228,0.050,1.4e-1,0.);
           next1418=true;
        }
     }
     if (next1560){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=90.61){
           nucltransK(1.189,0.050,1.1e-3,0.1e-4);
           next371=true;
        }
        else if (p<=90.62){
           nucltransK(0.563,0.050,6.0e-2,0.);
           next996=true;
        }
        else if (p<=92.94){
           nucltransK(0.296,0.050,1.6e-2,0.);
           next1264=true; 
        }
        else if (p<=98.83){
           nucltransK(0.283,0.050,1.7e-2,0.);
           next1277=true;
        }
        else{
           nucltransK(0.163,0.050,5.0e-1,0.);
           next1398=true;
        }
     }
     if (next1531){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=0.99){
           nucltransK(1.531,0.050,1.1e-3,0.4e-4);
           return;
        }
        else if (p<=4.80){
           nucltransK(1.408,0.050,3.7e-3,0.2e-4);
           next123=true;
        }
        else if (p<=12.06){
           nucltransK(1.161,0.050,2.2e-3,0.1e-4);
           next371=true;
        }
        else if (p<=52.17){
           nucltransK(0.851,0.050,3.9e-3,0.);
           next681=true;
        }
        else if (p<=82.25){
           nucltransK(0.716,0.050,1.3e-2,0.);
           next815=true;
        }
        else if (p<=90.87){
           nucltransK(0.535,0.050,2.5e-2,0.);
           next996=true;
        }
        else if (p<=91.69){
           nucltransK(0.484,0.050,1.6e-2,0.);
           next1048=true; 
        }
        else if (p<=95.70){
           nucltransK(0.404,0.050,2.8e-2,0.);
           next1128=true;
        }
        else if (p<=96.26){
           nucltransK(0.290,0.050,1.8e-2,0.);
           next1241=true; 
        }
        else if (p<=96.75){
           nucltransK(0.280,0.050,1.9e-2,0.);
           next1252=true;
        }
        else if (p<=98.55){
           nucltransK(0.267,0.050,8.0e-2,0.);
           next1264=true;
        }
        else if (p<=99.59){
           nucltransK(0.238,0.050,1.2e-1,0.);
           next1294=true;
        }
        else{
           nucltransK(0.117,0.050,2.0e-1,0.);
           next1414=true;
        }
     }
     if (next1510){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=20.){
           nucltransK(1.510,0.050,6.0e-4,0.5e-4);
           return;
        }
        else{
           nucltransK(1.387,0.050,6.5e-4,0.3e-4);
           next123=true;
        }
     }
     if (next1418){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=11.51){
           nucltransK(1.418,0.050,1.4e-3,0.2e-4);
           return;
        }
        else if (p<=20.34){
           nucltransK(1.295,0.050,1.6e-3,0.1e-4);
           next123=true;
        }
        else if (p<=68.31){
           nucltransK(1.047,0.050,2.5e-3,0.);
           next371=true;
        } 
        else if (p<=69.27){
           nucltransK(0.737,0.050,5.5e-3,0.);
           next681=true;
        }
        else if (p<=90.86){
           nucltransK(0.603,0.050,3.8e-2,0.);
           next815=true;
        }
        else if (p<=91.87){
           nucltransK(0.422,0.050,1.4e-1,0.);
           next996=true;
        }
        else if (p<=95.71){
           nucltransK(0.371,0.050,3.2e-2,0.);
           next1048=true;
        }
        else if (p<=96.56){
           nucltransK(0.290,0.050,6.5e-2,0.);
           next1128=true;
        }
        else if (p<=97.71){
           nucltransK(0.236,0.050,1.2e-1,0.);
           next1182=true;
        }
        else if (p<=98.30){
           nucltransK(0.177,0.050,6.0e-2,0.);
           next1241=true;
        }
        else if (p<=98.94){
           nucltransK(0.167,0.050,7.0e-2,0.);
           next1252=true;
        }
        else if (p<=99.45){
           nucltransK(0.124,0.050,1.2e-0,0.);
           next1294=true;
        }
        else{
           nucltransK(0.123,0.050,1.2e-0,0.);
           next1296=true;
        }
     }
     if (next1414){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=21.22){
           nucltransK(1.414,0.050,6.1e-4,0.4e-4);
           return; 
        } 
        else if (p<=97.01){
           nucltransK(1.291,0.050,7.2e-4,0.2e-4);
           next123=true;
        }
        else if (p<=99.68){
           nucltransK(0.599,0.050,3.5e-3,0.);
           next815=true;
        }
        else if (p<=99.79){
           nucltransK(0.167,0.050,4.0e-1,0.);
           next1252=true;
        }
        else{
           nucltransK(0.120,0.050,1.8e-1,0.);
           next1294=true;
        }
     }
     if (next1404){
        fThlev=0.;
        nucltransK(1.033,0.050,1.1e-3,0.);
        next371=true;
     } 
     if (next1398){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=96.877){
           nucltransK(1.274,0.050,7.4e-4,0.2e-4);
           next123=true;
        }
        else if (p<=99.338){
           nucltransK(0.582,0.050,3.4e-3,0.);
           next815=true;
        }
        else if (p<=99.876){
           nucltransK(0.401,0.050,7.0e-2,0.);
           next996=true;
        }
        else if (p<=99.896){
           nucltransK(0.270,0.050,2.2e-2,0.);
           next1128=true;
        }
        else if (p<=99.902){
           nucltransK(0.260,0.050,2.3e-2,0.);
           next1136=true; 
        }
        else if (p<=99.929){
           nucltransK(0.156,0.050,9.0e-2,0.);
           next1241=true;
        }
        else{
           nucltransK(0.146,0.050,9.5e-2,0.);
           next1252=true;
        }
     }
     if (next1296){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=0.73){
           particle(3,1.245,1.245,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.050,0.050,0.,pi,0.,twopi,0,0);
           return;
        }
        else if (p<=73.72){
           nucltransK(1.173,0.050,2.2e-3,0.);
           next123=true;
        }
        else if (p<=74.45){
           particle(3,0.565,0.565,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.050,0.050,0.,pi,0.,twopi,0,0);
           next681=true;
        }
        else if (p<=78.83){
           nucltransK(0.480,0.050,1.6e-2,0.);
           next815=true;
        }
        else{
           nucltransK(0.299,0.050,6.0e-2,0.);
           next996=true;
        }
     }
     if (next1294){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=0.51){
           nucltransK(1.294,0.050,1.8e-3,0.1e-4);
           return;
        } 
        else if (p<=47.79){
           nucltransK(1.171,0.050,2.2e-3,0.1e-4);
           next123=true;
        }
        else if (p<=83.25){
           nucltransK(0.923,0.050,3.5e-3,0.);
           next371=true;
        }
        else{
           nucltransK(0.112,0.050,1.7e-0,0.);
           next1182=true;
        }
     }
     if (next1277){
        fThlev=0.;
        nucltransK(0.906,0.050,1.4e-3,0.);
        next371=true;
     }
     if (next1264){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=30.09){
           nucltransK(1.141,0.050,2.1e-3,0.);
           next123=true;
        }
        else if (p<=95.92){
           nucltransK(0.893,0.050,3.7e-3,0.);
           next371=true;
        }
        else if (p<=97.76){
           nucltransK(0.546,0.050,1.2e-2,0.);
           next718=true;
        }
        else{
           nucltransK(0.267,0.050,9.5e-2,0.);
           next996=true;
        }
     } 
     if (next1252){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=3.76){
           nucltransK(1.252,0.050,3.5e-3,0.);
           return;
        }
        else if (p<=80.51){
           nucltransK(1.129,0.050,9.1e-4,0.1e-4);
           next123=true;
        }
        else{
           nucltransK(0.881,0.050,1.5e-3,0.);
           next371=true;
        }
     }
     if (next1241){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=54.62){
           nucltransK(1.241,0.050,7.7e-4,0.2e-4);
           return;
        }
        else if (p<=99.02){
           nucltransK(1.118,0.050,9.3e-4,0.1e-4);
           next123=true; 
        }
        else if (p<=99.68){
           nucltransK(0.561,0.050,4.0e-3,0.);
           next681=true;
        }
        else{
           nucltransK(0.426,0.050,7.0e-3,0.);
           next815=true;
        }
     }
     if (next1182){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=0.21){
           particle(3,1.132,1.132,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.050,0.050,0.,pi,0.,twopi,0,0);
           return;
        } 
        else if (p<=84.00){
           nucltransK(1.059,0.050,2.5e-3,0.);
           next123=true;
        }
        else if (p<=84.84){
           particle(3,0.451,0.451,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.050,0.050,0.,pi,0.,twopi,0,0);
           next681=true;
        }
        else{
           nucltransK(0.367,0.050,3.3e-2,0.);
           next815=true;
        }
     }
     if (next1136){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=71.94){
           nucltransK(1.136,0.050,2.2e-3,0.);
           return;
        }
        else{
           nucltransK(1.013,0.050,2.8e-3,0.);
           next123=true;
        }
     }
     if (next1128){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=79.629){
           nucltransK(1.005,0.050,2.8e-3,0.);
           next123=true;
        }
        else if (p<=99.855){
           nucltransK(0.757,0.050,5.2e-3,0.);
           next371=true; 
        } 
        else if (p<=99.938){
           nucltransK(0.312,0.050,5.5e-2,0.);
           next815=true;
        }
        else if (p<=99.987){
           nucltransK(0.132,0.050,9.5e-1,0.);
           next996=true;
        }
        else{
           nucltransK(0.080,0.050,6.0e-0,0.);
           next1048=true;
        }
     }
     if (next1048){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=23.98){
           nucltransK(0.925,0.050,3.3e-3,0.);
           next123=true;
        }  
        else if (p<=86.75){
           nucltransK(0.677,0.050,5.1e-2,0.);
           next371=true;
        }
        else if (p<=90.58){
           nucltransK(0.330,0.050,4.5e-2,0.);
           next718=true;
        }
        else{
           nucltransK(0.232,0.050,1.4e-1,0.);
           next815=true;
        }
     }
     if (next996){
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=46.43){
           nucltransK(0.996,0.050,2.8e-3,0.);
           return;
        }
        else if (p<=98.59){
           nucltransK(0.873,0.050,3.7e-3,0.);
           next123=true;
        }
        else if (p<=99.95){
           nucltransK(0.625,0.050,8.0e-3,0.);
           next371=true;
        }
        else if (p<=99.98){
           nucltransK(0.315,0.050,5.2e-2,0.);
           next681=true;
        }
        else{
           nucltransK(0.181,0.050,3.5e-1,0.);
           next815=true;
        }
     }
     if (next815){
        fThlev=6.9e-12;
        p=100.*GetRandom();
        if (p<=17.86){
           nucltransK(0.816,0.050,4.3e-3,0.);
           return;
         }
         else if (p<=80.18){
           nucltransK(0.692,0.050,4.6e-2,0.);
           next123=true;
         }
         else if (p<=99.75){
           nucltransK(0.444,0.050,1.9e-2,0.);
           next371=true;
         } 
         else{
           nucltransK(0.135,0.050,8.7e-1,0.);
           next681=true;
         }
     }
     if (next718){
        fThlev=7.8e-12;
        nucltransK(0.347,0.050,3.9e-2,0.);
        next371=true;
     }
     if (next681){
        fThlev=4.0e-12;
        p=100.*GetRandom();
        if (p<=2.06){
           particle(3,0.631,0.631,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.050,0.050,0.,pi,0.,twopi,0,0);
           return;
        }
        else{
           nucltransK(0.558,0.050,1.1e-2,0.);
           next123=true;
        }
     }
     if (next371){
        fThlev=45.2e-12;
        nucltransK(0.248,0.050,1.1e-1,0.);
        next123=true;
     }
     if (next123){
        fThlev=1.186e-9;
        nucltransK(0.123,0.050,1.2e-0,0.);
        return;
     }
  }
}
//-----------------------------------------------

void Decay::Gd146(float tcnuc)
{
  // Scheme of 146Gd decay ("Table of Isotopes", 7th ed., 1978).
  // VIT, 10.03.1996.
  fThnuc=4.173e6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float pdecay,p;
  bool next385=false, next230=false, next115=false; 
  pdecay=100.*GetRandom();
  if (pdecay<= 0.3){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=36.){
        nucltransK(0.576,0.049,1.8e-2,0.);
        next230=true;
     } 
     else{
        nucltransK(0.421,0.049,4.5e-2,0.);
        next385=true;
     }
  }
  else if (pdecay<= 0.8){
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     fThlev=0.;
     nucltransK(0.743,0.049,1.0e-2,0.);
     return;     
  }
  else if (pdecay<= 77.2){
     p=100.*GetRandom();
     if (p<=99.91) particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     else beta(0.343,0.,0.); 
     next385=true;
  }
  else{
     particle(1,0.049,0.049,0.,pi,0.,twopi,0,0);
     next230=true;
  }
  if (next385){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=0.4){
        nucltransK(0.269,0.049,8.0e-2,0.);
        next115=true;
     }
     else{
        nucltransK(0.155,0.049,6.5e-1,0.);
        next230=true;
     }
  }
  if (next230){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=0.2){
        nucltransK(0.230,0.049,1.4e-1,0.);
        return;
     }     
     else{
        nucltransK(0.115,0.049,1.5e-0,0.);
        next115=true;
     }
  }
  if (next115){
     fThlev=0.;
     nucltransK(0.115,0.049,1.5e-0,0.);
     return;
  }  
} 
//-----------------------------------------------

void Decay::I126(float tcnuc)
{
  // Model for scheme of I126 decay(J.Katakura et al.,Nucl.Data Sheets  // 97(2002)765).
  // VIT, 25.11.2003. 
  fThnuc=1.1172e6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  float pdecay,pbeta,p,pec;
  bool next666=false, next389=false;
  pdecay=100.*GetRandom();
  if (pdecay<=47.300){     // 47.300% beta- to 126Xe
     fZdtr=54;
     pbeta=100.*GetRandom();
     if (pbeta<=7.65){     // 7.65%
        beta(0.378,0.,0.);
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=20.51){
           nucltransK(0.880,0.035,2.2e-3,0.);
           return;
        } 
        else{
           nucltransK(0.491,0.035,9.5e-3,0.);
           next389=true; 
        }
     }
     else if (pbeta<=78.22){// 70.57%
        beta(0.869,0.,0.);
        next389=true;
     }
     else{
        beta(1.258,0.,0.);
        return;
     }
     if (next389){
        fThlev=41.3e-12;
        nucltransK(0.389,0.035,1.9e-2,0.);
        return; 
     } 
  } 
  else if(pdecay<=98.992){ // 51.692% EC    to 126Te 
     particle(1,0.032,0.032,0.,pi,0.,twopi,0,0);
     pec=100.*GetRandom(); 
     if (pec<=0.014){      // 0.014%
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=66.00){
           nucltransK(2.045,0.032,3.0e-4,2.3e-4);
           return;
        }
        else{
          nucltransK(1.379,0.032,1.2e-3,0.7e-4);
          next666=true;
        }
     }
     else if (pec<= 0.015){//  0.001%
        fThlev=0.;
        nucltransK(1.207,0.032,1.1e-3,0.1e-4);
        next666=true;
     }
     else if (pec<= 8.630){//  8.615%
        fThlev=0.;
        p=100.*GetRandom();
        if (p<=6.83){
           nucltransK(1.420,0.032,7.0e-4,0.2e-4);
           return;
        }
        else{
           nucltransK(0.754,0.032,2.8e-3,0.);
           next666=true;
        }
     }
     else if (pec<=63.800){// 55.170%
        next666=true; 
     }
     else return;          // 36.200%

     if (next666){
        fThlev=0.;
        nucltransK(0.666,0.032,3.8e-3,0.);
        return;
     }
  }
  else{                    // 1.008% beta+ to 126Te
     pbeta=100.*GetRandom();
     fZdtr=-52;
     if (pbeta<=19.64){
        beta(0.467,0.,0.);
        fThlev=0.;
        nucltransK(0.666,0.032,3.8e-3,0.);
        return;
     }
     else{
        beta(1.133,0.,0.);
        return;
     }
  }
}
//-----------------------------------------------

void Decay::I133(float tcnuc)
{
  // Scheme of I133 decay in accordance with S.Rab, NDS 75(1995)491.
  // VIT, 13.12.2003.
  fThnuc=74880.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  bool next1386=false, next1052=false, next911=false, next875=false;
  bool  next680=false,  next530=false, next263=false, next233=false;
  float pbeta,p;
  pbeta=100.*GetRandom();
  if (pbeta<=0.42){
     beta(0.181,0.,0.);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=0.72){
        nucltransK(1.590,0.035,1.0e-3,1.1e-4);
        return;
     }
     else if (p<=33.78){
        nucltransK(1.060,0.035,1.7e-3,0.);
        next530=true;
     }
     else if (p<=85.06){
        nucltransK(0.910,0.035,2.0e-2,0.);
        next680=true;
     }
     else if (p<=90.33){
        nucltransK(0.679,0.035,5.5e-3,0.);
        next911=true;
     }
     else if (p<=98.95){
        nucltransK(0.538,0.035,8.5e-3,0.);
        next1052=true;
     }
     else{
        nucltransK(0.204,0.035,1.5e-1,0.);
        next1386=true;
     }
  }
  else if (pbeta<= 1.68){
     beta(0.385,0.,0.);
     next1386=true;
  }
  else if (pbeta<= 2.08){
     beta(0.421,0.,0.);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=37.48){
        nucltransK(1.350,0.035,9.9e-4,0.6e-4);
        return;
     }
     else if (p<=40.53){
        nucltransK(1.088,0.035,1.3e-3,0.);
        next263=true;
     }
     else if (p<=79.26){
        nucltransK(0.821,0.035,3.1e-3,0.);
        next530=true;
     }
     else if (p<=90.00){
        nucltransK(0.670,0.035,5.0e-3,0.);
        next680=true;
     }
     else{
        nucltransK(0.439,0.035,1.4e-2,0.);
        next911=true;
     }
  }
  else if (pbeta<= 5.89){
     beta(0.473,0.,0.);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=62.29){
        nucltransK(1.298,0.035,9.6e-4,0.5e-4);
        return;
     }
     else if (p<=62.52){
        nucltransK(1.036,0.035,1.5e-3,0.);
        next263=true;
     }
     else if (p<=74.71){
        nucltransK(0.768,0.035,3.5e-3,0.);
        next530=true;
     }
     else if (p<=89.13){
        nucltransK(0.618,0.035,5.9e-3,0.);
        next680=true;
     }
     else if (p<=97.51){
        nucltransK(0.423,0.035,1.5e-2,0.);
        next875=true;
     }
     else if (p<=99.07){
        nucltransK(0.387,0.035,2.0e-2,0.);
        next911=true;
     }
     else{
        nucltransK(0.246,0.035,8.0e-2,0.);
        next1052=true;
     }
  }
  else if (pbeta<= 9.10){
     beta(0.535,0.,0.);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=47.94){
        nucltransK(1.236,0.035,1.0e-3,0.1e-4);
        return;
     }
     else if (p<=95.88){
        nucltransK(0.707,0.035,4.3e-3,0.);
        next530=true;
     }
     else if (p<=96.51){
        nucltransK(0.556,0.035,7.0e-3,0.);
        next680=true;
     }
     else{
        nucltransK(0.361,0.035,3.0e-2,0.);
        next875=true;
     }
  }
  else if (pbeta<= 9.73){
     beta(0.719,0.,0.);
     next1052=true;
  }
  else if (pbeta<=13.94){
     beta(0.896,0.,0.);
     next875=true;
  } 
  else if (pbeta<=15.77){
     beta(1.027,0.,0.); 
     fThlev=0.;
     nucltransK(0.510,0.035,9.0e-3,0.);
     next233=true;  
  }
  else if (pbeta<=98.96){
     beta(1.241,0.,0.);
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=99.87){
        nucltransK(0.530,0.035,8.2e-3,0.);
        return;
     }
     else{
        nucltransK(0.267,0.035,6.0e-2,0.);
        next263=true;
     }
  }
  else{
      beta(1.538,0.,0.);
      next233=true;
  }
  if (next1386){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=0.72){
        nucltransK(1.386,0.035,3.5e-3,0.);
        return;
     }
     else{
        nucltransK(0.856,0.035,2.6e-3,0.);
        next530=true;
     }
  }
  if (next1052){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=90.26){
        nucltransK(1.052,0.035,2.0e-3,0.);
        return;
     }
     else if (p<=98.37){
        nucltransK(0.790,0.035,3.0e-3,0.);
        next263=true;
     }
     else{
        nucltransK(0.372,0.035,2.5e-2,0.);
        next680=true;
     }
  }
  if (next911){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=31.08){
        nucltransK(0.911,0.035,3.0e-3,0.);
        return;
     }
     else if (p<=69.59){
        nucltransK(0.649,0.035,6.0e-3,0.);
        next263=true;
     }
     else{
        nucltransK(0.382,0.035,2.0e-2,0.);
        next530=true; 
     }
  }
  if (next875){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=97.75){
        nucltransK(0.875,0.035,2.2e-3,0.);
        return;
     }
     else{
        nucltransK(0.345,0.035,3.0e-2,0.);
        next530=true;
     }
  }
  if (next680){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=77.66){
        nucltransK(0.680,0.035,5.4e-3,0.);
        return;
     }
     else if (p<=96.42){
        nucltransK(0.418,0.035,1.6e-2,0.);
        next263=true;
     }
     else{
        nucltransK(0.150,0.035,3.0e-1,0.);
        next530=true;
     }
  }
  if (next530){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=99.87){
        nucltransK(0.530,0.035,8.2e-3,0.);
        return;
     }
     else{
        nucltransK(0.267,0.035,6.0e-2,0.);
        next263=true;
     }
  }
  if (next263){
     fThlev=0.;
     nucltransK(0.263,0.035,5.8e-2,0.);
     return; 
  }
  if (next233){
      fThlev=189216.;
      nucltransK(0.233,0.035,8.8,0.);
      return;
     
  }
}
//-----------------------------------------------

void Decay::I134(float tcnuc)
{
  // Scheme of I134 decay in accordance with Yu.V.Sergeenkov, 
  // NDS 71(1994)557.
  // VIT, 8.10.2002.
  fThnuc=3168.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fThlev=0.;
  fTclev=0.;
  bool next3314=false, next2654=false, next2588=false, next2548=false, next2408=false;
  bool next2353=false, next2302=false, next2272=false, next2137=false, next1920=false;
  bool next1731=false, next1614=false, next847=false;
  float pbeta,p; 
  pbeta=100.*GetRandom();
  if (pbeta<=0.368){
     beta(0.678,0.,0.); 
     p=100.*GetRandom();
     if (p<=5.17){
        nucltransK(2.646,0.035,1.5e-4,5.2e-4);
        next847=true; 
     }
     else{
        nucltransK(1.190,0.035,2.0e-3,0.2e-4);
        next2302=true;
     }
  }
  else if (pbeta<= 0.435){
     beta(0.693,0.,0.);
     p=100.*GetRandom();
     if (p<=18.77){
        nucltransK(2.630,0.035,1.5e-4,5.1e-4);
        next847=true; 
     }
     else{
        nucltransK(0.162,0.035,2.5e-1,0.);
        next3314=true;
     }
  }
  else if (pbeta<= 2.216){
     beta(0.795,0.,0.);
     p=100.*GetRandom();
     if (p<=21.79){
        nucltransK(1.644,0.035,1.5e-3,1.2e-4);
        next1731=true;
     }
     else if (p<=33.52){
        nucltransK(1.239,0.035,1.7e-3,0.4e-4);
        next2137=true;
     }
     else if (p<=78.21){
        nucltransK(1.103,0.035,2.1e-3,0.1e-4);
        next2272=true;
     }
     else{
        nucltransK(0.967,0.035,2.5e-3,0.);
        next2408=true;
     }
  }
  else if (pbeta<= 3.300){
     beta(0.810,0.,0.);
     p=100.*GetRandom();
     if (p<=6.16){
        nucltransK(2.513,0.035,2.0e-4,4.6e-4);
        next847=true;
     }
     else if (p<=23.64){
        nucltransK(1.629,0.035,9.0e-4,1.1e-4);
        next1731=true;
     }
     else{
        nucltransK(0.707,0.035,5.0e-3,0.);
        next2654=true;
     }
  }
  else if (pbeta<= 3.509){
     beta(0.856,0.,0.);
     next3314=true;
  }
  else if (pbeta<= 3.567){
     beta(0.870,0.,0.);
     nucltransK(2.453,0.035,2.0e-4,4.4e-4);
     next847=true;
  }
  else if (pbeta<= 3.786){
     beta(0.914,0.,0.);
     p=100.*GetRandom();
     if (p<=35.78){
        nucltransK(2.409,0.035,2.1e-4,4.1e-4);
        next847=true;
     }
     else{
        nucltransK(1.336,0.035,1.2e-3,0.6e-4);
        next1920=true;
     }
  }
  else if (pbeta<= 5.129){
     beta(1.086,0.,0.);
     p=100.*GetRandom();
     if (p<=3.95){
        nucltransK(2.237,0.035,3.0e-4,3.2e-4);
        next847=true;
     }
     else if (p<=59.80){
        nucltransK(1.470,0.035,1.0e-3,0.9e-4);
        next1614=true;
     }
     else if (p<=90.32){
        nucltransK(1.353,0.035,1.2e-3,0.6e-4);
        next1731=true;
     }
     else{
        nucltransK(1.164,0.035,1.6e-3,0.2e-4);
        next1920=true;
     }
  }
  else if (pbeta<=35.372){
     beta(1.303,0.,0.);
     p=100.*GetRandom();
     if (p<=0.62){
        nucltransK(2.021,0.035,3.0e-4,2.3e-4);
        next847=true;
     } 
     else if (p<=30.54){
        nucltransK(1.136,0.035,1.6e-3,0.3e-4);
        next1731=true;
     }
     else if (p<=43.73){
        nucltransK(0.948,0.035,2.4e-3,0.);
        next1920=true;
     }
     else if (p<=49.75){
        nucltransK(0.731,0.035,3.4e-3,0.);
        next2137=true; 
     }
     else if (p<=86.25){
        nucltransK(0.595,0.035,7.2e-3,0.);
        next2272=true;
     }
     else if (p<=93.65){
        nucltransK(0.514,0.035,8.0e-3,0.);
        next2353=true;
     }
     else if (p<=97.99){
        nucltransK(0.459,0.035,1.1e-2,0.);
        next2408=true;
     }
     else if (p<=99.50){
        nucltransK(0.320,0.035,3.7e-2,0.);
        next2548=true;
     }
     else{
        nucltransK(0.279,0.035,5.2e-2,0.);
        next2588=true; 
     }
  }
  else if (pbeta<=35.889){
     beta(1.397,0.,0.); 
     p=100.*GetRandom();
     if (p<=34.62){
        nucltransK(1.926,0.035,4.0e-4,2.0e-4);
        next847=true;
     }
     else{
        nucltransK(1.159,0.035,7.0e-4,0.3e-4);
        next1614=true;
     }
  }
  else if (pbeta<=42.754){
     beta(1.516,0.,0.);
     next2654=true;
  }
  else if (pbeta<=58.871){
     beta(1.582,0.,0.);
     next2588=true;
  }
  else if (pbeta<=62.950){
     beta(1.622,0.,0.);
     next2548=true;
  }
  else if (pbeta<=69.516){
     beta(1.762,0.,0.);
     next2408=true;
  }
  else if (pbeta<=80.460){
     beta(1.817,0.,0.);
     next2353=true;
  }
  else if (pbeta<=82.689){
     beta(1.868,0.,0.);
     next2302=true;
  }
  else if (pbeta<=84.281){
     beta(1.898,0.,0.);
     next2272=true;
  }
  else if (pbeta<=87.564){
     beta(2.250,0.,0.);
     next1920=true;
  } 
  else{
     beta(2.439,0.,0.);
     next1731=true;
  }

  if (next3314){
     p=100.*GetRandom();
     if (p<=63.81){
        nucltransK(2.467,0.035,2.0e-4,4.5e-4);
        next847=true;
     }
     else{
        nucltransK(1.395,0.035,1.2e-3,0.7e-4);
        next1920=true;
     }
  }
  if (next2654){
     p=100.*GetRandom();
     if (p<=71.82){
        nucltransK(1.807,0.035,6.0e-4,1.7e-4);
        next847=true;
     }
     else if (p<=98.18){
        nucltransK(1.040,0.035,2.0e-3,0.);
        next1614=true;
     }
     else{
        nucltransK(0.922,0.035,2.1e-3,0.);
        next1731=true;
     }
  }
  if (next2588){
     p=100.*GetRandom();
     if (p<=15.66){
        nucltransK(1.741,0.035,6.0e-4,0.8e-4);
        next847=true;
     }
     else if (p<=44.90){
        nucltransK(0.975,0.035,1.7e-3,0.);
        next1614=true;
     }
     else if (p<=85.87){
        nucltransK(0.857,0.035,2.9e-3,0.);
        next1731=true;
     }
     else{
        nucltransK(0.235,0.035,8.5e-2,0.);
        next2353=true;
     }
  }
  if (next2548){
     p=100.*GetRandom();
     if (p<=13.78){
        nucltransK(0.816,0.035,4.0e-3,0.);
        next1731=true;
     }
     else if (p<=63.11){
        nucltransK(0.628,0.035,4.9e-3,0.);
        next1920=true;
     }
     else if (p<=75.78){
        nucltransK(0.411,0.035,2.0e-2,0.);
        next2137=true;
     }
     else{
        nucltransK(0.139,0.035,4.5e-1,0.);
        next2408=true;
     }
  }
  if (next2408){
     p=100.*GetRandom();
     if (p<=84.57){
        nucltransK(0.677,0.035,5.3e-3,0.);
        next1731=true; 
     }
     else{
        nucltransK(0.489,0.035,9.6e-3,0.);
        next1920=true;
     }
  }
  if (next2353){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=4.42){
        nucltransK(0.739,0.035,3.3e-3,0.);
        next1614=true;
     }
     else if (p<=73.01){
        nucltransK(0.622,0.035,6.1e-3,0.);
        next1731=true;
     }
     else{
        nucltransK(0.433,0.035,1.6e-2,0.);
        next1920=true;
     }
  }
  if (next2302){
     p=100.*GetRandom();
     if (p<=88.08){
        nucltransK(1.455,0.035,1.1e-3,0.8e-4);
        next847=true;
     }
     else{
        nucltransK(0.571,0.035,9.0e-3,0.);
        next1731=true;
     }
  }
  if (next2272){
     p=100.*GetRandom();
     if (p<=57.01){
        nucltransK(0.541,0.035,7.7e-3,0.);
        next1731=true;
     }
     else{
        nucltransK(0.135,0.035,3.5e-1,0.);
        next2137=true;
     }
  }
  if (next2137){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=96.64){
        nucltransK(0.405,0.035,1.8e-2,0.);
        next1731=true;
     }
     else{
        nucltransK(0.217,0.035,1.2e-1,0.);
        next1920=true;
     }
  }
  if (next1920){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=94.36){
        nucltransK(1.073,0.035,1.9e-3,0.1e-4);
        next847=true;
     }
     else{
        nucltransK(0.188,0.035,1.7e-1,0.);
        next1731=true; 
     }
  }
  if (next1731){
     fThlev=0.;
     nucltransK(0.884,0.035,2.2e-3,0.);
     next847=true;
  }
  if (next1614){
     fThlev=0.;
     p=100.*GetRandom();
     if (p<=50.83){
        nucltransK(1.614,0.035,3.0e-4,0.5e-4);
        return; 
     } 
     else{
        nucltransK(0.767,0.035,3.2e-3,0.);
        next847=true;
     }
  }
  if (next847){
     fThlev=1.9e-12;
     nucltransK(0.847,0.035,2.4e-3,0.);
     return;
  }
} 
//-----------------------------------------------

void Decay::I135(float tcnuc)
{
  // Scheme of I135 decay (Yu.V.Sergeenkov et al., Nucl. Data Sheets 84(1998)115).
  fThnuc=23652.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fThlev=0.;
  fTclev=0.;
  bool next2093=false, next2049=false, next1927=false, next1791=false;
  bool next1781=false,  next1678=false, next1565=false, next1544=false;
  bool next1458=false, next1448=false, next1260=false, next1132=false;
  bool next527=false,  next288=false;
  float pbeta,p; 
  pbeta=100.*GetRandom();
  if (pbeta<=0.013){
     beta(0.170,0.,0.);
     p=100.*GetRandom();
     if (p<=9.72){
        nucltransK(2.477,0.035,2.0e-4,4.5e-4);
        return; 
     } 
     else{
        nucltransK(2.189,0.035,3.0e-4,3.1e-4);
        next288=true;
     }
  }
  else if (pbeta<= 0.155){
      beta(0.173,0.,0.);
      p=100.*GetRandom();
      if (p<=45.00){
         nucltransK(1.948,0.035,2.4e-4,2.1e-4);
         next527=true;
      }
      else{
         nucltransK(1.344,0.035,4.0e-4,0.7e-4);
         next1132=true; 
      }
  }
  else if (pbeta<= 0.282){
      beta(0.182,0.,0.);
      p=100.*GetRandom();
      if (p<=56.69){
         nucltransK(2.466,0.035,4.0e-4,5.3e-4);
         return; 
      }
      else if (p<=81.89){
         nucltransK(1.335,0.035,1.2e-3,0.6e-4);
         next1132=true;
      }
      else{
         nucltransK(0.685,0.035,4.0e-3,0.);
         next1781=true;
      }
  }
  else if (pbeta<= 0.422){
      beta(0.201,0.,0.);
      p=100.*GetRandom();
      if (p<=46.81){
         nucltransK(1.316,0.035,4.5e-4,0.5e-4);
         next1132=true;
      }
      else{
         nucltransK(0.656,0.035,1.6e-3,0.);
         next1791=true;
      }
  }
  else if (pbeta<= 1.449){
      beta(0.239,0.,0.);
      p=100.*GetRandom();
      if (p<=93.48){
         nucltransK(2.409,0.035,3.5e-4,4.1e-4);
         return;
      }
      else if (p<=96.40){
         nucltransK(0.960,0.035,2.6e-3,0.);
         next1448=true;
      }
      else{
         nucltransK(0.617,0.035,7.0e-3,0.);
         next1791=true;
      }
  }
  else if (pbeta<= 2.366){
      beta(0.276,0.,0.);
      p=100.*GetRandom();
      if (p<=0.66){
         nucltransK(1.845,0.035,2.6e-4,1.8e-4);
         next527=true;
      }
      else{
         nucltransK(1.240,0.035,1.4e-3,0.5e-4);
         next1132=true;
      }
  }
  else if (pbeta<= 3.752){
      beta(0.291,0.,0.);
      p=100.*GetRandom();
      if (p<=41.91){
         nucltransK(1.831,0.035,2.6e-4,1.6e-4);
         next527=true;
      }
      else if (p<=45.02){
         nucltransK(1.226,0.035,1.4e-3,0.4e-4);
         next1132=true;
      }
      else if (p<=51.45){
         nucltransK(1.097,0.035,1.4e-3,0.);
         next1260=true;
      }
      else if (p<=55.42){
         nucltransK(0.679,0.035,6.0e-3,0.);
         next1678=true;
      }
      else if (p<=64.74){
         nucltransK(0.576,0.035,8.5e-3,0.);
         next1781=true;
      }
      else if (p<=86.41){
         nucltransK(0.430,0.035,1.8e-2,0.);
         next1927=true;
      }
      else{
         nucltransK(0.264,0.035,6.5e-2,0.);
         next2093=true;
      }
  }
  else if (pbeta<= 8.518){
      beta(0.393,0.,0.);
      p=100.*GetRandom();
      if (p<=12.88){
         nucltransK(2.255,0.035,3.7e-4,3.8e-4);
         return;
      }
      else if (p<=89.34){
         nucltransK(1.124,0.035,1.7e-3,0.2e-4);
         next1132=true;
      }
      else if (p<=92.51){
         nucltransK(0.995,0.035,1.9e-3,0.);
         next1260=true;
      } 
      else if (p<=93.48){
         nucltransK(0.807,0.035,2.6e-3,0.);
         next1448=true; 
      }
      else if (p<=97.07){
         nucltransK(0.798,0.035,4.0e-3,0.);
         next1458=true;
      }
      else if (p<=99.79){
         nucltransK(0.690,0.035,5.5e-3,0.);
         next1565=true;
      }
      else{
         nucltransK(0.163,0.035,2.2e-1,0.);
         next2093=true; 
      }
  }
  else if (pbeta<=15.896){
      beta(0.415,0.,0.);
      p=100.*GetRandom();
      if (p<=55.75){
         nucltransK(1.706,0.035,2.8e-4,1.0e-4);
         next527=true;
      }
      else if (p<=77.64){
         nucltransK(1.102,0.035,1.4e-3,0.1e-4);
         next1132=true;
      }
      else if (p<=94.09){
         nucltransK(0.973,0.035,1.8e-3,0.);
         next1260=true;
      }
      else if (p<=98.39){
         nucltransK(0.452,0.035,1.6e-2,0.);
         next1781=true;
      }
      else if (p<=99.68){
         nucltransK(0.306,0.035,4.0e-2,0.);
         next1927=true;
      }
      else{
         nucltransK(0.184,0.035,4.0e-2,0.);
         next2049=true;
      }
  }
  else if (pbeta<=15.919){
      beta(0.496,0.,0.);
      nucltransK(2.152,0.035,2.0e-4,2.9e-4);
      return;
  }
  else if (pbeta<=15.989){
      beta(0.536,0.,0.); 
      nucltransK(2.112,0.035,2.0e-4,2.7e-4);
      return;
  }
  else if (pbeta<=17.574){
      beta(0.555,0.,0.);
      next2093=true;
  }
  else if (pbeta<=18.681){
      beta(0.602,0.,0.);
      p=100.*GetRandom();
      if (p<=79.31){
         nucltransK(2.046,0.035,4.5e-4,2.3e-4);
         return;
      }
      else if (p<=93.16){
         nucltransK(0.785,0.035,3.6e-3,0.);
         next1260=true;
      }
      else if (p<=97.90){
         nucltransK(0.588,0.035,8.5e-3,0.);
         next1458=true;
      }
      else{
         nucltransK(0.255,0.035,6.5e-2,0.);
         next1791=true;
      }
  }
  else if (pbeta<=26.657){
      beta(0.680,0.,0.);
      p=100.*GetRandom();
      if (p<=0.21){
         nucltransK(1.442,0.035,0.8e-4,1.0e-3);
         next527=true;
      }
      else if (p<=84.72){
         nucltransK(0.837,0.035,2.5e-3,0.);
         next1132=true;
      }
      else if (p<=93.05){
         nucltransK(0.708,0.035,5.0e-3,0.);
         next1260=true;
      }
      else if (p<=95.98){
         nucltransK(0.403,0.035,1.8e-2,0.);
         next1565=true;
      }
      else{
         nucltransK(0.290,0.035,4.6e-2,0.);
         next1678=true;
      }
  }
  else if (pbeta<=26.707){
      beta(0.721,0.,0.);
      next1927=true;
  }
  else if (pbeta<=27.315){
      beta(0.754,0.,0.);
      nucltransK(1.368,0.035,9.0e-4,0.2e-4);
      next527=true;  
  }
  else if (pbeta<=36.089){
      beta(0.857,0.,0.);
      next1791=true;
  }
  else if (pbeta<=57.825){
     beta(0.970,0.,0.);
     next1678=true;
  }
  else if (pbeta<=65.801){
     beta(1.083,0.,0.);
     next1565=true;
  }
  else if (pbeta<=73.279){
     beta(1.190,0.,0.);
     next1458=true;
  }
  else if (pbeta<=96.810){
     beta(1.388,0.,0.);
     next1260=true;
  }
  else if (pbeta<=98.106){
     beta(1.516,0.,0.);
     next1132=true;
  }
  else{
     beta(2.121,0.,0.);
     next527=true;
  }
 
  if (next2093){
      p=100.*GetRandom();
      if (p<=72.79){
         nucltransK(1.566,0.035,2.8e-4,0.7e-4);
         next527=true;
      }
      else if (p<=81.26){
         nucltransK(0.961,0.035,2.6e-3,0.);
         next1132=true;
      }
      else if (p<=98.25){
         nucltransK(0.415,0.035,2.0e-2,0.);
         next1678=true; 
      }
      else{
         nucltransK(0.166,0.035,2.2e-1,0.);
         next1927=true;
      } 
  }
  if (next2049){
     nucltransK(1.522,0.035,9.0e-4,0.9e-4);
     next527=true;
  }
  if (next1927){
      p=100.*GetRandom();
      if (p<= 58.50){ 
         nucltransK(1.927,0.035,5.0e-4,1.4e-4);
         return;
      }
      else if (p<=63.05){
         nucltransK(0.796,0.035,4.0e-3,0.);
         next1132=true;
      }
      else{
         nucltransK(0.362,0.035,2.8e-2,0.);
         next1565=true;
      }
  }
  if (next1791){
      p=100.*GetRandom();
      if (p<=86.68){
        nucltransK(1.791,0.035,6.5e-4,1.6e-4);
        return;
     }
     else if (p<=98.81){
        nucltransK(1.503,0.035,7.5e-4,0.3e-4);
        next288=true;
     }
     else if (p<=99.17){
        nucltransK(0.531,0.035,9.3e-3,0.);
        next1260=true;
     }
     else if (p<=99.18){
        nucltransK(0.343,0.035,3.0e-2,0.);
        next1448=true;
     }
     else if (p<=99.60){
        nucltransK(0.334,0.035,3.5e-2,0.);
        next1458=true;
     }
     else if (p<=99.92){
        nucltransK(0.248,0.035,8.0e-3,0.);
        next1544=true; 
     }
     else{
        nucltransK(0.113,0.035,6.0e-1,0.);
        next1678=true;
     }
  }
  if (next1781){
     p=100.*GetRandom();
      if (p<=1.29){
         nucltransK(1.255,0.035,4.7e-4,0.3e-4);
         next527=true;
      }
      else{ 
         nucltransK(0.645,0.035,4.5e-3,0.);
         next1132=true;
      }
  }
  if (next1678){ 
     p=100.*GetRandom();
     if (p<=42.57){
        nucltransK(1.678,0.035,5.6e-4,0.7e-4);
        return;
     }
     else if (p<=74.58){
        nucltransK(0.547,0.035,9.1e-3,0.);
        next1132=true;
     }
     else if (p<=90.20){
        nucltransK(0.418,0.035,1.6e-2,0.);
        next1260=true;
     }
     else if (p<=91.38){
        nucltransK(0.230,0.035,1.0e-1,0.);
        next1448=true;
     }
     else if (p<=99.94){
        nucltransK(0.221,0.035,1.0e-1,0.);
        next1458=true;
     }
     else{
        nucltransK(0.113,0.035,6.0e-1,0.);
        next1565=true;
     }
  }
  if (next1565){
     p=100.*GetRandom();
     if (p<=93.00){
        nucltransK(1.039,0.035,6.4e-4,0.);
        next527=true; 
     }
     else if (p<=99.62){
        nucltransK(0.434,0.035,1.6e-2,0.);
        next1132=true;
     }
     else{
        nucltransK(0.305,0.035,4.0e-2,0.);
        next1260=true;
     }
  }
  if (next1544){
     p=100.*GetRandom();
     if (p<=81.25){
        nucltransK(1.544,0.035,9.0e-4,1.0e-4);
        return; 
     }
     else{
        nucltransK(1.255,0.035,1.2e-3,0.5e-4);
        next288=true;
     }
  }
  if (next1458){
     p=100.*GetRandom();
     if (p<=90.44){
        nucltransK(1.458,0.035,7.0e-4,0.9e-4);
        return;
     }
     else if (p<=99.59){
        nucltransK(1.169,0.035,1.3e-3,0.);
        next288=true;
     } 
     else if (p<=99.61){
        nucltransK(0.326,0.035,3.5e-2,0.);
        next1132=true;
     }
     else{
        nucltransK(0.197,0.035,1.5e-1,0.);
        next1260=true;
     }
  } 
  if (next1448){
     p=100.*GetRandom();
     if (p<=75.65){
        nucltransK(1.448,0.035,7.0e-4,0.9e-4);
        return;
     }
     else{
        nucltransK(1.160,0.035,1.3e-3,0.);
        next288=true;
     }
  }
  if (next1260){
     p=100.*GetRandom();
     if (p<=96.99){
        nucltransK(1.260,0.035,1.2e-3,0.1e-4);
        return;
     }
     else{
        nucltransK(0.972,0.035,1.8e-3,0.);
        next288=true;
     }
  }
  if (next1132){
     nucltransK(1.132,0.035,1.3e-3,0.);
     return;
  }
  if (next527){
     fThlev=917.4;
     nucltransK(0.527,0.035,2.4e-1,0.);
     return; 
  }
  if (next288){
     nucltransK(0.288,0.035,4.7e-2,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::K40(float tcnuc)
{
  // Scheme of K40 decay (J.A.Cameron et al., ENSDF 29.10.2002).
  // 3rd forbidden unique shape for beta decay, VIT, 27.10.2006
  // in accordance with H.Daniel, RMP 40(1968)659 and W.H.Kelly et al., 
  // Nucl. Phys. 11(1959)492
  fThnuc=3.94775e+16;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fZdtr=20;
  float pdecay=100.*GetRandom();
  if (pdecay<=89.140){ 
     beta2f(1.311,0.,0.,3,7.,7.,1.,0.);
  }
  else if (pdecay<=99.800){  // 10.660% ec 40Ar(1461)
     particle(1,0.003,0.003,0.,pi,0.,twopi,0,0);
     fThlev=1.12e-12;
     nucltransK(1.461,0.003,3.0e-5,7.3e-5);
  }
  else if (pdecay<=99.999){  // 0.199% ec 40Ar(gs)
     particle(1,0.003,0.003,0.,pi,0.,twopi,0,0);
  }
  else{                      // 0.001% b+ 40Ar(gs) 
     fZdtr=-18;
     beta(0.483,0.,0.);
  }
  return;
}
//-----------------------------------------------

void Decay::K42(float tcnuc)
{
  // Scheme of K42 decay (B.Singh et al., NDS 92(2001)1
  // with additional information from ToI-1978 and ToI-1998).
  // VIT, 31.03.2006; 29.10.2006.
  fThnuc=44496.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  bool next2753=false, next2424=false, next1525=false;
  float pbeta,p;
  pbeta=100.*GetRandom();
  if (pbeta<=0.07){
     beta(0.080,0.,0.);
     fThlev=250.e-15;  // ToI-1998
     p=100.*GetRandom();
     if (p<=63.66){
        nucltransK(1.921,0.004,1.2e-5,5.8e-4);
        next1525=true;
     }
     else if (p<=94.88){
        nucltransK(1.021,0.004,3.9e-5,0.);
        next2424=true;
     }
     else{
        nucltransK(0.692,0.004,8.7e-5,0.);
        next2753=true;
     }
  }
  else if (pbeta<= 0.12){
     beta(1.101,0.,0.); 
     next2424=true;
  }
  else if (pbeta<= 0.46){
     beta1fu(1.688,0.,0.,0.,0.,0.,0.); 
     fThlev=0.33e-9;
     p=100.*GetRandom();
     if (p<=2.1){  // ToI-1978
        p=100.*GetRandom();
        if (p<=90.){
           pair(0.815);
        }
        else{
           particle(3,1.833,1.833,0.,pi,0.,twopi,fTclev,fThlev);
           particle(1,0.004,0.004,0.,pi,0.,twopi,0,0);
        }
        return;
     }
     else{
        nucltransK(0.313,0.004,3.4e-3,0.);
        next1525=true;
     }
  }
  else if (pbeta<=18.10){
     beta1f(2.000,0.,0.,0.81,0.15,-0.02,0.);
     next1525=true;
  }
  else{
     beta1fu(3.525,0.,0.,0.,-0.01,0.,0.);
     return;
  }
  if (next2753){
     fThlev=3.0e-12; // ToI-1998
     nucltransK(1.228,0.004,5.6e-5,1.4e-5);
     next1525=true;
  }
  if (next2424){
     fThlev=140.e-15; // ToI-1998
     p=100.*GetRandom();
     if (p<=27.78){
        nucltransK(2.424,0.004,1.5e-5,5.2e-4);
        return;
     }
     else{
        nucltransK(0.900,0.004,8.3e-5,0.);
        next1525=true;
     }
  }
  if (next1525){
     fThlev=0.82e-12;  // ToI-1998
     nucltransK(1.525,0.004,3.6e-5,9.8e-5);
     return;
  }
}
//-----------------------------------------------

void Decay::Pa234(float tcnuc)
{
  // Model (not the exact) scheme of Pa234m decay ("Table of Isotopes",
  // 7th ed., 1978): decays of Pa234m to excited levels of U234 with energies
  // greater than 1.045 MeV are not considered (p=0.20%).
  // VIT, 18.08.1992, 22.10.1995;
  // VIT, 3.11.2006 (update to NDS 71(1994)181 -however, levels 
  // with E_exc>1.045 MeV are still not considered; 
  // change of beta shape of the gs-to-gs beta decay to the 1st forbidden 
  // with exp. corr.).
  fThnuc=70.2;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next927=false, next852=false, next849=false, next810=false;
  bool next786=false, next143=false, next43=false;
  float pdecay, pbeta, p;
  pdecay=100.*GetRandom();
  if (pdecay<=0.16){  // IT to Pa234
     nucltransK(0.074,0.021,1.1e+1,0.);
     return;
  }
  pbeta=100.*GetRandom();  // beta decay to U234
  if (pbeta<=1.008){
     beta(1.224,0.,0.); 
     p=100.*GetRandom();
     if (p<=83.7){
        nucltransK(1.001,0.116,1.1e-2,0.);
        next43=true;
     }
     else if (p<=91.3){
        nucltransK(0.258,0.116,5.5e-2,0.);
        next786=true;
     }
     else if (p<=99.9){
        particle(3,0.120,0.120,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.116,0.116,0.,pi,0.,twopi,0.,0.);
        next810=true;
     }
     else{
        nucltransK(0.193,0.116,8.8e-1,0.);
        next852=true;
     }
  }
  else if (pbeta<=1.012){
     beta(1.280,0.,0.); 
     p=100.*GetRandom();
     if (p<=44.8){
        nucltransK(0.946,0.116,4.1e-3,0.);
        next43=true;
     }
     else if (p<=56.1){
        nucltransK(0.203,0.116,1.5e0,0.);
        next786=true;
     }
     else if (p<=92.3){
        nucltransK(0.140,0.116,5.5e0,0.);
        next849=true;
     }
     else{
        nucltransK(0.063,0.022,4.3e-1,0.);
        next927=true;
     }
  }
  else if (pbeta<=1.022){
     beta(1.417,0.,0.);
     next852=true;
  }
  else if (pbeta<=1.712){
     beta(1.459,0.,0.);
     next810=true;
  }
  else if (pbeta<=1.744){
     beta(1.483,0.,0.);
     next786=true;
  }
  else{
     beta1f(2.269,0.,0.,0.,-0.09,0.,0.);
     return;
  }
 
  if (next927){
     p=100.*GetRandom();
     if (p<=40.4){
        nucltransK(0.927,0.116,1.3e-2,0.);
        return;
     }
     else if (p<=98.7){
        nucltransK(0.883,0.116,1.4e-2,0.);
        next43=true;
     }
     else{
        nucltransK(0.783,0.116,1.8e-2,0.);
        next143=true;
     }
  } 
  if (next852){
     p=100.*GetRandom();
     if (p<=18.6){ 
        nucltransK(0.852,0.116,1.5e-2,0.);
        return;
     }
     else if (p<=64.6){
        nucltransK(0.808,0.116,4.2e0,0.);
        next43=true;
     }
     else{
        nucltransK(0.042,0.022,1.0e+1,0.);
        next810=true;
     }
  }
  if (next849){
     p=100.*GetRandom();
     if (p<=51.8){
        nucltransK(0.806,0.116,5.5e-3,0.);
        next43=true;
     }
     else{
        nucltransK(0.706,0.116,7.0e-3,0.);
        next143=true;
     }
  }
  if (next810){
     p=100.*GetRandom();
     if (p<=63.0){
        particle(3,0.694,0.694,0.,pi,0.,twopi,fTclev,fThlev);
        particle(1,0.116,0.116,0.,pi,0.,twopi,0.,0.);
        return;
     }
     else{
        nucltransK(0.766,0.116,1.9e-2,0.);
        next43=true;
     }
  }
  if (next786){
     p=100.*GetRandom();
     if (p<=37.6){
        nucltransK(0.786,0.116,5.8e-3,0.);
        return;
     }
     else{
        nucltransK(0.743,0.116,6.4e-3,0.);
        next43=true;
     }
  }
  if (next143){
     nucltransK(0.100,0.022,1.4e+1,0.);
     next43=true;
  }
  if (next43){
     fThlev=0.25e-9;
     nucltransK(0.043,0.022,7.2e+2,0.);
     return; 
  }
}
//-----------------------------------------------

void Decay::Pb211(float tcnuc)
{
  // Scheme of Pb211 decay ("Table of Isotopes", 7th ed., 1978).
  // VIT, 14.08.1992, 22.10.1995;
  // VIT, 31.10.2006 (updated to NDS 103(2004)183)
  fThnuc=2166.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next1014=false, next767=false, next405=false;
  float pbeta, p;   
  pbeta=100.*GetRandom();  // beta decay to U234
  if (pbeta<=0.0186){
     beta(0.096,0.,0.);
     p=100.*GetRandom();
     if (p<=36.76){
        nucltransK(1.271,0.091,1.0e-2,0.2e-4);
        return;
     }
     else if (p<=68.65){
        nucltransK(0.866,0.091,2.8e-2,0.);
        next405=true; 
     }
     else{
        nucltransK(0.504,0.091,1.1e-1,0.);
        next767=true;
     }
  }
  else if (pbeta<=0.0199){
     beta(0.133,0.,0.); 
     nucltransK(1.234,0.091,1.1e-2,0.1e-4);
     return;
  }
  else if (pbeta<=0.0369){
     beta(0.171,0.,0.);
     p=100.*GetRandom();
     if (p<=62.96){
        nucltransK(1.196,0.091,1.2e-2,0.1e-4);
        return;
     }
     else{
        nucltransK(0.430,0.091,1.8e-1,0.);
        next767=true;
     }
  }
  else if (pbeta<=0.8660){
     beta(0.258,0.,0.);
     p=100.*GetRandom();
     if (p<=13.9){
        nucltransK(1.109,0.091,1.5e-2,0.4e-6);
        return;
     }
     else if (p<=72.8){
        nucltransK(0.705,0.091,5.0e-2,0.);
        next405=true;
     }
     else if (p<=77.0){
        nucltransK(0.343,0.091,3.2e-1,0.);
        next767=true;
     }
     else{
        nucltransK(0.095,0.016,9.6e+0,0.);
        next1014=true;
     }
  }
  else if (pbeta<=0.8706){
     beta(0.263,0.,0.);
     nucltransK(1.104,0.091,1.5e-2,0.3e-6);
     return;
  }
  else if (pbeta<=0.9265){
     beta(0.286,0.,0.);
     p=100.*GetRandom();
     if (p<=21.8){
        nucltransK(1.081,0.091,1.6e-2,0.1e-4);
        return;
     }
     else if (p<=44.9){
        nucltransK(0.677,0.091,5.3e-2,0.);
        next405=true;
     }
     else{
        nucltransK(0.314,0.091,4.1e-1,0.);
        next767=true;
     }
  }
  else if (pbeta<=0.9485){
     beta(0.416,0.,0.);
     nucltransK(0.951,0.091,2.2e-2,0.);
     return;
  }
  else if (pbeta<=7.2616){
     beta(0.535,0.,0.);
     p=100.*GetRandom();
     if (p<=57.4){
        nucltransK(0.832,0.091,2.9e-2,0.);
        return;
     } 
     else if (p<=90.3){
        nucltransK(0.427,0.091,1.9e-1,0.);
        next405=true;
     } 
     else{
        nucltransK(0.065,0.016,6.9e+0,0.);
        next767=true; 
     }
  }
  else if (pbeta<=8.7999){
      beta(0.962,0.,0.);
      next405=true;
  }
  else{
     beta(1.367,0.,0.);
     return; 
  }

  if (next1014){
     p=100.*GetRandom();
     if (p<=28.7){
        nucltransK(1.014,0.091,1.9e-2,0.);
        return; 
     }
     else{
        nucltransK(0.609,0.091,7.0e-2,0.);
        next405=true;
     }
  }
  if (next767){
     p=100.*GetRandom();
     if (p<=57.4){
        nucltransK(0.767,0.091,4.0e-2,0.);
        return;
     } 
     else{
        nucltransK(0.362,0.091,2.8e-1,0.);
        next405=true;
     }
  }
  if (next405){
      fThlev=0.317e-9;
      nucltransK(0.405,0.091,1.3e-1,0.);
      return;
  }
}
//-----------------------------------------------

void Decay::Pb212(float tcnuc)
{
  // Scheme of Pb212 decay ("Table of Isotopes", 7th ed., 1978).
  // VIT, 5.08.1992, 22.10.1995.
  fThnuc=38304.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next239=false, next115=false; 
  float pbeta, p;
  pbeta=100.*GetRandom();  // beta decay to U234
  if (pbeta<=5.){
     beta(0.158,0.,0.);
     p=100.*GetRandom();
     if (p<=0.5){
        nucltransK(0.415,0.091,0.24,0.);
        return;
     }
     else if (p<=98.5){
        nucltransK(0.300,0.091,0.55,0.);
        next115=true;
     }
     else{
        nucltransK(0.177,0.091,2.4,0.);
        next239=true;
     }
  }
  else if (pbeta<=88.){
      beta(0.334,0.,0.);
      next239=true;
  }
  else{
     beta(0.573,0.,0.);
     return; 
  }
  
  if (next239){
     fThlev=1.e-12;
     nucltransK(0.239,0.091,1.1,0.);
     return; 
  }
  if (next115){
     fThlev=8.e-12;
     nucltransK(0.115,0.091,8.0,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Pb214(float tcnuc)
{
  // Scheme of 214Pb decay ("Table of Isotopes", 7th ed., 1978).
  // VIT, 11.08.1992, 22.10.1995.
  // VIT, 4.01.2007: updated to NDS 76(1995)127, 
  // conversion from K, L, M shells is introduced.
  fThnuc=1608.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next534=false, next377=false, next352=false, next295=false;
  bool next259=false, next63=false, next53=false;
  float pbeta, p;
  pbeta=100.*GetRandom();  // beta decay to U234
  if (pbeta<=0.034){
     beta(0.135,0.,0.);
     nucltransK(0.511,0.091,6.0e-2,0.);
     next377=true;
  }
  else if (pbeta<= 2.766){
     beta(0.184,0.,0.);
     p=100.*GetRandom();
     if (p<=21.3){
        nucltransKLM(0.839,0.091,2.99e-3,0.016,4.7e-4,0.004,1.5e-4,0.);
        return;
     }
     else if (p<=60.1){
        nucltransKLM(0.786,0.091,3.38e-3,0.016,5.3e-4,0.004,1.7e-4,0.);
        next53=true;
     }
     else if (p<=72.9){
        nucltransKLM(0.580,0.091,6.06e-3,0.016,9.7e-4,0.004,3.2e-4,0.);
        next259=true;
     }
     else if (p<=75.4){
        nucltransKLM(0.544,0.091,6.90e-3,0.016,1.11e-3,0.004,3.7e-4,0.);
        next295=true;
     }
     else if (p<=90.8){
        nucltransKLM(0.487,0.091,8.65e-3,0.016,1.41e-3,0.004,4.4e-4,0.);
        next352=true;
     }
     else if (p<=98.8){
        nucltransKLM(0.462,0.091,9.64e-3,0.016,1.58e-3,0.004,4.8e-4,0.);
        next377=true;
     }
     else{
        nucltransKLM(0.305,0.091,2.40e-2,0.016,4.1e-3,0.004,1.3e-3,0.);
        next534=true;
     }
  }
  else if (pbeta<= 2.787){
     beta(0.226,0.,0.);
     nucltransK(0.538,0.091,8.5e-3,0.);
     next259=true;
  }
  else if (pbeta<= 3.951){
     beta(0.489,0.,0.);
     next534=true;
  }
  else if (pbeta<=52.172){
     beta(0.671,0.,0.);
     next352=true;
  }
  else if (pbeta<=93.787){
     beta(0.728,0.,0.); 
     next295=true;
  }
  else{
     beta(1.023,0.,0.);
     return;
  }

  if (next534){
     p=100.*GetRandom();
     if (p<=16.3){
        nucltransKL(0.534,0.091,5.0e-2,0.016,1.0e-2,0.);
        return;
     }
     else if (p<=46.0){
        nucltransKLM(0.480,0.091,1.22e-1,0.016,1.9e-2,0.004,6.0e-3,0.);
        next53=true;
     }
     else{
        nucltransKLM(0.275,0.091,2.9e-1,0.016,7.3e-2,0.004,2.4e-2,0.);
        next259=true;
     }
  }
  if (next377){
     p=100.*GetRandom();
     if (p<=26.2){
        nucltransK(0.324,0.091,2.1e-1,0.);
        next53=true;
     } 
     else{
        nucltransK(0.314,0.091,2.3e-1,0.);
        next63=true;
     }
  }
  if (next352){
     nucltransKLM(0.352,0.091,2.55e-1,0.016,4.41e-2,0.004,1.38e-2,0.);
     return;
  }
  if (next295){
     p=100.*GetRandom();
     if (p<=67.10){
        nucltransKLM(0.295,0.091,3.8e-1,0.016,6.9e-2,0.004,2.2e-2,0.);
        return;
     }
     else{
        nucltransKLM(0.242,0.091,7.13e-1,0.016,1.23e-1,0.004,3.88e-2,0.);
        next53=true;
     }
  }
  if (next259){
     p=100.*GetRandom();
     if (p<=81.4){
        nucltransKLM(0.259,0.091,5.92e-1,0.016,1.03e-1,0.004,3.2e-2,0.);
        return;
     } 
     else if (p<=83.8){
        nucltransKLM(0.206,0.091,1.12e+0,0.016,1.95e-1,0.004,6.1e-2,0.);
        next53=true;
     }
     else{
       nucltransKLM(0.196,0.091,1.28e+0,0.016,2.23e-1,0.004,7.0e-2,0.);
       next63=true;
     }
  }
  if (next63){
     nucltransK(0.063,0.016,1.0e+1,0.);
     return; 
  }
  if (next53){
     nucltransKL(0.053,0.016,9.69e+0,0.004,3.05e+0,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Ra228(float tcnuc)
{
  // Scheme of Ra228 decay in accordance with NDS 80(1997)723 and
  // ENSDF database at NNDC site on 8.08.2007.
  fThnuc=1.814512e8;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  float pbeta, p;
  bool next202=false, next67=false;
  pbeta=100.*GetRandom();
  if (pbeta<=30.){
     beta(0.0128,0.,0.);
     p=100*GetRandom();
     if (p<=50.){
        nucltransK(0.0264,0.0198,2.1e2,0.);
        next67=true;
     }
     else{
        nucltransK(0.0128,0.0050,8.7e0,0.);
        next202=true;
     }
  }
  else if (pbeta<=50.){
      beta(0.0257,0.,0.); 
      next202=true;
  }
  else if (pbeta<=90.){
      beta(0.0392,0.,0.); 
      next67=true;
  }
  else{
     beta(0.0396,0.,0.);
     nucltransK(0.0063,0.0050,7.1e6,0.);
     return;
  }
  if (next202){
      nucltransK(0.0135,0.0050,6.1e0,0.);
      next67=true;
  }
  if (next67){
     nucltransK(0.0067,0.0050,1.6e6,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Rh106(float tcnuc)
{
  // Approximate scheme of 106Rh decay ("Table of Isotopes", 7th ed., 1978)
  // (beta decays to excited levels of 106Pd not higher 2.002 MeV, 99.32% of decay)
  // VIT, 17.12.1995.
  fThnuc=29.80;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next1562=false, next1134=false, next1128=false, next512=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<=0.45){
     beta(1.539,0.,0.);
     p=100*GetRandom();
     if (p<=0.4){
        nucltransK(1.490,0.024,4.0e-4,0.3e-4);
        next512=true;
     }
     else if (p<=98.6){
        nucltransK(0.874,0.024,1.3e-3,0.);
        next1128=true;
     } 
     else{
        nucltransK(0.440,0.024,8.5e-3,0.);
        next1562=true; 
     }
  }
  else if (pbeta<= 2.32){
      beta(1.979,0.,0.);
      next1562=true;
  }
  else if (pbeta<=12.32){
      beta(2.407,0.,0.);
      next1134=true;
  }
  else if (pbeta<=19.32){
      beta(3.029,0.,0.);
      next512=true;
  }
  else{
      beta(3.541,0.,0.);
      return;
  }

  if (next1562){
      p=100*GetRandom();
      if (p<=9.1){
         nucltransK(1.562,0.024,3.5e-4,0.4e-4);
         return; 
      } 
      else if (p<=95.6){
         nucltransK(1.051,0.024,8.5e-4,0.);
         next512=true;
      }
      else if (p<=96.8){
         nucltransK(0.434,0.024,8.5e-3,0.);
         next1128=true;
      }
      else{
         nucltransK(0.428,0.024,8.5e-3,0.);
         next1134=true;
      }
  }
  if (next1134){
      fThlev=7.0e-12;
      nucltransK(0.622,0.024,3.3e-3,0.);
      next512=true;
  }
  if (next1128){
     fThlev=3.2e-12;
     p=100*GetRandom();
     if (p<=34.){
        nucltransK(1.128,0.024,7.0e-4,0.);
        return;
     }
     else{
        nucltransK(0.616,0.024,3.0e-3,0.);
        next512=true;
     }
  }
  if (next512){
      fThlev=11.0e-12;
      nucltransK(0.512,0.024,5.5e-3,0.);
      return;
  }
}
//-----------------------------------------------

void Decay::Sb125(float tcnuc)
{
  // Scheme of 125Sb decay (NDS 86(1999)955 + NNDC on 7.02.2010).
  // VIT, 13.02.2010.
  fThnuc=8.705115e7;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next539=false, next525=false, next463=false, next444=false;
  bool next321=false, next145=false, next35=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<=13.446){
     beta(0.096,0.,0.);
     fThlev=1.26e-12;
     p=100*GetRandom();
     if (p<=13.347){
        nucltransK(0.671,0.032,3.8e-3,0.);
        return; 
     }
     else if (p<=96.921){
        nucltransK(0.636,0.032,5.3e-3,0.);
        next35=true;
     }
     else if (p<=97.980){
        nucltransK(0.228,0.032,8.4e-2,0.);
        next444=true;
     }
     else if (p<=99.993){
        nucltransK(0.208,0.032,9.1e-2,0.);
        next463=true;
     }
     else{
        nucltransK(0.133,0.032,3.1e-1,0.);
        next539=true;
     }
  }
  else if (pbeta<=13.502){
     beta(0.114,0.,0.); 
     p=100*GetRandom();
     if (p<=4.85){
        nucltransK(0.653,0.032,4.0e-3,0.);
        return;
     }
     else if (p<=14.56){
        nucltransK(0.617,0.032,5.6e-3,0.);
        next35=true; 
     }
     else if (p<=19.11){
        nucltransK(0.332,0.032,2.8e-2,0.);
        next321=true;
     }
     else{
        nucltransK(0.209,0.032,8.9e-2,0.);
        next444=true;
     }
  }
  else if (pbeta<=19.265){
     beta(0.125,0.,0.);
     fThlev=70e-12;
     p=100*GetRandom();
     if (p<=86.63){
        nucltransK(0.607,0.032,4.9e-3,0.);
        next35=true;
     }
     else if (p<=86.70){
        nucltransK(0.497,0.032,3.2e-2,0.);
        next145=true;
     }
     else if (p<=93.91){
        nucltransK(0.321,0.032,7.9e-3,0.);
        next321=true;
     }
     else if (p<=94.17){
        nucltransK(0.199,0.032,1.5e-1,0.);
        next444=true;
     }
     else if (p<=94.86){
        nucltransK(0.179,0.032,1.8e-1,0.);
        next463=true;
     }
     else{
        nucltransK(0.117,0.032,1.3e-1,0.);
        next525=true;
     }
  }
  else if (pbeta<=37.180){
     beta(0.131,0.,0.);
     fThlev=40e-12;
     p=100*GetRandom();
     if (p<=98.716){
        nucltransK(0.601,0.032,5.0e-3,0.);
        next35=true;
     }
     else if (p<=98.743){
        nucltransK(0.491,0.032,3.2e-2,0.);
        next145=true;
     }
     else if (p<=98.766){
        nucltransK(0.315,0.032,8.3e-3,0.);
        next321=true;
     }
     else if (p<=99.994){
        nucltransK(0.173,0.032,1.5e-1,0.);
        next463=true;
     }
     else{
        nucltransK(0.111,0.032,1.5e-1,0.);
        next525=true; 
     }
  }
  else if (pbeta<=38.791){
     beta(0.242,0.,0.);
     next525=true;
  }
  else if (pbeta<=79.199){
     beta(0.304,0.,0.);
     next463=true;
  }
  else if (pbeta<=79.247){
     beta(0.323,0.,0.);
     next444=true;
  }
  else if (pbeta<=79.268){
     beta(0.365,0.,0.);
     p=100*GetRandom();
     if (p<=29.72){
        nucltransK(0.402,0.032,1.9e-1,0.);
        return;
     }
     else if (p<=67.45){
        nucltransK(0.367,0.032,2.0e-2,0.);
        next35=true;
    }
    else{
       nucltransK(0.081,0.032,3.6e-1,0.);
       next321=true;
    }
  }
  else if (pbeta<=86.464){
     beta(0.446,0.,0.);
     next321=true;
  }
  else {
     beta(0.622,0.,0.);
     next145=true;
  }

  if (next539){
     fThlev=0.;
     p=100*GetRandom();
     if (p<=26.42){
        nucltransK(0.539,0.032,7.8e-3,0.);
        return;
     }
     else{
        nucltransK(0.503,0.032,9.3e-3,0.);
        next35=true;
     }
  }
  if (next525){ 
     fThlev=160e-12;
     p=100*GetRandom();
     if (p<=0.07){
        nucltransK(0.490,0.032,3.3e-2,0.);
        next35=true;
     }
     else if (p<=81.12){
        nucltransK(0.380,0.032,1.8e-2,0.);
        next145=true;
     } 
     else if (p<=99.89){
        nucltransK(0.204,0.032,1.3e-1,0.);
        next321=true;
     }
     else{
        nucltransK(0.062,0.032,7.4e-1,0.);
        next463=true;
     }
  }
  if (next463){
     fThlev=13.2e-12;
     p=100*GetRandom();
     if (p<=26.08){
        nucltransK(0.463,0.032,1.0e-2,0.);
        return;
     }
     else if (p<=99.39){
        nucltransK(0.428,0.032,1.4e-2,0.);
        next35=true;
     }
     else{
        nucltransK(0.020,0.005,11.1,0.);
        next444=true;
     }
  }
  if (next444){
     fThlev=19.1e-12;
     p=100*GetRandom();
     if (p<=62.40){
        nucltransK(0.444,0.032,1.2e-2,0.);
        return; 
     }
     else{
        nucltransK(0.408,0.032,1.5e-2,0.);
        next35=true;
     }
  }
  if (next321){
     fThlev=0.673e-9;
     nucltransK(0.176,0.032,1.7e-1,0.);
     next145=true;
  }
  if (next145){ 
     fThlev=4.959e6;
     nucltransK(0.109,0.032,355.,0.);
     next35=true;
  }
  if (next35){
     fThlev=1.48e-9;
     nucltransK(0.035,0.032,14.,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Sb126(float tcnuc)
{
  // Scheme of 126Sb decay (J.Katakura et al., NDS 97(2002)765).
  // VIT, 27.11.2003. Corrected 2.12.2003, thanks F.Capella.
  fThnuc=1.0670e6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next2840=false, next2812=false, next2766=false, next2515=false, next2497=false, next2396=false, next2218=false, next1776=false, next667=false, next1362=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<=0.48){
     beta(0.200,0.,0.);
     nucltransK(0.958,0.032,1.5e-3,0.);
     next2515=true;
  }
  else if (pbeta<= 2.48){
     beta(0.223,0.,0.);
     p=100*GetRandom();
     if (p<=57.14){
        nucltransK(0.954,0.032,1.5e-3,0.);
        next2497=true;
     }
     else{
        nucltransK(0.639,0.032,3.6e-3,0.);
        next2812=true;
     }
  }
  else if (pbeta<=30.28){
     beta(0.479,0.,0.);
     nucltransK(0.697,0.032,3.0e-3,0.);
     next2497=true;
  } 
  else if (pbeta<=35.94){ 
     beta(0.502,0.,0.); 
     p=100*GetRandom();
     if (p<=62.82){
        nucltransK(0.675,0.032,3.2e-3,0.);
        next2497=true;
     }
     else{
        nucltransK(0.656,0.032,3.4e-3,0.);
        next2515=true;
     }
  }
  else if (pbeta<=43.99){
     beta(0.602,0.,0.);
     p=100*GetRandom();
     if (p<=79.86){
        nucltransK(0.574,0.032,5.0e-3,0.);
        next2497=true; 
     }
     else{
        nucltransK(0.556,0.032,5.0e-3,0.);
        next2515=true;
     }
  }
  else if (pbeta<=48.02){
     beta(0.684,0.,0.);
     p=100*GetRandom();
     if (p<=57.18){
        nucltransK(1.213,0.032,1.0e-3,0.1e-4);
        next1776=true;
     }
     else if (p<=90.43){
        nucltransK(0.224,0.032,9.0e-2,0.);
        next2766=true;
     }
     else{
        nucltransK(0.149,0.032,4.0e-1,0.);
        next2840=true;
     }
  }
  else if (pbeta<=48.50){
     beta(0.699,0.,0.);
     nucltransK(0.209,0.032,1.3e-1,0.);
     next2766=true;
  }
  else if (pbeta<=49.27){
     beta(0.833,0.,0.);
     next2840=true;
  }
  else if (pbeta<=57.04){
     beta(0.861,0.,0.);
     next2812=true;
  }
  else if (pbeta<=61.74){
     beta(0.907,0.,0.); 
     next2766=true;
  }
  else if (pbeta<=77.08){
     beta(1.176,0.,0.);
     next2497=true;
  }
  else if (pbeta<=77.94){
     beta(1.277,0.,0.);
     next2396=true;
  }
  else if (pbeta<=80.82){
     beta(1.455,0.,0.);
     next2218=true;
  }
  else{
     beta(1.897,0.,0.);
     next1776=true;
  }
 
  if (next2840){
     p=100*GetRandom();
     if (p<=23.73){
        nucltransK(1.477,0.032,4.0e-4,0.3e-4);
        next1362=true;
     }
     else{
        nucltransK(1.064,0.032,1.2e-3,0.);
        next1776=true; 
     }
  }
  if (next2812){
     p=100*GetRandom();
     if (p<=83.33){
        nucltransK(0.593,0.032,5.0e-3,0.);
        next2218=true;
     }
     else if (p<=94.44){
        nucltransK(0.415,0.032,1.2e-2,0.);
        next2396=true; 
     }
     else{
        nucltransK(0.297,0.032,3.0e-2,0.);
        next2515=true;
     }
  }
  if (next2766){
     nucltransK(0.990,0.032,1.5e-3,0.);
     next1776=true; 
  }
  if (next2515){
     nucltransK(0.297,0.032,4.0e-2,0.);
     next2218=true; 
  }
  if (next2497){
     fThlev=0.152e-9;
     p=100*GetRandom();
     if (p<=95.56){
        nucltransK(0.721,0.032,1.0e-2,0.);
        next1776=true;
     }     
     else{ 
        nucltransK(0.278,0.032,4.9e-2,0.);
        next2218=true;
     }
  }
  if (next2396){
     p=100*GetRandom();
     if (p<=52.63){
        nucltransK(1.036,0.032,1.3e-3,0.);
        next1362=true;
     }
     else{
        nucltransK(0.620,0.032,5.0e-3,0.);
        next1776=true;
     }
  }
  if (next2218){
     fThlev=0.;
     nucltransK(0.857,0.032,8.4e-4,0.);
     next1362=true;
  } 
  if (next1776){
     fThlev=68.e-12;
     nucltransK(0.415,0.032,1.4e-2,0.);
     next1362=true;
  }
  if (next1362){
     fThlev=0.;
     nucltransK(0.695,0.032,3.4e-3,0.);
     next667=true;
  }
  if (next667){
     fThlev=0.;
     nucltransK(0.667,0.032,3.8e-3,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Sb133(float tcnuc)
{
  // Scheme of 133Sb decay (S.Rab, NDS 75(1995)491)
  // VIT, 11.12.2003. 
  fThnuc=150;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next2332=false, next2211=false, next2024=false, next1913=false;
  bool next1729=false, next1640=false, next1501=false, next1265=false;
  bool next1096=false, next334=false, next308=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<=27.84){
     beta(1.247,0.,0.);
     p=100*GetRandom();
     if (p<=45.15){
        nucltransK(2.755,0.032,1.5e-4,5.7e-4);
        return;
     }
     else if (p<=49.81){
        nucltransK(2.447,0.032,1.5e-4,5.2e-4);
        next308=true;
     }
     else if (p<=57.90){
        nucltransK(1.659,0.032,4.0e-4,1.3e-4);
        next1096=true;
     }
     else if (p<=67.22){
        nucltransK(1.490,0.032,5.0e-4,0.9e-4);
        next1265=true;
     }
     else if (p<=87.10){
        nucltransK(1.026,0.032,1.5e-3,0.);
        next1729=true; 
     }
     else{
        nucltransK(0.423,0.032,1.2e-2,0.);
        next2332=true;
     } 
  }
  else if (pbeta<=53.70){
     beta(1.253,0.,0.);
     p=100*GetRandom();
     if (p<=32.79){
        nucltransK(2.416,0.032,8.0e-5,5.0e-4);
        next334=true;
     }
     else if (p<=39.00){
        nucltransK(1.654,0.032,4.0e-4,1.3e-4);
        next1096=true; 
     }
     else if (p<=42.36){
        nucltransK(1.250,0.032,1.0e-3,0.4e-4);
        next1501=true;
     }
     else if (p<=49.24){
        nucltransK(1.111,0.032,1.5e-3,0.1e-4);
        next1640=true;
     }
     else if (p<=92.58){
        nucltransK(0.837,0.032,8.0e-4,0.);
        next1913=true;
     }
     else{
        nucltransK(0.539,0.032,7.0e-3,0.);
        next2211=true;
     }
  }
  else if (pbeta<=58.28){
     beta(1.792,0.,0.);
     next2211=true;
  }
  else if (pbeta<=58.58){
     beta(1.979,0.,0.);
     next2024=true;
  }
  else if (pbeta<=59.58){
     beta(2.027,0.,0.);
     nucltransK(1.976,0.032,3.0e-4,3.2e-4);
     return;
  }
  else if (pbeta<=68.54){
     beta(2.090,0.,0.);
     next1913=true;
  }
  else if (pbeta<=75.31){
     beta(2.274,0.,0.);
     next1729=true;
  }     
  else if (pbeta<=76.31){
     beta(2.297,0.,0.);
     nucltransK(1.706,0.032,4.5e-4,1.5e-4);
     return;
  }
  else if (pbeta<=80.19){
     beta(2.361,0.,0.);
     nucltransK(1.642,0.032,5.0e-4,1.2e-4);
     return;
  }
  else if (pbeta<=83.97){
     beta(2.363,0.,0.);
     next1640=true;
  }
  else if (pbeta<=88.15){
     beta(2.451,0.,0.);
     nucltransK(1.552,0.032,5.0e-4,1.0e-4);
     return;
  }
  else if (pbeta<=88.95){
     beta(2.502,0.,0.);
     nucltransK(0.404,0.032,1.3e-2,0.);
     next1096=true;
  }
  else if (pbeta<=91.14){
     beta(2.582,0.,0.);
     next1501=true;
  }
  else if (pbeta<=94.03){
     beta(2.738,0.,0.);
     next1265=true;
  }
  else{
     beta(2.907,0.,0.); 
     next1096=true;
  }

  if (next2332){
     p=100*GetRandom();
     if (p<=65.82){
        nucltransK(1.236,0.032,1.0e-3,0.4e-4);
        next1096=true;
     }
     else{
        nucltransK(0.308,0.032,2.5e-2,0.);
        next2024=true;
     }
  }
  if (next2211){ 
     p=100*GetRandom();
     if (p<=23.75){
        nucltransK(1.877,0.032,3.5e-4,1.9e-4);
        next334=true;
     }
     else if (p<=97.38){
        nucltransK(1.115,0.032,1.5e-3,0.1e-4);
        next1096=true;
     }
     else{
        nucltransK(0.572,0.032,5.5e-3,0.);
        next1640=true;
     }
  }
  if (next2024){
     p=100*GetRandom();
     if (p<=81.61){
        nucltransK(0.928,0.032,1.8e-3,0.);
        next1096=true;
     }
     else{
        nucltransK(0.523,0.032,7.0e-3,0.);
        next1501=true;
     }
  }
  if (next1913){
     p=100*GetRandom();
     if (p<=7.81){
        nucltransK(1.579,0.032,4.5e-4,0.4e-4);
        next334=true;
     }
     else if (p<=98.08){
        nucltransK(0.818,0.032,8.0e-4,0.);
        next1096=true;
     }
     else{
        nucltransK(0.413,0.032,1.2e-2,0.);
        next1501=true;
     }
  }
  if (next1729){
     p=100*GetRandom();
     if (p<=68.47){
        nucltransK(1.729,0.032,4.5e-4,1.5e-4);
        return;
     }
     else{
        nucltransK(0.632,0.032,4.5e-3,0.);
        next1096=true;
     }
  }
  if (next1640){ 
     nucltransK(1.305,0.032,9.0e-4,0.5e-4);
     next334=true;
  }
  if (next1501){ 
     p=100*GetRandom();
     if (p<=13.44){
        nucltransK(1.421,0.032,6.5e-4,0.7e-4);
        return;
     }
     else{
        nucltransK(1.113,0.032,1.3e-3,0.);
        next308=true;
     }
  }
  if (next1265){
     nucltransK(1.265,0.032,9.0e-4,0.5e-4);
     return;
  }
  if (next1096){
     nucltransK(1.096,0.032,1.1e-3,0.);
     return;
  }
  if (next334){
     return; // creation of isomeric 133mTe with E_exc=334 keV 
             // and T1/2=55.4 m
  }
  if (next308){
     nucltransK(0.308,0.032,2.5e-2,0.);
     return; // creation of 133Te in g.s. (T1/2=12.5 m)
  }
}
//-----------------------------------------------

void Decay::Ta182(float tcnuc)
{
  // Scheme of 182Ta decay ("Table of Isotopes", 7th ed., 1978).
  // VIT, 6.03.1996.
  fThnuc=9.936e6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next1488=false, next1443=false, next1374=false, next1331=false;
  bool next1289=false, next1257=false, next1221=false, next329=false;
  bool next100=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<=28.586){
     beta(0.258,0.,0.);
     fThlev=1.23e-9;
     p=100*GetRandom();
     if (p<=0.25){
        nucltransK(1.453,0.070,4.5e-3,0.1e-4);
        next100=true;
     }
     else if (p<=1.45){
        nucltransK(1.223,0.070,2.5e-3,0.1e-4);
        next329=true;
     }
     else if (p<=22.47){
        nucltransK(0.264,0.070,1.4e-1,0.);
        next1289=true;
     }
     else if (p<=65.50){
        nucltransK(0.222,0.070,5.0e-2,0.);
        next1331=true;
     }
     else if (p<=83.50){
        nucltransK(0.179,0.070,7.5e-1,0.);
        next1374=true;
     }
     else if (p<=84.00){
        nucltransK(0.110,0.070,3.0e-1,0.);
        next1443=true;
     }
     else{
        nucltransK(0.066,0.012,3.0e-0,0.);
        next1488=true;
     }
  }
  else if (pbeta<=28.716){
     beta(0.301,0.,0.);
     fThlev=0.; 
     p=100*GetRandom();
     if (p<=34.){
        nucltransK(1.410,0.070,2.4e-3,0.1e-4);
        next100=true;
     }
     else{
        nucltransK(1.181,0.070,1.4e-3,0.1e-4);
        next329=true;
     }
  }
  else if (pbeta<=31.416){
     beta(0.323,0.,0.);
     next1488=true;
  }
  else if (pbeta<=32.076){
     beta(0.368,0.,0.);
     next1443=true;
  }
  else if (pbeta<=52.066){
     beta(0.437,0.,0.);
     next1374=true;
  }
  else if (pbeta<=54.366){
     beta(0.480,0.,0.);
     next1331=true;
  }
  else if (pbeta<=94.346){
     beta(0.522,0.,0.);
     next1289=true;
  }
  else if (pbeta<=94.846){
     beta(0.554,0.,0.);
     next1257=true;
  }
  else if (pbeta<=99.846){
     beta(0.590,0.,0.);
     next1221=true;
  }
  else if (pbeta<=99.942){
     beta(1.482,0.,0.);
     next329=true;
  }
  else{
     beta(1.711,0.,0.);
     next100=true;
  }


  if (next1488){
     fThlev=0.; 
     p=100*GetRandom();
     if (p<=1.1){
        nucltransK(1.387,0.070,5.0e-3,0.5e-5);
        next100=true;
     }
     else if (p<=6.1){
        nucltransK(1.158,0.070,1.5e-3,0.5e-5);
        next329=true;
     }
     else if (p<=29.3){
        nucltransK(0.198,0.070,3.2e-1,0.);
        next1289=true;
     }
     else if (p<=70.8){
        nucltransK(0.156,0.070,1.2e-1,0.);
        next1331=true;
     }
     else{
        nucltransK(0.114,0.070,3.8e-0,0.);
        next1374=true;
     }
  }
  if (next1443){
     fThlev=0.; 
     p=100*GetRandom();
     if (p<=40.){
        nucltransK(1.343,0.070,2.8e-3,0.1e-4);
        next100=true;
     }
     else{
        nucltransK(1.113,0.070,6.0e-3,0.);
        next329=true;
     }
  }
  if (next1374){
     fThlev=0.08e-9;
     p=100*GetRandom();
     if (p<=2.0){
        nucltransK(1.374,0.070,5.5e-3,0.5e-5);
        return;
     } 
     else if (p<=7.7){
        nucltransK(1.274,0.070,3.0e-3,0.1e-4);
        next100=true;
     }
     else if (p<=9.8){
        nucltransK(1.044,0.070,6.6e-3,0.);
        next329=true;
     }
     else if (p<=71.3){
        nucltransK(0.152,0.070,1.2e-1,0.);
        next1221=true; 
     }
     else if (p<=75.1){
        nucltransK(0.116,0.070,2.6e-1,0.);
        next1257=true;
     }
     else if (p<=97.9){
        nucltransK(0.085,0.070,8.5e-0,0.);
        next1289=true; 
     } 
     else{
        nucltransK(0.043,0.012,7.0e-1,0.);
        next1331=true;
     }
  }
  if (next1331){ 
     fThlev=0.;
     p=100*GetRandom();
     if (p<=85.){
        nucltransK(1.231,0.070,3.0e-3,0.1e-4);
        next100=true;
     }
     else{
        nucltransK(1.002,0.070,4.7e-3,0.);
        next329=true;
     }
  }
  if (next1289){
     fThlev=1.12e-9;
     p=100*GetRandom();
     if (p<=2.35){
        nucltransK(1.289,0.070,1.3e-2,0.3e-4);
        return;
     }  
     else if (p<=29.75){
        nucltransK(1.189,0.070,5.3e-3,0.1e-4);
        next100=true;
     }
     else if (p<=30.34){
        nucltransK(0.960,0.070,1.3e-2,0.);
        next329=true;
     }
     else if (p<=99.00){
        nucltransK(0.068,0.012,2.0e-1,0.);
        next1221=true;
     } 
     else{
        nucltransK(0.032,0.012,1.6e-0,0.);
        next1257=true;
     }
  }
  if (next1257){ 
     fThlev=1.7e-12;
     p=100*GetRandom();
     if (p<=54.5){
        nucltransK(1.257,0.070,3.0e-3,0.1e-4);
        return;
     }
     else if (p<=77.8){
        nucltransK(1.157,0.070,5.3e-3,0.);
        next100=true;
     }
     else{
        nucltransK(0.928,0.070,5.5e-3,0.);
        next329=true; 
     }
  }
  if (next1221){
     fThlev=0.37e-12;
     p=100*GetRandom();
     if (p<=44.00){
        nucltransK(1.221,0.070,3.0e-3,0.1e-4);
        return;
     }
     else if (p<=99.92){
        nucltransK(1.121,0.070,3.5e-3,0.1e-5);
        next100=true;
     }
     else{
        nucltransK(0.892,0.070,6.0e-3,0.);
        next329=true;
     }
  }
  if (next329){ 
     fThlev=64.e-12;
     nucltransK(0.229,0.070,2.4e-1,0.);
     next100=true;
  }
  if (next100){
     fThlev=1.38e-9;
     nucltransK(0.100,0.012,4.0e-0,0.);
     return;  
  }
}
//-----------------------------------------------

void Decay::Te133(float tcnuc)
{
  // Model for scheme of Te133 decay(S.Rab, Nucl.Data Sheets 75(1995)491).
  fThnuc=750.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next2597=false, next2467=false, next2284=false, next2225=false;
  bool next2210=false, next2054=false, next2040=false, next2025=false;
  bool next1718=false, next1671=false, next1564=false, next1374=false;
  bool next1333=false, next1313=false, next1307=false, next1240=false;
  bool  next915=false,  next913=false,  next787=false,  next720=false;
  bool  next312=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<=0.13){
     beta(0.001,0.,0.);
     p=100*GetRandom();
     if (p<=72.87){
        nucltransK(2.624,0.033,1.5e-4,5.1e-4);
        next312=true;
     }
     else if (p<=92.25){
        nucltransK(2.148,0.033,3.0e-4,2.9e-4);
        next787=true;
     }
     else{
        nucltransK(1.372,0.033,1.0e-3,0.7e-4);
        next1564=true;
     }
  }
  else if (pbeta<= 0.60){
     beta(0.054,0.,0.);
     p=100*GetRandom();
     if (p<=74.63){
        nucltransK(2.554,0.033,1.8e-4,4.9e-4);
        next312=true;
     }
     else if (p<=94.67){
        nucltransK(2.079,0.033,3.2e-4,2.7e-4);
        next787=true;
     }
     else{
        nucltransK(1.493,0.033,7.5e-4,0.9e-4);
        next1374=true;
     }
  }
  else if (pbeta<= 0.86){
     beta(0.095,0.,0.);
     p=100*GetRandom();
     if (p<=60.94){
        nucltransK(2.825,0.033,1.3e-4,5.9e-4);
        return; 
     }
     else{
        nucltransK(2.106,0.033,3.2e-4,2.7e-4);
        next720=true;
     }
  }
  else if (pbeta<= 1.29){
     beta(0.112,0.,0.);
     p=100*GetRandom();
     if (p<=46.23){
        nucltransK(2.496,0.033,2.0e-4,4.7e-4);
        next312=true;
     }
     else if (p<=63.26){
        nucltransK(1.244,0.033,1.1e-3,0.4e-4);
        next1564=true;
     }
     else if (p<=92.46){
        nucltransK(1.137,0.033,1.5e-3,0.2e-4);
        next1671=true;
     }
     else{
        nucltransK(0.341,0.033,3.0e-2,0.);
        next2467=true;
     }
  }
  else if (pbeta<= 2.34){
     beta(0.152,0.,0.); 
     p=100*GetRandom();
     if (p<=24.57){
        nucltransK(2.456,0.033,2.0e-4,4.4e-4);
        next312=true;
     }
     else if (p<=28.07){
        nucltransK(2.049,0.033,1.2e-3,2.5e-4);
        next720=true;
     }
     else if (p<=42.25){
        nucltransK(1.455,0.033,1.0e-3,0.8e-4);
        next1313=true;
     }
     else if (p<=71.55){
        nucltransK(0.743,0.033,5.0e-3,0.);
        next2025=true;
     }
     else if (p<=82.89){
        nucltransK(0.544,0.033,2.0e-4,0.);
        next2225=true;
     }
     else if (p<=88.56){
        nucltransK(0.485,0.033,1.0e-2,0.);
        next2284=true;
     }
     else if (p<=91.49){
        nucltransK(0.302,0.033,4.5e-2,0.);
        next2467=true;
     }
     else{
        nucltransK(0.171,0.033,2.0e-1,0.);
        next2597=true;
     }
  }
  else if (pbeta<= 2.52){
     beta(0.259,0.,0.); 
     p=100*GetRandom();
     if (p<=43.10){
        nucltransK(2.661,0.033,1.5e-4,5.2e-4);
        return;
     }
     else if (p<=47.70){
        nucltransK(2.349,0.033,2.2e-4,3.9e-4);
        next312=true;
     }
     else if (p<=82.18){
        nucltransK(0.943,0.033,2.5e-3,0.);
        next1718=true;
     }
     else{
        nucltransK(0.620,0.033,6.5e-3,0.);
        next2040=true;
     }
  }
  else if (pbeta<= 2.70){
     beta(0.323,0.,0.);
     next2597=true;
  }
  else if (pbeta<= 4.86){
     beta(0.378,0.,0.);
     p=100*GetRandom();
     if (p<=23.06){
        nucltransK(2.542,0.033,1.6e-4,4.8e-4);
        return;
     }
     else if (p<=63.19){
        nucltransK(2.230,0.033,2.5e-4,3.2e-4);
        next312=true;
     }
     else if (p<=73.34){
        nucltransK(1.822,0.033,5.0e-4,1.8e-4);
        next720=true;
     }
     else if (p<=75.37){
        nucltransK(1.755,0.033,5.5e-4,1.6e-4);
        next787=true;
     }
     else if (p<=77.40){
        nucltransK(1.302,0.033,1.0e-3,0.5e-4);
        next1240=true;
     }
     else if (p<=84.78){
        nucltransK(1.209,0.033,1.4e-3,0.3e-4);
        next1333=true;
     }
     else if (p<=90.32){
        nucltransK(0.978,0.033,2.5e-3,0.);
        next1564=true;
     }
     else if (p<=94.47){
        nucltransK(0.488,0.033,1.2e-2,0.);
        next2054=true;
     }
     else{
        nucltransK(0.332,0.033,3.5e-2,0.);
        next2210=true; 
     }
  } 
  else if (pbeta<= 5.24){
     beta(0.394,0.,0.);
     p=100*GetRandom();
     if (p<=3.97){
        nucltransK(2.526,0.033,1.7e-4,5.7e-4);
        return;
     }
     else if (p<=37.30){
        nucltransK(2.214,0.033,2.5e-4,3.1e-4);
        next312=true;
     }
     else if (p<=78.57){
        nucltransK(1.807,0.033,5.0e-4,1.8e-4);
        next720=true;
     }
     else if (p<=86.51){
        nucltransK(1.738,0.033,5.5e-4,1.5e-4);
        next787=true;
     }
     else if (p<=90.48){
        nucltransK(1.286,0.033,1.0e-3,0.5e-4);
        next1240=true;
     }
     else{
        nucltransK(0.854,0.033,3.5e-3,0.);
        next1671=true;
     }
  } 
  else if (pbeta<= 5.60){
     beta(0.427,0.,0.);
     p=100*GetRandom();
     if (p<=8.93){
        nucltransK(2.181,0.033,2.5e-4,3.4e-4);
        next312=true;
     }
     else if (p<=48.72){
        nucltransK(1.773,0.033,5.0e-4,0.9e-4);
        next720=true;
     }
     else if (p<=65.81){
        nucltransK(1.706,0.033,6.0e-4,1.3e-4);
        next787=true;
     }
     else{
        nucltransK(0.928,0.033,2.8e-3,0.);
        next1564=true;
     }
  }
  else if (pbeta<= 6.20){
     beta(0.453,0.,0.);
     next2467=true;
  }
  else if (pbeta<= 6.75){
     beta(0.503,0.,0.); 
     p=100*GetRandom();
     if (p<=34.48){
        nucltransK(2.417,0.033,2.0e-4,4.2e-4);
        return;
     }
     else if (p<=52.63){
        nucltransK(2.106,0.033,3.0e-4,2.7e-4);
        next312=true;
     }
     else if (p<=63.52){
        nucltransK(1.697,0.033,5.0e-4,1.3e-4);
        next720=true;
     }
     else if (p<=68.06){
        nucltransK(1.630,0.033,5.2e-4,1.2e-4);
        next787=true;
     }
     else if (p<=73.69){
        nucltransK(1.503,0.033,7.0e-4,0.3e-4);
        next915=true;
     }
     else if (p<=95.46){
        nucltransK(1.110,0.033,1.8e-3,0.1e-4);
        next1307=true;
     }
     else{
        nucltransK(0.207,0.033,1.2e-1,0.);
        next2210=true;
     }
  }
  else if (pbeta<= 6.86){
     beta(0.527,0.,0.); 
     p=100*GetRandom();
     if (p<=11.32){
        nucltransK(2.393,0.033,2.0e-4,4.1e-4);
        return; 
     }
     else if (p<=64.16){
        nucltransK(2.081,0.033,3.0e-4,2.6e-4);
        next312=true;
     }
     else if (p<=82.08){
        nucltransK(0.722,0.033,5.0e-3,0.);
        next1671=true;
     }
     else{
        nucltransK(0.183,0.033,1.7e-1,0.);
        next2210=true; 
     }
  }
  else if (pbeta<= 7.00){
     beta(0.556,0.,0.); 
     p=100*GetRandom();
     if (p<=18.52){
        nucltransK(2.363,0.033,2.0e-4,4.0e-4);
        return;
     }
     else if (p<=62.96){
        nucltransK(1.124,0.033,1.8e-3,0.3e-4);
        next1240=true;
     }
     else{
        nucltransK(1.051,0.033,2.0e-3,0.);
        next1313=true;
     }
  }
  else if (pbeta<= 7.26){
     beta(0.636,0.,0.);
     next2284=true;
  }
  else if (pbeta<= 7.72){
     beta(0.654,0.,0.); 
     p=100*GetRandom();
     if (p<=52.17){
        nucltransK(2.266,0.033,2.5e-4,3.5e-4);
        return;
     }
     else if (p<=64.34){
        nucltransK(1.027,0.033,2.2e-3,0.);
        next1240=true;
     }
     else if (p<=90.43){
        nucltransK(0.934,0.033,2.8e-3,0.);
        next1333=true;
     }
     else{
        nucltransK(0.702,0.033,5.5e-3,0.);
        next1564=true;
     }
  }
  else if (pbeta<= 8.96){
     beta(0.665,0.,0.);
     p=100*GetRandom();
     if (p<=16.92){
        nucltransK(2.255,0.033,3.0e-4,3.8e-4);
        return;
     }
     else if (p<=23.45){
        nucltransK(1.944,0.033,3.5e-4,2.1e-4);
        next312=true;
     }
     else if (p<=40.37){
        nucltransK(1.535,0.033,6.0e-4,1.0e-4);
        next720=true; 
     }
     else if (p<=44.40){
        nucltransK(1.468,0.033,6.5e-4,0.9e-4);
        next787=true;
     }
     else if (p<=54.07){
        nucltransK(1.015,0.033,2.3e-3,0.);
        next1240=true;
     }
     else if (p<=79.05){
        nucltransK(0.942,0.033,2.7e-3,0.);
        next1313=true;
     }
     else if (p<=88.72){
        nucltransK(0.922,0.033,2.8e-3,0.);
        next1333=true;
     }
     else{
        nucltransK(0.691,0.033,5.5e-3,0.);
        next1564=true;
     }
  }
  else if (pbeta<= 9.95){
     beta(0.695,0.,0.);
     next2225=true;
  }
  else if (pbeta<=11.27){
     beta(0.710,0.,0.);
     next2210=true;
  }
  else if (pbeta<=13.53){
     beta(0.726,0.,0.);
     p=100*GetRandom();
     if (p<=25.07){
        nucltransK(2.194,0.033,3.0e-4,3.1e-4);
        return;
     }
     else if (p<=78.72){
        nucltransK(1.882,0.033,4.5e-4,1.9e-4);
        next312=true;
     }
     else if (p<=92.79){
        nucltransK(1.474,0.033,7.5e-4,0.9e-4);
        next720=true;
     }
     else if (p<=94.72){
        nucltransK(0.886,0.033,3.0e-3,0.);
        next1307=true;
     }
     else if (p<=97.36){
        nucltransK(0.881,0.033,3.0e-3,0.);
        next1313=true;
     }
     else{
        nucltransK(0.860,0.033,3.4e-3,0.);
        next1333=true;
     }
  }
  else if (pbeta<=15.90){
     beta(0.784,0.,0.);
     p=100*GetRandom();
     if (p<=57.73){
        nucltransK(2.136,0.033,3.0e-4,2.9e-4);
        return;
     }
     else if (p<=69.61){
        nucltransK(1.824,0.033,5.0e-4,1.7e-4);
        next312=true;
     } 
     else if (p<=75.10){
        nucltransK(1.417,0.033,1.0e-3,0.8e-4);
        next720=true;
     }
     else if (p<=79.32){
        nucltransK(1.350,0.033,1.2e-3,0.7e-4);
        next787=true;
     }
     else if (p<=80.12){
        nucltransK(1.222,0.033,1.0e-3,0.1e-4);
        next915=true;
     }
     else if (p<=82.23){
        nucltransK(0.897,0.033,3.0e-3,0.);
        next1240=true;
     }
     else if (p<=86.03){
        nucltransK(0.829,0.033,3.7e-3,0.);
        next1307=true;
     }
     else if (p<=89.41){
        nucltransK(0.824,0.033,3.7e-3,0.);
        next1313=true;
     }
     else if (p<=94.94){
        nucltransK(0.803,0.033,4.0e-3,0.);
        next1333=true;
     }
     else{
        nucltransK(0.763,0.033,3.2e-3,0.);
        next1374=true; 
     }
  }
  else if (pbeta<=16.97){
     beta(0.866,0.,0.);
     next2054=true;
  }
  else if (pbeta<=17.31){
     beta(0.880,0.,0.);
     next2040=true;
  }
  else if (pbeta<=18.70){
     beta(0.895,0.,0.);
     next2025=true;
  }
  else if (pbeta<=28.99){ 
     beta(1.202,0.,0.);
     next1718=true;
  }
  else if (pbeta<=30.31){
     beta(1.249,0.,0.);
     next1671=true;
  }
  else if (pbeta<=33.78){
     beta(1.356,0.,0.);
     next1564=true;
  }
  else if (pbeta<=34.97){
     beta(1.546,0.,0.);
     next1374=true;
  }
  else if (pbeta<=47.93){
     beta(1.587,0.,0.);
     next1333=true;
  }
  else if (pbeta<=51.20){
     beta(1.607,0.,0.);
     next1313=true;
  }
  else if (pbeta<=51.40){
     beta(1.680,0.,0.);
     next1240=true;
  }
  else if (pbeta<=79.41){
     beta(2.200,0.,0.);
     next720=true;
  }
  else{
     beta(2.608,0.,0.);
     next312=true;
  }
 
  if (next2597){
     p=100*GetRandom();
     if (p<=20.00){
        nucltransK(2.597,0.033,1.6e-4,5.1e-4);
        return;
     }
     else if (p<=23.21){
        nucltransK(2.286,0.033,2.3e-4,3.5e-4);
        next312=true;
     }
     else if (p<=73.21){
        nucltransK(1.683,0.033,5.0e-4,0.6e-4);
        next915=true;
     }
     else if (p<=80.00){
        nucltransK(1.290,0.033,1.0e-3,0.5e-4);
        next1307=true;
     }
     else if (p<=88.93){
        nucltransK(1.285,0.033,1.0e-3,0.5e-4);
        next1313=true;
     }
     else if (p<=91.07){
        nucltransK(1.224,0.033,1.0e-3,0.1e-4);
        next1374=true;
     }
     else{
        nucltransK(0.572,0.033,9.0e-3,0.);
        next2025=true;
     }
  }
  if (next2467){
     p=100*GetRandom();
     if (p<=61.38){
        nucltransK(2.467,0.033,1.8e-4,5.2e-4);
        return;
     }
     else if (p<=65.12){
        nucltransK(2.155,0.033,3.0e-4,3.0e-4);
        next312=true;
     }
     else if (p<=78.59){
        nucltransK(1.680,0.033,5.0e-4,1.3e-4);
        next787=true;
     }
     else if (p<=95.36){
        nucltransK(1.228,0.033,1.2e-3,0.4e-4);
        next1240=true;
     }
     else{
        nucltransK(0.242,0.033,8.0e-2,0.);
        next2225=true;
     }
  }
  if (next2284){ 
     p=100*GetRandom();
     if (p<=5.77){
        nucltransK(1.972,0.033,3.5e-4,1.7e-4);
        next312=true;
     }
     else if (p<=34.62){
        nucltransK(1.564,0.033,6.0e-4,0.4e-4);
        next720=true;
     }
     else if (p<=53.85){
        nucltransK(0.971,0.033,2.5e-3,0.);
        next1313=true;
     }
     else{
        nucltransK(0.910,0.033,3.0e-3,0.);
        next1374=true;
     }
  }
  if (next2225){ 
     p=100*GetRandom();
     if (p<=19.00){
        nucltransK(2.225,0.033,3.0e-4,3.2e-4);
        return;
     }
     else if (p<=29.28){
        nucltransK(1.913,0.033,4.5e-4,2.0e-4);
        next312=true;
     }
     else if (p<=35.32){
        nucltransK(1.505,0.033,7.5e-4,0.9e-4);
        next720=true; 
     }
     else if (p<=35.84){
        nucltransK(1.438,0.033,8.5e-4,0.8e-4);
        next787=true;
     }
     else if (p<=47.67){
        nucltransK(1.310,0.033,9.0e-4,0.1e-4);
        next915=true;
     }
     else if (p<=53.02){
        nucltransK(0.912,0.033,3.0e-3,0.);
        next1313=true;
     }
     else if (p<=79.79){
        nucltransK(0.852,0.033,2.5e-3,0.);
        next1374=true;
     }
     else if (p<=84.97){
        nucltransK(0.554,0.033,9.5e-3,0.);
        next1671=true;
     }
     else if (p<=96.80){
        nucltransK(0.507,0.033,1.2e-2,0.);
        next1718=true;
     }
     else{
        nucltransK(0.200,0.033,1.4e-1,0.);
        next2025=true;
     }
  }
  if (next2210){
     p=100*GetRandom();
     if (p<=46.12){
        nucltransK(2.210,0.033,2.9e-4,3.6e-4);
        return;
     }
     else if (p<=53.21){
        nucltransK(1.898,0.033,4.5e-4,2.0e-4);
        next312=true;
     }
     else if (p<=61.23){
        nucltransK(1.490,0.033,7.5e-4,0.9e-4);
        next720=true;
     }
     else if (p<=73.93){
        nucltransK(0.903,0.033,3.0e-3,0.);
        next1307=true;
     }
     else{
        nucltransK(0.646,0.033,6.5e-3,0.);
        next1564=true;
     }
  }
  if (next2054){  
     p=100*GetRandom();
     if (p<=12.21){
        nucltransK(2.054,0.033,4.0e-4,2.5e-4);
        return;
     }
     else if (p<=24.15){
        nucltransK(1.742,0.033,5.0e-4,1.5e-4);
        next312=true;
     }
     else if (p<=31.12){
        nucltransK(1.334,0.033,1.2e-3,0.6e-4);
        next720=true;
     }
     else if (p<=47.69){
        nucltransK(1.267,0.033,1.3e-3,0.5e-4);
        next787=true;
     }
     else if (p<=58.15){
        nucltransK(0.813,0.033,3.8e-3,0.);
        next1240=true; 
     }
     else if (p<=67.74){
        nucltransK(0.746,0.033,4.5e-3,0.);
        next1307=true;
     }
     else if (p<=84.31){
        nucltransK(0.741,0.033,4.5e-3,0.);
        next1313=true; 
     }
     else if (p<=94.77){
        nucltransK(0.720,0.033,5.0e-3,0.);
        next1333=true;
     }
     else{
        nucltransK(0.680,0.033,4.0e-3,0.);
        next1374=true;
     }
  }
  if (next2040){
     p=100*GetRandom();
     if (p<=6.87){
        nucltransK(1.320,0.033,9.5e-4,0.1e-4);
        next720=true;
     }
     else if (p<=12.09){
        nucltransK(1.254,0.033,1.2e-3,0.5e-4);
        next787=true;
     }
     else if (p<=23.08){
        nucltransK(0.727,0.033,5.0e-3,0.);
        next1313=true;
     }
     else if (p<=75.28){
        nucltransK(0.667,0.033,6.0e-3,0.);
        next1374=true;
     }
     else{
        nucltransK(0.369,0.033,2.5e-2,0.);
        next1671=true;
     }
  }
  if (next2025){
     p=100*GetRandom();
     if (p<=4.60){
        nucltransK(2.025,0.033,3.5e-4,2.2e-4);
        return;
     }
     else if (p<=25.61){
        nucltransK(1.713,0.033,5.0e-4,1.4e-4);
        next312=true;
     }
     else if (p<=40.94){
        nucltransK(1.306,0.033,1.2e-3,0.6e-4);
        next720=true;
     }
     else if (p<=47.75){
        nucltransK(1.239,0.033,1.3e-3,0.5e-4);
        next787=true;
     }
     else if (p<=54.00){
        nucltransK(0.718,0.033,1.4e-3,0.);
        next1307=true;
     }
     else if (p<=64.79){
        nucltransK(0.713,0.033,5.0e-3,0.);
        next1313=true; 
     }
     else{
        nucltransK(0.461,0.033,1.5e-2,0.);
        next1564=true;
     }
  } 
  if (next1718){ 
     p=100*GetRandom();
     if (p<=30.02){
        nucltransK(1.718,0.033,6.0e-4,1.4e-4);
        return;
     }
     else if (p<=35.59){
        nucltransK(1.406,0.033,1.2e-3,0.7e-4);
        next312=true;
     }
     else if (p<=45.41){
        nucltransK(0.998,0.033,2.4e-3,0.);
        next720=true;
     }
     else if (p<=81.39){
        nucltransK(0.931,0.033,2.7e-3,0.);
        next787=true;
     }
     else if (p<=81.97){
        nucltransK(0.803,0.033,1.1e-3,0.);
        next915=true;
     }
     else if (p<=85.56){
        nucltransK(0.478,0.033,1.3e-2,0.);
        next1240=true;
     }
     else if (p<=94.43){
        nucltransK(0.410,0.033,2.0e-2,0.);
        next1307=true; 
     }
     else if (p<=96.88){
        nucltransK(0.405,0.033,2.0e-2,0.);
        next1313=true;
     }
     else if (p<=99.43){
        nucltransK(0.384,0.033,2.5e-2,0.);
        next1333=true; 
     }
     else{
        nucltransK(0.344,0.033,3.0e-2,0.);
        next1374=true;
     }
  }
  if (next1671){
     p=100*GetRandom();
     if (p<= 9.59){
        nucltransK(1.671,0.033,5.0e-4,0.7e-4);
        return;
     }
     else if (p<=15.22){
        nucltransK(1.359,0.033,1.2e-3,0.7e-4);
        next312=true;
     }
     else if (p<=28.40){
        nucltransK(0.952,0.033,2.6e-3,0.);
        next720=true;
     }
     else if (p<=71.54){
        nucltransK(0.884,0.033,3.2e-3,0.);
        next787=true;
     }
     else if (p<=78.73){
        nucltransK(0.432,0.033,1.8e-2,0.);
        next1240=true;
     }
     else if (p<=83.94){
        nucltransK(0.359,0.033,3.0e-2,0.);
        next1313=true;
     }
     else{
        nucltransK(0.338,0.033,3.5e-2,0.);
        next1333=true;
     }
  }
  if (next1564){
     p=100*GetRandom();
     if (p<=28.18){
        nucltransK(1.252,0.033,1.3e-3,0.5e-4);
        next312=true;
     }
     else if (p<=92.95){
        nucltransK(0.844,0.033,3.5e-3,0.);
        next720=true; 
     }
     else if (p<=96.86){
        nucltransK(0.778,0.033,4.2e-3,0.);
        next787=true;
     }
     else if (p<=97.85){
        nucltransK(0.324,0.033,3.7e-2,0.);
        next1240=true;
     }
     else if (p<=98.46){
        nucltransK(0.251,0.033,7.5e-2,0.);
        next1313=true;
     } 
     else if (p<=98.83){
        nucltransK(0.231,0.033,8.5e-2,0.);
        next1333=true;
     }
     else{
        nucltransK(0.191,0.033,1.5e-1,0.);
        next1374=true;
     }
  }
  if (next1374){
     p=100*GetRandom();
     if (p<=56.13){
        nucltransK(1.062,0.033,1.5e-3,0.);
        next312=true; 
     }
     else if (p<=14.62){
        nucltransK(0.654,0.033,4.5e-3,0.);
        next720=true;  
     }
     else{
        nucltransK(0.587,0.033,8.5e-3,0.);
        next787=true;
     } 
  }
  if (next1333){
     p=100*GetRandom();
     if (p<=74.49){
        nucltransK(1.333,0.033,1.2e-3,0.6e-4);
        return;
     }
     else if (p<=94.05){
        nucltransK(1.021,0.033,2.3e-3,0.);
        next312=true;
     }
     else if (p<=96.28){
        nucltransK(0.614,0.033,7.0e-3,0.);
        next720=true;
     }
     else if (p<=99.83){
        nucltransK(0.546,0.033,9.5e-3,0.);
        next787=true;
     }
     else{
        nucltransK(0.418,0.033,1.5e-2,0.);
        next915=true;
     }
  }
  if (next1313){
     p=100*GetRandom();
     if (p<=17.45){
        nucltransK(1.313,0.033,9.0e-4,0.);
        return;   
     }
     else if (p<=91.78){
        nucltransK(1.001,0.033,2.5e-3,0.);
        next312=true;
     }
     else if (p<=95.27){
        nucltransK(0.593,0.033,8.5e-3,0.);
        next720=true;
     }
     else{
        nucltransK(0.526,0.033,1.2e-2,0.);
        next787=true;
     }
  }
  if (next1307){
     p=100*GetRandom();
     if (p<=33.13){
        nucltransK(1.307,0.033,1.2e-3,0.5e-4);
        return;
     }
     else if (p<=75.46){
        nucltransK(0.995,0.033,2.5e-3,0.);
        next312=true;
     }
     else if (p<=81.59){
        nucltransK(0.588,0.033,8.5e-3,0.);
        next720=true;
     }
     else if (p<=82.76){
        nucltransK(0.520,0.033,8.5e-3,0.);
        next787=true;
     }
     else if (p<=84.66){
        nucltransK(0.394,0.033,2.0e-2,0.);
        next913=true;
     }
     else{
        nucltransK(0.392,0.033,2.2e-2,0.);
        next915=true;
     }
  }
  if (next1240){
     p=100*GetRandom();
     if (p<=24.03){
        nucltransK(1.240,0.033,1.2e-3,0.);
        return;
     }
     else if (p<=81.24){
        nucltransK(0.928,0.033,2.7e-3,0.);
        next312=true;
     }
     else if (p<=86.27){
        nucltransK(0.520,0.033,1.2e-2,0.);
        next720=true;
     }
     else{
        nucltransK(0.453,0.033,1.5e-2,0.);
        next787=true;
     }
  }
  if (next915){
     nucltransK(0.915,0.033,2.8e-3,0.);
     return;
  } 
  if (next913){
     nucltransK(0.913,0.033,2.0e-3,0.);
     return;
  }
  if (next787){ 
     p=100*GetRandom();
     if (p<=78.49){
        nucltransK(0.787,0.033,3.0e-3,0.);
        return;
     }
     else if (p<=91.28){
        nucltransK(0.475,0.033,1.2e-2,0.);
        next312=true;
     }
     else{
        nucltransK(0.067,0.033,4.9,0.);
        next720=true;
     }
  } 
  if (next720){
     p=100*GetRandom();
     if (p<=24.72){
        nucltransK(0.720,0.033,5.0e-3,0.);
        return;
     }
     else{
        nucltransK(0.408,0.033,2.0e-2,0.);
        next312=true;
     }
  } 
  if (next312){
     nucltransK(0.312,0.033,3.5e-2,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Te133m(float tcnuc)
{
  // Model for scheme of Te133m decay(S.Rab,Nucl.Data Sheets 75(1995)491; 
  // E_exc=334 keV).
  fThnuc=3324.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next2596=false, next2556=false, next2516=false, next2500=false;
  bool next2445=false, next2372=false, next2262=false, next2249=false;
  bool next2212=false, next2049=false, next2005=false, next1991=false;
  bool next1975=false, next1943=false, next1893=false, next1886=false;
  bool next1817=false, next1799=false, next1797=false, next1777=false;
  bool next1729=false, next1707=false, next1704=false, next1647=false;
  bool next1634=false, next1560=false, next1516=false, next1455=false;
  bool next1307=false,  next915=false,  next913=false,  next312=false; 
  float pdecay, pbeta, p;
  pdecay=100.*GetRandom(); 
  if (pdecay<=17.5){  // 17.5% IT to 133Te(g.s.)
     nucltransK(0.334,0.033,1.431,0.);
     return;
  }
  else{               // 82.5% beta decay to 133I
     pbeta=100.*GetRandom(); 
     if (pbeta<=1.77){
        beta(0.203,0.,0.); 
        p=100*GetRandom();
        if (p<=21.62){
           nucltransK(3.051,0.033,2.0e-5,18.1e-4);
           return;
        } 
        else if (p<=29.05){
           nucltransK(1.405,0.033,3.5e-4,0.3e-4);
           next1647=true;
        }
        else{
           nucltransK(0.535,0.033,8.0e-3,0.);
           next2516=true;
        }
     }
     else if (pbeta<= 2.70){
        beta(0.225,0.,0.);
        p=100*GetRandom();
        if (p<=35.90){
            nucltransK(1.574,0.033,2.6e-4,0.7e-4);
            next1455=true;
         }
         else if (p<=78.21){
            nucltransK(1.252,0.033,4.0e-4,0.2e-4);
            next1777=true;
         }
         else{
            nucltransK(1.054,0.033,7.0e-4,0.);
            next1975=true; 
         }
     }
     else if (pbeta<= 3.03){
        beta(0.279,0.,0.);
        p=100*GetRandom();
        if (p<=31.25){
           nucltransK(2.062,0.033,1.4e-4,4.8e-4);
           next913=true;
        }
        else{
           nucltransK(1.198,0.033,5.0e-4,0.1e-4);
           next1777=true;
        }
     }
     else if (pbeta<= 3.16){
        beta(0.286,0.,0.); 
        nucltransK(2.968,0.033,8.0e-5,15.1e-4);
        return;
     }
     else if (pbeta<= 4.25){
        beta(0.373,0.,0.); 
        p=100*GetRandom();
        if (p<=17.89){
           nucltransK(1.968,0.033,1.7e-4,3.2e-4);
           next913=true;
        }
        else if (p<=47.36){
           nucltransK(0.890,0.033,3.1e-3,0.);
           next1991=true;
        }
        else if (p<=76.83){
           nucltransK(0.632,0.033,7.0e-3,0.);
           next2249=true;
        } 
        else{
           nucltransK(0.285,0.033,5.0e-2,0.);
           next2596=true;
        }
     }
     else if (pbeta<= 5.22){
        beta(0.428,0.,0.);
        p=100*GetRandom();
        if (p<=16.10){
           nucltransK(2.826,0.033,9.0e-5,13.7e-4);
           return;
        }
        else if (p<=23.00){
           nucltransK(1.914,0.033,1.8e-4,2.7e-4);
           next913=true;
        }
        else if (p<=55.18){
           nucltransK(1.372,0.033,3.5e-4,0.3e-4);
           next1455=true;
        }
        else if (p<=67.82){
           nucltransK(0.852,0.033,1.0e-3,0.);
           next1975=true;
        }
        else{
           nucltransK(0.326,0.033,3.6e-2,0.);
           next2500=true;
        }
     }
     else if (pbeta<= 6.02){
        beta(0.446,0.,0.);
        p=100*GetRandom();
        if (p<=16.42){
           nucltransK(1.104,0.033,6.0e-4,0.);
           next1704=true;
        }
        else if (p<=58.21){
           nucltransK(0.308,0.033,1.0e-2,0.);
           next2500=true;
        }
        else{
           nucltransK(0.251,0.033,7.5e-2,0.);
           next2556=true;
        }
     } 
     else if (pbeta<= 7.59){
        beta(0.470,0.,0.);
        p=100*GetRandom();
        if (p<=42.97){ 
           nucltransK(1.871,0.033,1.8e-4,2.2e-4);
           next913=true;
        }
        else if (p<=94.53){
           nucltransK(1.008,0.033,7.0e-4,0.);
           next1777=true;
        }
        else{
           nucltransK(0.734,0.033,1.4e-3,0.);
           next2049=true;
        }
     }
     else if (pbeta<=10.97){
        beta(0.568,0.,0.);
        p=100*GetRandom();
        if (p<=23.49){
           nucltransK(1.773,0.033,1.9e-4,1.4e-4);
           next913=true;
        }
        else if (p<=62.64){
           nucltransK(0.801,0.033,1.1e-3,0.);
           next1886=true;
        }
        else if (p<=70.47){
           nucltransK(0.637,0.033,1.7e-3,0.);
           next2049=true;
        }
        else if (p<=74.38){
           nucltransK(0.475,0.033,3.5e-3,0.);
           next2212=true;
        }
        else if (p<=88.26){
           nucltransK(0.314,0.033,1.0e-2,0.);
           next2372=true;
        }
        else{
           nucltransK(0.241,0.033,2.0e-2,0.);
           next2445=true;
        }
     }
     else if (pbeta<=31.35){
        beta(0.658,0.,0.);
        next2596=true;
     }
     else if (pbeta<=36.78){
        beta(0.698,0.,0.);
        next2556=true;
     }
     else if (pbeta<=37.87){
        beta(0.738,0.,0.); 
        next2516=true;
     }
     else if (pbeta<=39.00){
        beta(0.748,0.,0.);
        p=100*GetRandom();
        if (p<=64.89){
           nucltransK(0.945,0.033,2.6e-3,0.);
           next1560=true;
        } 
        else{
           nucltransK(0.244,0.033,7.5e-2,0.);
           next2262=true;
        }
     }
     else if (pbeta<=41.05){
        beta(0.754,0.,0.);
        next2500=true;
     }
     else if (pbeta<=41.77){
        beta(0.771,0.,0.); 
        p=100*GetRandom();
        if (p<=11.66){
           nucltransK(2.483,0.033,2.2e-4,4.5e-4);
           return;
        }
        else if (p<=31.09){
           nucltransK(1.570,0.033,7.0e-4,1.1e-4);
           next913=true;
        } 
        else{
           nucltransK(1.174,0.033,1.5e-3,0.3e-4);
           next1307=true;
        }
     }
     else if (pbeta<=45.63){
        beta(0.787,0.,0.); 
        p=100*GetRandom();
        if (p<=5.30){
           nucltransK(1.552,0.033,7.0e-4,1.1e-4);
           next915=true;
        }
        else if (p<=27.73){
           nucltransK(0.574,0.033,8.5e-3,0.);
           next1893=true;
        }
        else if (p<=51.72){
           nucltransK(0.493,0.033,1.2e-2,0.);
           next1975=true; 
        }
        else{
           nucltransK(0.462,0.033,1.4e-2,0.);
           next2005=true;
        }
     }
     else if (pbeta<=46.96){
        beta(0.809,0.,0.);
        next2445=true;
     }
     else if (pbeta<=47.36){
        beta(0.827,0.,0.);
        nucltransK(0.178,0.033,1.7e-1,0.);
        next2249=true;
     }
     else if (pbeta<=48.10){
        beta(0.835,0.,0.);
        p=100*GetRandom();
        if (p<= 45.91){
           nucltransK(1.506,0.033,9.0e-4,0.9e-4);
           next913=true; 
        }
        else if (p<=63.94){
           nucltransK(0.859,0.033,2.9e-3,0.);
           next1560=true;
        } 
        else if (p<=81.97){
           nucltransK(0.415,0.033,2.0e-2,0.);
           next2005=true;           
        }
        else{
           nucltransK(0.158,0.033,2.5e-1,0.);
           next2262=true;
        }
     }
     else if (pbeta<=55.22){
        beta(0.882,0.,0.);
        next2372=true;
     }
     else if (pbeta<=56.31){
        beta(0.992,0.,0.);
        next2262=true;
     }
     else if (pbeta<=57.64){
        beta(1.005,0.,0.);
        next2249=true;
     }
     else if (pbeta<=59.21){
        beta(1.042,0.,0.);
        next2212=true;
     }
     else if (pbeta<=61.14){
        beta(1.112,0.,0.);
        p=100*GetRandom();
        if (p<=12.99){
           nucltransK(1.230,0.033,1.5e-3,0.4e-4);
           next913=true;
        }
        else if (p<=23.03){
           nucltransK(1.228,0.033,1.5e-3,0.4e-4);
           next915=true;
        }
        else if (p<=24.98){
           nucltransK(0.249,0.033,7.5e-2,0.);
           next1893=true;
        }
        else if (p<=63.96){
           nucltransK(0.151,0.033,1.8e-2,0.);
           next1991=true;
        }
        else if (p<=72.82){
           nucltransK(0.137,0.033,4.0e-1,0.);
           next2005=true;
        }
        else{
           nucltransK(0.092,0.033,1.3,0.);
           next2049=true;
        }
     }
     else if (pbeta<=63.07){
        beta(1.205,0.,0.);
        next2049=true;
     }
     else if (pbeta<=66.21){
        beta(1.249,0.,0.);
        next2005=true;
     }
     else if (pbeta<=72.72){
        beta(1.263,0.,0.); 
        next1991=true;
     }
     else if (pbeta<=73.32){
        beta(1.279,0.,0.);
        next1975=true;
     }
     else if (pbeta<=74.41){
        beta(1.311,0.,0.); 
        next1943=true;
     }
     else if (pbeta<=78.03){
        beta(1.368,0.,0.);
        next1886=true;
     }
     else if (pbeta<=86.59){
        beta(1.477,0.,0.);
        next1777=true;
     }
     else if (pbeta<=88.28){
        beta(1.607,0.,0.);
        next1647=true;
     }
     else if (pbeta<=91.06){
        beta(1.694,0.,0.);
        next1560=true;
     }
     else if (pbeta<=91.66){
        beta(1.738,0.,0.);
        next1516=true;
     }
     else if (pbeta<=92.87){
        beta(1.799,0.,0.);
        next1455=true;
     }
     else if (pbeta<=93.96){
        beta(2.339,0.,0.);
        next915=true;
     }
     else if (pbeta<=95.17){
        beta(2.341,0.,0.);
        next913=true;
     }
     else{
        beta(3.254,0.,0.); 
        return;
     }
    
     if (next2596){ 
        p=100*GetRandom();
        if (p<=24.04){
           nucltransK(1.683,0.033,2.5e-4,1.0e-4);
           next913=true;
        }
        else if (p<=27.23){
           nucltransK(1.080,0.033,6.0e-4,0.);
           next1516=true;
        }
        else if (p<=31.05){
           nucltransK(0.949,0.033,8.0e-4,0.);
           next1647=true;
        }
        else if (p<=37.13){
           nucltransK(0.891,0.033,9.0e-4,0.);
           next1704=true; 
        }
        else if (p<=41.94){
           nucltransK(0.889,0.033,9.0e-4,0.);
           next1707=true;
        }
        else if (p<=42.92){
           nucltransK(0.819,0.033,3.7e-3,0.);
           next1777=true; 
        }
        else if (p<=47.09){
           nucltransK(0.710,0.033,1.4e-3,0.);
           next1886=true;
        }
        else if (p<=61.18){
           nucltransK(0.703,0.033,1.4e-3,0.);
           next1893=true;
        }
        else if (p<=64.71){
           nucltransK(0.653,0.033,1.7e-3,0.);
           next1943=true;
        }
        else if (p<=67.61){
           nucltransK(0.621,0.033,1.8e-3,0.);
           next1975=true;
        }
        else if (p<=74.97){
           nucltransK(0.605,0.033,7.7e-3,0.);
           next1991=true;
        }
        else if (p<=75.95){
           nucltransK(0.384,0.033,2.5e-2,0.);
           next2212=true;
        }
        else if (p<=79.77){
           nucltransK(0.347,0.033,3.0e-2,0.);
           next2249=true;
        }
        else if (p<=99.02){
           nucltransK(0.334,0.033,3.4e-2,0.); 
           next2262=true;
        }
        else{
           nucltransK(0.224,0.033,1.0e-1,0.);
           next2372=true;
        }
     }
     if (next2556){
        p=100*GetRandom();
        if (p<=6.82){
           nucltransK(1.644,0.033,2.6e-4,0.9e-4);
           next913=true;
        }
        else if (p<=15.08){
           nucltransK(0.996,0.033,7.0e-4,0.);
           next1560=true;
        }
        else if (p<=26.44){
           nucltransK(0.827,0.033,3.7e-3,0.);
           next1729=true;
        }
        else if (p<=63.01){
           nucltransK(0.780,0.033,4.2e-3,0.);
           next1777=true;
        }
        else if (p<=65.28){
           nucltransK(0.663,0.033,6.0e-3,0.);
           next1893=true;
        }
        else if (p<=75.61){
           nucltransK(0.581,0.033,8.2e-3,0.);
           next1975=true;
        }
        else if (p<=77.06){
           nucltransK(0.565,0.033,8.3e-3,0.);
           next1991=true;
        }
        else if (p<=91.94){
           nucltransK(0.344,0.033,3.2e-2,0.);
           next2212=true;
        }
        else if (p<=96.49){
           nucltransK(0.295,0.033,4.8e-2,0.);
           next2262=true;
        }
        else{
           nucltransK(0.185,0.033,1.6e-1,0.);
           next2372=true;
        }
     }
     if (next2516){
        p=100*GetRandom();
        if (p<=30.50){
           nucltransK(0.740,0.033,4.8e-3,0.);
           next1777=true;
        }
        else if (p<=72.00){
           nucltransK(0.719,0.033,4.9e-3,0.);
           next1797=true;
        }
        else if (p<=86.00){
           nucltransK(0.623,0.033,6.0e-3,0.);
           next1893=true;
        }
        else{
           nucltransK(0.526,0.033,1.1e-2,0.);
           next1991=true;
        }
     }
     if (next2500){
        p=100*GetRandom();
        if (p<=62.53){
           nucltransK(1.588,0.033,7.0e-4,1.1e-4);
           next913=true;
        }
        else if (p<=67.31){
           nucltransK(0.796,0.033,1.1e-3,0.);
           next1704=true;
        }
        else if (p<=72.09){
           nucltransK(0.793,0.033,1.1e-3,0.);
           next1707=true;
        }
        else if (p<=84.25){
           nucltransK(0.724,0.033,5.0e-3,0.);
           next1777=true;
        }
        else if (p<=91.63){
           nucltransK(0.607,0.033,2.0e-3,0.);
           next1893=true;
        }
        else{
           nucltransK(0.495,0.033,3.0e-3,0.);
           next2005=true;
        }
     }
     if (next2445){
        p=100*GetRandom();
        if (p<=75.00){
           nucltransK(0.885,0.033,9.5e-4,0.);
           next1560=true;
        }
        else{
           nucltransK(0.629,0.033,7.0e-3,0.);
           next1817=true;
        }
     }
     if (next2372){
        p=100*GetRandom();
        if (p<= 2.55){
           nucltransK(1.459,0.033,2.5e-4,0.5e-4);
           next913=true;
        } 
        else if (p<= 4.20){
           nucltransK(1.456,0.033,2.5e-4,0.5e-4);
           next915=true;
        }
        else if (p<= 5.85){
           nucltransK(0.724,0.033,1.3e-2,0.);
           next1647=true;
        }
        else if (p<=19.06){
           nucltransK(0.642,0.033,5.0e-3,0.);
           next1729=true;
        }
        else if (p<=37.37){
           nucltransK(0.574,0.033,9.0e-3,0.);
           next1797=true;
        }
        else if (p<=39.02){
           nucltransK(0.555,0.033,7.0e-3,0.);
           next1817=true;
        }
        else if (p<=53.13){
           nucltransK(0.479,0.033,3.5e-3,0.);
           next1893=true;
        }
        else if (p<=86.29){
           nucltransK(0.429,0.033,4.5e-3,0.);
           next1943=true;
        }
        else if (p<=97.10){
           nucltransK(0.397,0.033,5.2e-3,0.);
           next1975=true;
        }
        else if (p<=98.75){
           nucltransK(0.322,0.033,9.0e-3,0.);
           next2049=true;
        }
        else{
           nucltransK(0.110,0.033,1.6e-1,0.);
           next2262=true;
        }
     }
     if (next2262){
        p=100*GetRandom();
        if (p<=29.62){
           nucltransK(1.349,0.033,5.0e-4,0.3e-4);
           next913=true;
        }
        else if (p<=47.12){
           nucltransK(0.532,0.033,1.1e-2,0.);
           next1729=true;
        }
        else if (p<=52.69){
           nucltransK(0.464,0.033,3.5e-3,0.);
           next1797=true;
        }
        else if (p<=93.44){
           nucltransK(0.445,0.033,1.6e-2,0.);
           next1817=true;
        }
        else if (p<=95.63){
           nucltransK(0.369,0.033,6.0e-3,0.);
           next1893=true;
        }
        else{
           nucltransK(0.319,0.033,4.0e-2,0.);
           next1943=true;
        }
     }
     if (next2249){
        p=100*GetRandom();
        if (p<=34.30){
           nucltransK(0.472,0.033,1.3e-2,0.);
           next1777=true;
        }
        else if (p<=54.96){
           nucltransK(0.363,0.033,7.0e-3,0.);
           next1886=true;
        }
        else if (p<=81.82){
           nucltransK(0.355,0.033,7.0e-3,0.);
           next1893=true;
        }
        else{
           nucltransK(0.258,0.033,6.5e-2,0.);
           next1991=true;
        }
     }
     if (next2212){ 
        p=100*GetRandom();
        if (p<=7.42){
           nucltransK(1.299,0.033,4.0e-4,0.3e-4);
           next913=true;
        }
        else if (p<=60.70){
           nucltransK(0.435,0.033,1.7e-2,0.);
           next1777=true;
        }
        else if (p<=89.52){
           nucltransK(0.413,0.033,1.6e-2,0.);
           next1799=true;
        }
        else{
           nucltransK(0.221,0.033,1.0e-1,0.);
           next1991=true;
        }
     }
     if (next2049){
        p=100*GetRandom();
        if (p<=50.00){
           nucltransK(2.049,0.033,1.5e-4,2.2e-4);
           return;
        }
        else if (p<=61.48){
           nucltransK(1.137,0.033,1.7e-3,0.3e-4);
           next913=true;
        }
        else if (p<=75.00){
           nucltransK(1.135,0.033,1.7e-3,0.3e-4);
           next915=true;
        }
        else if (p<=90.98){
           nucltransK(0.743,0.033,3.5e-3,0.);
           next1307=true;
        }
        else{
           nucltransK(0.346,0.033,3.0e-2,0.);
           next1704=true;
        }
     }
     if (next2005){
        p=100*GetRandom();
        if (p<=72.63){
           nucltransK(2.005,0.033,8.0e-4,2.3e-4);
           return;
        }
        else if (p<=72.86){
           nucltransK(1.693,0.033,4.5e-4,0.7e-4);
           next312=true; 
        }
        else if (p<=75.21){
           nucltransK(1.091,0.033,2.0e-3,0.);
           next915=true;
        }
        else if (p<=95.30){
           nucltransK(0.698,0.033,5.5e-3,0.);
           next1307=true;
        }
        else if (p<=97.65){
           nucltransK(0.119,0.033,5.0e-1,0.);
           next1886=true;
        }
        else{
           nucltransK(0.112,0.033,5.5e-1,0.);
           next1893=true;
        } 
     } 
     if (next1991){ 
        p=100*GetRandom();
        if (p<= 1.96){
           nucltransK(1.078,0.033,6.0e-4,0.);
           next913=true;
        }
        else if (p<=91.68){
           nucltransK(0.262,0.033,6.5e-2,0.);
           next1729=true;
        }
        else if (p<=98.47){
           nucltransK(0.193,0.033,3.5e-2,0.);
           next1797=true;
        }
        else{
           nucltransK(0.098,0.033,2.2e-1,0.);
           next1893=true;
        }
     }
     if (next1975){ 
        p=100*GetRandom();
        if (p<= 1.13){
           nucltransK(1.975,0.033,4.0e-4,1.9e-4);
           return;
        } 
        else if (p<=49.41){
           nucltransK(1.062,0.033,2.0e-3,0.1e-4);
           next913=true;
        }
        else if (p<=51.15){
           nucltransK(1.060,0.033,2.0e-3,0.1e-4);
           next915=true;
        }
        else if (p<=59.29){
           nucltransK(0.520,0.033,8.0e-3,0.);
           next1455=true;
        }
        else if (p<=62.49){
           nucltransK(0.458,0.033,1.5e-2,0.);
           next1516=true;
        }
        else if (p<=67.43){
           nucltransK(0.198,0.033,3.2e-2,0.);
           next1777=true;
        }
        else if (p<=73.83){
           nucltransK(0.177,0.033,1.9e-1,0.);
           next1797=true; 
        }
        else{
           nucltransK(0.082,0.033,1.8,0.);
           next1893=true;
        }
     }
     if (next1943){
        p=100*GetRandom();
        if (p<=35.36){
           nucltransK(1.030,0.033,2.1e-3,0.);
           next913=true;
        }
        else if (p<=97.97){
           nucltransK(0.213,0.033,2.6e-2,0.);
           next1729=true; 
        }
        else{
           nucltransK(0.050,0.033,1.4e+1,0.);
           next1893=true;
        }
     }
     if (next1893){
        p=100*GetRandom();
        if (p<=2.13){
           nucltransK(1.893,0.033,4.0e-4,1.4e-4);
           return;
        } 
        else if (p<=23.29){
           nucltransK(0.980,0.033,2.5e-3,0.);
           next913=true; 
        }
        else if (p<=92.90){
           nucltransK(0.978,0.033,2.5e-3,0.);
           next915=true;
        }
        else if (p<=96.02){
           nucltransK(0.377,0.033,2.5e-2,0.);
           next1516=true;
        }
        else{
           nucltransK(0.116,0.033,1.5e-1,0.);
           next1777=true;
        }
     } 
     if (next1886){ 
        p=100*GetRandom();
        if (p<=18.33){
           nucltransK(1.886,0.033,2.0e-4,1.2e-4);
           return;
        }
        else if (p<=28.52){
           nucltransK(0.973,0.033,2.5e-3,0.);
           next913=true;
        }
        else if (p<=34.63){
           nucltransK(0.971,0.033,2.5e-3,0.);
           next915=true; 
        }
        else if (p<=36.67){
           nucltransK(0.369,0.033,2.6e-2,0.);
           next1516=true;
        }
        else if (p<=40.74){
           nucltransK(0.177,0.033,2.0e-1,0.);
           next1707=true;
        }
        else{
           nucltransK(0.088,0.033,1.4,0.);
           next1797=true;
        }
     }
     if (next1817){
        p=100*GetRandom();
        if (p<=52.17){
           nucltransK(0.040,0.033,3.2e+1,0.);
           next1777=true;
        }
        else{
           nucltransK(0.018,0.005,2.0e+1,0.);
           next1799=true;
        }
     }
     if (next1799){
        nucltransK(0.164,0.033,2.5e-1,0.);
        next1634=true;
     }
     if (next1797){
        p=100*GetRandom();
        if (p<=3.06){
           nucltransK(1.797,0.033,4.5e-4,0.9e-4);
           return;
        }
        else if (p<=19.90){
           nucltransK(0.885,0.033,3.0e-3,0.);
           next913=true;
        }
        else if (p<=57.49){
           nucltransK(0.883,0.033,3.0e-3,0.);
           next915=true;
        }
        else if (p<=65.99){
           nucltransK(0.343,0.033,3.0e-2,0.);
           next1455=true;
        }
        else if (p<=67.86){
           nucltransK(0.281,0.033,5.5e-2,0.);
        }
        else if (p<=75.34){
           nucltransK(0.151,0.033,3.3e-1,0.);
           next1647=true;
        }
        else{
           nucltransK(0.021,0.005,2.6,0.);
           next1777=true;
        }
     }
     if (next1777){ 
        p=100*GetRandom();
        if (p<=98.61){
           nucltransK(0.864,0.033,2.5e-3,0.);
           next913=true;
        }
        else{
           nucltransK(0.047,0.033,1.7e+1,0.);
           next1729=true;
        }
     }
     if (next1729){
        fThlev=170.e-9;
        p=100*GetRandom();
        if (p<=37.93){
           nucltransK(0.169,0.033,4.7e-2,0.);
           next1560=true;
        }
        else{
           nucltransK(0.095,0.033,2.1,0.);
           next1634=true;
        }
     }
     if (next1707){
        p=100*GetRandom();
        if (p<=90.91){
           nucltransK(0.795,0.033,4.0e-3,0.);
           next913=true;
        }
        else{
           nucltransK(0.793,0.033,3.0e-3,0.);
           next915=true;
        }
     }
     if (next1704){
        p=100*GetRandom();
        if (p<=52.17){
           nucltransK(1.704,0.033,6.0e-4,1.4e-4);
           return;
        }
        else if (p<=60.14){
           nucltransK(1.392,0.033,8.5e-4,0.2e-4);
           next312=true;
        }
        else if (p<=68.11){
           nucltransK(0.792,0.033,4.0e-3,0.);
           next913=true; 
        }
        else{
           nucltransK(0.790,0.033,4.0e-3,0.);
           next915=true;
        }
     }
     if (next1647){
        p=100*GetRandom();
        if (p<=9.88){
           nucltransK(1.646,0.033,5.0e-4,0.6e-4);
           return;
        }
        else if (p<=72.34){
           nucltransK(0.734,0.033,4.6e-3,0.);
           next913=true;
        }
        else if (p<=93.86){
           nucltransK(0.732,0.033,4.6e-3,0.);
           next915=true;
        }
        else{
           nucltransK(0.087,0.033,2.9,0.);
           next1560=true;
        }
     }
     if (next1634){
        fThlev=9.;
        nucltransK(0.074,0.033,2.4e+1,0.);
        next1560=true;
     }
     if (next1560){ 
        fThlev=0.;
        nucltransK(0.648,0.033,4.6e-3,0.);
        next913=true;
     }
     if (next1516){ 
        fThlev=0.;
        p=100*GetRandom();
        if (p<=78.54){
           nucltransK(1.516,0.033,8.0e-4,0.9e-4);
           return;
        }
        else if (p<=92.15){
           nucltransK(1.204,0.033,1.1e-3,0.1e-4);
           next312=true;
        }
        else{
           nucltransK(0.602,0.033,8.0e-3,0.);
           next915=true;
        }
     }
     if (next1455){
        fThlev=0.;
        p=100*GetRandom();
        if (p<=30.90){
           nucltransK(1.455,0.033,1.0e-3,0.9e-4);
           return;
        }
        else if (p<=87.98){
           nucltransK(1.143,0.033,1.8e-3,0.3e-4);
           next312=true;
        }
        else{ 
           nucltransK(0.540,0.033,1.0e-2,0.);
           next915=true;
        }
     }
     if (next1307){
        fThlev=0.;
        p=100*GetRandom();
        if (p<=39.00){
           nucltransK(1.307,0.033,1.2e-3,0.5e-4);
           return;
        }
        else if (p<=89.00){
           nucltransK(0.995,0.033,2.5e-3,0.);
           next312=true; 
        }
        else{
           nucltransK(0.392,0.033,2.1e-2,0.);
           next915=true; 
        }
     }
     if (next915){
        fThlev=0.;
        p=100*GetRandom();
        if (p<=99.84){
           nucltransK(0.915,0.033,2.7e-3,0.);
           return;
        }
        else{
           nucltransK(0.602,0.033,5.0e-3,0.);
           next312=true;
        }
     }
     if (next913){
        fThlev=0.;
        nucltransK(0.913,0.033,2.1e-3,0.);
        return;
     }
     if (next312){
        fThlev=0.;
        nucltransK(0.312,0.033,4.0e-2,0.);
        return;
     }
  }
}
//-----------------------------------------------

void Decay::Te134(float tcnuc)
{
  // Model for scheme of Te134 decay(Yu.V.Sergeenkov,Nucl.Data Sheets 71(1994)557).
  fThnuc=2508.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next923=false, next847=false, next645=false, next210=false;
  bool next181=false, next79=false, next44=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<=14.){
     beta(0.454,0.,0.);
     p=100*GetRandom();
     if (p<= 3.33){
        nucltransK(1.027,0.033,1.5e-3,0.);
        next79=true;
     }
     else if (p<=14.45){
        nucltransK(0.926,0.033,2.5e-3,0.);
        next181=true;
     }
     else if (p<=17.78){
        nucltransK(0.896,0.033,3.0e-3,0.);
        next210=true;
     }
     else if (p<=91.15){
        nucltransK(0.461,0.033,1.5e-2,0.);
        next645=true;
     }
     else if (p<=94.70){
        nucltransK(0.260,0.033,6.0e-2,0.);
        next847=true; 
     }
     else{
        nucltransK(0.183,0.033,1.8e-1,0.);
        next923=true;
     }
  }
  else if (pbeta<=58.){
     beta(0.637,0.,0.);
     next923=true;
  }
  else{
     beta(0.713,0.,0.);
     next847=true;
  }

  if (next923){
     p=100*GetRandom();
     if (p<= 2.71){
        nucltransK(0.844,0.033,2.5e-3,0.);
        next79=true;
     }
     else if (p<=37.31){
        nucltransK(0.743,0.033,4.5e-3,0.);
        next181=true;
     }
     else if (p<=47.94){
        nucltransK(0.713,0.033,5.0e-3,0.);
        next210=true;
     }
     else if (p<=98.37){
        nucltransK(0.278,0.033,4.9e-2,0.);
        next645=true;
     }
     else{
        nucltransK(0.077,0.033,1.61,0.);
        next847=true;
     }
  }
  if (next847){
     p=100*GetRandom();
     if (p<=69.71){
        nucltransK(0.767,0.033,3.3e-3,0.);
        next79=true;
     }
     else if (p<=72.49){
        nucltransK(0.666,0.033,6.0e-3,0.);
        next181=true;
     }
     else if (p<=76.45){
        nucltransK(0.636,0.033,7.0e-3,0.);
        next210=true; 
     }
     else{
        nucltransK(0.201,0.033,1.3e-1,0.);
        next645=true; 
     }
  } 
  if (next645){
     p=100*GetRandom();
     if (p<= 2.03){
        nucltransK(0.645,0.033,5.0e-3,0.);
        return;
     }
     else if (p<=45.01){
        nucltransK(0.566,0.033,9.0e-3,0.);
        next79=true;
     }
     else if (p<=55.87){
        nucltransK(0.465,0.033,1.5e-2,0.);
        next181=true;
     }
     else{
        nucltransK(0.435,0.033,1.4e-2,0.);
        next210=true;
     }
  }
  if (next210){
     fThlev=0.15e-9;
     p=100*GetRandom();
     if (p<=98.94){
        nucltransK(0.210,0.033,1.1e-1,0.);
        return;
     }
     else{
        nucltransK(0.131,0.033,5.2e-1,0.);
        next79=true;
     }
  } 
  if (next181){
     fThlev=0.1e-9;
     p=100*GetRandom();
     if (p<=95.83){
        nucltransK(0.181,0.033,1.8e-1,0.);
        return;
     }
     else if (p<=96.45){
        nucltransK(0.137,0.033,5.8e-1,0.);
        next44=true;
     }
     else{
        nucltransK(0.101,0.033,1.2,0.);
        next79=true;
     }
  }
  if (next79){
     fThlev=1.62e-9;
     nucltransK(0.079,0.033,1.50,0.);
     return; 
  } 
  if (next44){
     fThlev=0.;
     nucltransK(0.044,0.033,7.97,0.);
     return; 
  }
}
//-----------------------------------------------

void Decay::Th234(float tcnuc)
{ 
  // Scheme of Th234 decay to Pa234m (not Pa234!) in accordance with NDS 
  // 108(2007)681 and ENSDF database at NNDC site on 8.08.2007.
  fThnuc=2.08224e6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next167=false, next103=false, next74=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<= 1.501){
     beta(0.086,0.,0.);
     p=100*GetRandom();
     if (p<=17.0){
        nucltransK(0.113,0.021,2.3e-1,0.);
        next74=true;
     }
     else if (p<=21.7){
        nucltransK(0.083,0.021,2.0e-1,0.);
        next103=true;
     }
     else{
        nucltransK(0.020,0.005,2.4e2,0.);
        next167=true;
     }
  }
  else if (pbeta<= 1.516){
     beta(0.096,0.,0.); 
     nucltransK(0.103,0.021,3.8e0,0.);
     next74=true;
  }
  else if (pbeta<= 7.921){
     beta(0.106,0.,0.);
     next167=true; 
  }
  else if (pbeta<=21.933){
     beta(0.106,0.,0.);
     fThlev=0.;
     p=100*GetRandom();
     if (p<=97.0){
        nucltransK(0.092,0.021,5.3e0,0.);
        next74=true;
     }
     else{
        nucltransK(0.062,0.021,2.5e1,0.);
        next103=true;
     }
  }
  else{
     beta(0.199,0.,0.);
     next74=true;
  }
  
  if (next167){
     fThlev=0.55e-9;
     p=100*GetRandom();
     if (p<=32.1){
        nucltransK(0.093,0.021,1.5e-1,0.);
        next74=true;
     }
     else{
        nucltransK(0.063,0.021,4.1e-1,0.);
        next103=true;
     }
  }
  if (next103){
     fThlev=0.;
     nucltransK(0.029,0.021,4.4e3,0.);
     next74=true; 
  }
  if (next74){
     /// below is creation of Pa234m with T1/2=1.159 m which mainly 
     /// beta- decays to U234 (IT to Pa234 is only 0.16%); 
     /// decay of Pa234m should be generated independently
     return;
  } 

}
//-----------------------------------------------

void Decay::Tl208(float tcnuc)
{
  // Scheme of Tl208 decay ("Table of Isotopes", 7th ed., 1978).
  // VIT, 27.07.1992, 22.10.1995.
  // VIT, 11.05.2005, updated to NDS 47(1986)797 (last NDS issue for A=208).
  // VIT, 4.02.2009, updated to NDS 108(2007)1583; also were added:
  // LM conversion electrons; more complex emission of X rays 
  // emitted in K conversion.
  fThnuc=183.18;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next3708=false, next3475=false, next3198=false, next2615=false;
  float pbeta, p;
  pbeta=100.*GetRandom(); 
  if (pbeta<=0.052){
     beta(0.518,0.,0.); 
     nucltransKLM_Pb(1.283,0.088,7.75e-3,0.015,1.27e-3,0.004,0.41e-3,2.3e-5);
     next3198=true;
  }
  else if (pbeta<= 0.069){
     beta(0.616,0.,0.);
     nucltransKLM_Pb(1.185,0.088,9.49e-3,0.015,1.56e-3,0.004,0.47e-3,4.9e-6);
     next3198=true;
  }
  else if (pbeta<= 0.113){
     beta(0.641,0.,0.);
     p=100*GetRandom();
     if (p<= 4.55){
        nucltransKLM_Pb(1.744,0.088,3.56e-3,0.015,0.58e-3,0.004,0.17e-3,2.6e-4);
        next2615=true;
     }
     else if (p<=29.55){
        nucltransKLM_Pb(1.161,0.088,9.99e-3,0.015,1.64e-3,0.004,0.51e-3,2.6e-6);
        next3198=true;
     }
     else{
        nucltransKLM_Pb(0.883,0.088,20.13e-3,0.015,3.33e-3,0.004,1.02e-3,0.);
        next3475=true;
     }
  }
  else if (pbeta<= 0.118){
     beta(0.676,0.,0.);
     nucltransKLM_Pb(1.126,0.088,1.69e-3,0.015,0.26e-3,0.004,0.08e-3,2.1e-6);
     next3198=true;
  }
  else if (pbeta<= 0.219){
     beta(0.703,0.,0.); 
     p=100*GetRandom();
     if (p<=87.23){
        nucltransKLM_Pb(0.821,0.088,2.43e-2,0.015,0.40e-2,0.004,0.12e-2,0.);
        next3475=true;
     }
     else{
        nucltransKLM_Pb(0.588,0.088,5.78e-2,0.015,0.97e-2,0.004,0.29e-2,0.);
        next3708=true;
     }
  }
  else if (pbeta<= 0.221){
     beta(0.737,0.,0.);
     nucltransKLM_Pb(1.648,0.088,4.11e-3,0.015,0.67e-3,0.004,0.20e-3,0.19e-3);
     next2615=true;
  }
  else if (pbeta<= 0.448){
     beta(0.819,0.,0.);
     p=100*GetRandom();
     if (p<=90.31){
        nucltransKLM_Pb(0.983,0.088,1.53e-2,0.015,0.25e-2,0.004,0.08e-2,0.);
        next3198=true;
     }
     else{
        nucltransKLM_Pb(0.705,0.088,3.60e-2,0.015,0.60e-2,0.004,0.18e-2,0.); 
        next3475=true;
     }
  }
  else if (pbeta<= 0.623){
     beta(0.874,0.,0.); 
     p=100*GetRandom();
     if (p<=96.15){
        nucltransKLM_Pb(0.928,0.088,1.77e-2,0.015,0.29e-2,0.004,0.10e-2,0.);
        next3198=true; 
     }
     else{
        nucltransKLM_Pb(0.650,0.088,4.45e-2,0.015,0.75e-2,0.004,0.22e-2,0.);
        next3475=true;
     }
  }
  else if (pbeta<= 0.630){
     beta(1.003,0.,0.);
     nucltransKLM_Pb(1.381,0.088,6.43e-3,0.015,1.05e-3,0.004,0.32e-3,0.05e-3);
     next2615=true;
  }
  else if (pbeta<= 3.810){
     beta(1.038,0.,0.);
     p=100*GetRandom();
     if (p<=51.25){
        nucltransKLM_Pb(0.763,0.088,2.93e-2,0.015,0.49e-2,0.004,0.14e-2,0.);
        next3198=true;
     }
     else if (p<=64.82){
        nucltransKLM_Pb(0.486,0.088,9.54e-2,0.015,1.61e-2,0.004,0.49e-2,0.);
        next3475=true;
     } 
     else{
        nucltransKLM_Pb(0.253,0.088,51.60e-2,0.015,8.83e-2,0.004,2.57e-2,0.);
        next3708=true;
     }
  }
  else if (pbeta<= 3.856){
     beta(1.053,0.,0.);
     nucltransKLM_Pb(0.749,0.088,3.08e-2,0.015,0.51e-2,0.004,0.16e-2,0.);
     next3198=true;
  }
  else if (pbeta<= 4.486){
     beta(1.079,0.,0.);
     p=100*GetRandom();
     if (p<=39.49){
        nucltransKLM_Pb(0.722,0.088,3.20e-2,0.015,0.54e-2,0.004,0.16e-2,0.);
        next3198=true;
     }
     else{
        nucltransKLM_Pb(0.211,0.088,9.22e-1,0.015,1.59e-1,0.004,0.45e-1,0.);
        next3708=true;
     }
  }
  else if (pbeta<=28.686){
     beta(1.291,0.,0.);
     next3708=true;
  }
  else if (pbeta<=50.886){
     beta(1.524,0.,0.);
     next3475=true;
  }
  else{
     beta(1.801,0.,0.);
     next3198=true;
  }

  if (next3708){
     p=100*GetRandom();
     if (p<=1.66){
        nucltransKLM_Pb(1.094,0.088,4.49e-3,0.015,0.84e-3,0.004,0.27e-3,0.);
        next2615=true;
     }
     else if (p<=97.95){
        nucltransKLM_Pb(0.511,0.088,8.42e-2,0.015,1.42e-2,0.004,0.43e-2,0.);
        next3198=true;
     }
     else{
        nucltransKLM_Pb(0.233,0.088,5.47e-1,0.015,1.16e-1,0.004,0.37e-1,0.); 
        next3475=true;
     }
  }
  if (next3475){ 
     fThlev=4.e-12;
     p=100*GetRandom();
     if (p<=55.95){
        nucltransKLM_Pb(0.861,0.088,2.17e-2,0.015,0.36e-2,0.004,0.11e-2,0.);
        next2615=true;
     }
     else{
        nucltransKLM_Pb(0.277,0.088,4.36e-1,0.015,0.75e-1,0.004,0.22e-1,0.);
        next3198=true;
     } 
  }
  if (next3198){ 
     fThlev=294.e-12;
     nucltransKLM_Pb(0.583,0.088,1.51e-2,0.015,0.41e-2,0.004,0.13e-2,0.);
     next2615=true;
  }
  if (next2615){
     fThlev=16.7e-12;
     nucltransKLM_Pb(2.615,0.088,1.71e-3,0.015,0.29e-3,0.004,0.10e-3,0.37e-3);
  }

}  
//-----------------------------------------------

void Decay::Xe133(float tcnuc)
{
  // Scheme of Xe133 decay ("Table of Isotopes", 8th ed., 1996).
  // VIT, 18.08.1997. Updated 5.12.2003 in accordance with NDS 75(1995)491.
  fThnuc=452995.2;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next161=false, next81=false;
  float pbeta, p;
  pbeta=100.*GetRandom(); 
  if (pbeta<=0.008){
     beta(0.044,0.,0.);
     fThlev=21.e-12;
     p=100.*GetRandom();
     if (p<=32.313){
        nucltransK(0.384,0.036,2.0e-2,0.);
        return; 
     } 
     else if (p<=98.259){
        nucltransK(0.303,0.036,4.4e-2,0.);
        next81=true;  
     }
     else{
        nucltransK(0.223,0.036,9.8e-2,0.); 
        next161=true;
     }
  }
  else if (pbeta<=0.818){
     beta(0.267,0.,0.);
     next161=true;
  }
  else{
     beta(0.346,0.,0.);
     next81=true;
  }

  if (next161){
     fThlev=172.e-12;
     p=100.*GetRandom();
     if (p<=10.287){
        nucltransK(0.161,0.036,3.0e-1,0.);
        return;
     }
     else{
        nucltransK(0.080,0.036,1.8e+0,0.);
        next81=true;
     }
  }
  if (next81){
     fThlev=6.28e-9;
     nucltransK(0.081,0.036,1.7e+0,0.);
     return;
  }
} 
//-----------------------------------------------

void Decay::Xe135(float tcnuc)
{
  // Model of Xe135 decay(Yu.V.Sergeenkov et al.,Nucl.Data Sheets 84(1998)115).
  // VIT, 9.10.2002.
  fThnuc=32904.;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next608=false, next408=false, next250=false;
  float pbeta, p;
  pbeta=100.*GetRandom(); 
  if (pbeta<=0.123){
     beta(0.089,0.,0.);
     p=100.*GetRandom();
     if (p<= 3.34){
        nucltransK(1.062,0.036,2.0e-3,0.);
        return; 
     }
     else if (p<=60.40){
        nucltransK(0.813,0.036,3.8e-3,0.);
        next250=true;
     }
     else if (p<=97.07){
        nucltransK(0.654,0.036,6.5e-3,0.);
        next408=true;
     }
     else{
        nucltransK(0.454,0.036,1.5e-2,0.);
        next608=true; 
     }
  }
  else if (pbeta<=0.198){
     beta(0.170,0.,0.);
     p=100.*GetRandom();
     if (p<=73.53){
        nucltransK(0.732,0.036,5.0e-3,0.);
        next250=true;
     }
     else if (p<=79.95){
        nucltransK(0.573,0.036,9.0e-3,0.);
        next408=true;
     }
     else{
        nucltransK(0.373,0.036,2.5e-2,0.);
        next608=true;
     }
  }
  else if (pbeta<=3.311){
     beta(0.543,0.,0.);  
     next608=true;
  }
  else if (pbeta<=3.902){
     beta(0.743,0.,0.);
     next408=true;
  }
  else{
     beta(0.901,0.,0.);
     next250=true;
  }

  if (next608){
     p=100.*GetRandom();
     if (p<= 92.42){
        nucltransK(0.608,0.036,7.5e-3,0.);
        return;
     } 
     else if (p<= 99.62){
        nucltransK(0.358,0.036,2.7e-2,0.);
        next250=true;
     }
     else{
        nucltransK(0.200,0.036,1.4e-1,0.);
        next408=true;
     }
  }
  if (next408){
     p=100.*GetRandom();
     if (p<=55.33){
        nucltransK(0.408,0.036,2.0e-2,0.);
        return;
     }
     else{
        nucltransK(0.158,0.036,2.5e-1,0.);
        next250=true;
     }
  }
  if (next250){
     fThlev=0.28e-9;
     nucltransK(0.250,0.036,7.6e-2,0.);
     return;
  }
}
//-----------------------------------------------

void Decay::Y88(float tcnuc)
{
  // Scheme of Y88 decay ("Table of Isotopes", 7th ed., 1978).
  // Accuracy in description of: decay branches       - 0.001%,
  //                           : deexcitation process - 0.001%.
  //  VIT, 20.07.1993, 22.10.1995.
  //  VIT, 12.11.2006 (update to NDS 105(2005)419 and change 
  //  of beta+ spectrum shape)
  fThnuc=9.2124864e+6;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  bool next2734=false, next1836=false;
  float pdecay, p;
  
  pdecay=100.*GetRandom();
  if (pdecay<=0.065){
     particle(1,0.016,0.016,0.,pi,0.,twopi,0.,0.);
     fThlev=0.14e-9;
     nucltransK(0.851,0.016,8.5e-4,0.); 
     next2734=true; 
  }
  else if (pdecay<=0.093){
     particle(1,0.016,0.016,0.,pi,0.,twopi,0.,0.);
     fThlev=0.13e-12;
     p=100.*GetRandom();
     if (p<=25.){
        nucltransK(3.219,0.016,6.0e-5,8.7e-4);
        return;
     }
     else{
        nucltransK(1.382,0.016,2.6e-4,4.8e-5);
        next1836=true;
     }
  }
  else if (pdecay<=94.490){
     particle(1,0.016,0.016,0.,pi,0.,twopi,0.,0.);
     next2734=true;
  }  
  else if (pdecay<=99.79){
     particle(1,0.016,0.016,0.,pi,0.,twopi,0.,0.);
     next1836=true;
  }
  else{  //b+ to Sr88
     beta2f(0.765,0.,0.,1,1.,0.,0.,0.);
     next1836=true;
  }

  if (next2734){
     fThlev=0.78e-12;
     p=100.*GetRandom();
     if (p<=0.75){
        nucltransK(2.734,0.016,1.2e-4,3.3e-4);
        return;
     }
     else{
        nucltransK(0.898,0.016,3.1e-4,0.);
        next1836=true;
     }
  }
  if (next1836){
     fThlev=0.148e-12;
     nucltransK(1.836,0.016,1.4e-4,2.3e-4);
     return;
  }
}
//-----------------------------------------------

void Decay::Zn65(float tcnuc)
{
  // Scheme of Zn65 decay (NDS 69(1993)209 and NNDC online corrections 
  // on 28.03.2007).
  fThnuc=2.1086784e+07;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  float pec, p;
  pec=100.*GetRandom(); 
  if (pec<=1.42){
     pair(0.329); // beta+ decay to g.s. of 65-Cu
     return;
  }
  else{  // X ray after EC to 65-Cu
     particle(1,0.009,0.009,0.,pi,0.,twopi,0,0);
     if (pec<=52.03){
        fThlev=0.285e-12;
        p=100.*GetRandom();
        if (p<=99.994){
           nucltransK(1.116,0.009,1.9e-4,1.0e-6);
           return;
        }
        else{
           nucltransK(0.345,0.009,6.8e-3,0.);
           fThlev=99.e-15;
           nucltransK(0.771,0.009,3.8e-4,0.);
           return;
        }
     }
     else return;  
  }  
}
//-----------------------------------------------

void Decay::Nb96(float tcnuc)
{
  // Scheme of Nb96 decay (NDS 68(1993)165 and ENSDF at 
  // NNDC site on 13.11.2007).
  // VIT, 7.05.1998; update 13.11.2007.
  fThnuc=8.406e4;
  fTdnuc=tcnuc-fThnuc/log(2.)*log(GetRandom());
  fTclev=0.;
  fThlev=0.;
  fZdtr=42;
  bool next2755=false, next2441=false, next2439=false, next2219=false;
  bool next1978=false, next1870=false, next1628=false, next1626=false;
  bool next1498=false, next778=false;
  float pbeta, p;
  pbeta=100.*GetRandom();
  if (pbeta<=0.024){
     beta(0.212,0.,0.);
     nucltransK(1.347,0.020,4.3e-4,0.3e-4);
     next1628=true;
  }
  else if (pbeta<=0.524){
     beta(0.311,0.,0.);
     p=100.*GetRandom();
     if (p<=93.26){
        nucltransK(0.435,0.020,5.4e-3,0.);
        next2441=true; 
     }
     else{
        nucltransK(0.120,0.020,1.5e-1,0.);
        next2755=true;
     }
  }
  else if (pbeta<=1.014){
     beta(0.432,0.,0.);
     next2755=true;
  }
  else if (pbeta<=3.314){
     beta(0.746,0.,0.);
     next2441=true;
  }
  else{
     beta(0.748,0.,0.);
     next2439=true;
  }
 
  if (next2755){
     p=100.*GetRandom();
     if (p<=75.63){
        nucltransK(1.127,0.020,5.2e-4,0.1e-5);
        next1628=true;
     }
     else if (p<=86.64){
        nucltransK(0.316,0.020,1.2e-2,0.);
        next2439=true;
     }
     else{
        nucltransKLM(0.314,0.020,1.1e-2,0.003,1.3e-3,0.001,4.2e-4,0.);
        next2441=true;
     }
  }
  if (next2441){  
     nucltransK(0.813,0.020,1.3e-3,0.);
     next1628=true;
  }
  if (next2439){
     p=100.*GetRandom();
     if (p<=11.20){
        nucltransKLM(0.811,0.020,1.1e-3,0.003,1.3e-4,0.001,2.2e-5,0.);
        next1628=true;
     }
     else if (p<=69.87){
        nucltransKLM(0.569,0.020,2.6e-3,0.003,2.8e-4,0.001,5.8e-5,0.);
        next1870=true;
     }
     else if (p<=96.89){
        nucltransKLM(0.460,0.020,5.3e-3,0.003,6.4e-4,0.001,1.3e-4,0.);
        next1978=true;
     }
     else{
        nucltransKLM(0.219,0.020,3.2e-2,0.003,4.0e-3,0.001,1.4e-3,0.);
        next2219=true;
     }
  }
  if (next2219){
     p=100.*GetRandom();
     if (p<=11.19){
        nucltransK(1.441,0.020,4.3e-4,0.5e-4);
        next778=true;
     }
     else if (p<=36.94){
        nucltransKLM(0.722,0.020,1.5e-3,0.003,1.7e-4,0.001,3.6e-5,0.);
        next1498=true;
     }
     else if (p<=44.75){
        nucltransKLM(0.593,0.020,2.3e-3,0.003,2.6e-4,0.001,5.2e-5,0.);
        next1626=true;
     }
     else if (p<=68.44){
        nucltransKLM(0.591,0.020,2.4e-3,0.003,2.8e-4,0.001,4.9e-5,0.);
        next1628=true;
     }
     else if (p<=80.66){
        nucltransKLM(0.350,0.020,1.0e-2,0.003,1.2e-3,0.001,5.0e-4,0.);
        next1870=true;
     }
     else{
        nucltransKLM(0.241,0.020,2.1e-2,0.003,2.4e-3,0.001,4.9e-4,0.);
        next1978=true;
     }
  }
  if (next1978){
     p=100.*GetRandom();
     if (p<=71.90){
        nucltransK(1.200,0.020,4.6e-4,7.7e-6);
        next778=true;
     } 
     else if (p<=92.93){
        nucltransK(0.481,0.020,4.5e-3,0.);
        next1498=true;
     } 
     else if (p<=95.95){
        nucltransKLM(0.353,0.020,1.0e-2,0.003,1.2e-3,0.001,7.0e-4,0.);
        next1626=true;
     }
     else if (p<=99.81){
        nucltransKLM(0.350,0.020,1.0e-2,0.003,1.2e-3,0.001,5.0e-4,0.);
        next1628=true;
     }
     else{
        nucltransKLM(0.109,0.020,1.7e-1,0.003,2.0e-2,0.001,4.2e-3,0.);
        next1870=true;
     }
  }
  if (next1870){
     fThlev=6.4e-12;
     p=100.*GetRandom();
     if (p<=88.61){
        nucltransK(1.091,0.020,5.6e-4,0.);
        next778=true; 
     }
     else if (p<=93.45){
        nucltransKLM(0.372,0.020,1.0e-2,0.003,1.3e-3,0.001,4.3e-4,0.);
        next1498=true;
     }
     else{
        nucltransKLM(0.241,0.020,2.1e-2,0.003,2.4e-3,0.001,7.8e-4,0.);
        next1628=true;
     }
  } 
  if (next1628){
     fThlev=1.2e-12;
     nucltransKLM(0.850,0.020,1.0e-3,0.003,1.1e-4,0.001,2.3e-5,0.);
     next778=true;
  }
  if (next1626){
     p=100.*GetRandom();
     if (p<=11.72){
        nucltransK(1.626,0.020,2.8e-4,1.3e-4);
        return;
     }
     else if (p<=98.63){
        nucltransKLM(0.848,0.020,1.0e-3,0.003,1.1e-4,0.001,2.3e-5,0.);
        next778=true;
     }
     else{
        nucltransKLM(0.128,0.020,1.1e-1,0.003,1.3e-2,0.001,2.7e-3,0.);
        next1498=true;
     }
  }
  if (next1498){
     fThlev=0.78e-12;
     p=100.*GetRandom();
     if (p<=32.43){
        nucltransK(1.498,0.020,3.3e-4,0.8e-4);
        return;
     } 
     else{
        nucltransKLM(0.720,0.020,1.5e-3,0.003,1.7e-4,0.001,3.5e-5,0.);
        next778=true;
     }
  }
  if (next778){
     fThlev=3.67e-12;
     nucltransKLM(0.778,0.020,1.2e-3,0.003,1.4e-3,0.001,2.9e-5,0.);
     return;
  }
}
//-----------------------------------------------
//-----------------------------------------------
void Decay::particle(int np,float E1,float E2,float teta1,float teta2,float phi1,float phi2,float tclev,float thlev)
{
  //Generation of isotropical emission of particle in the range of  energies and angles.
  //  E1,E2       - range of kinetic energy of particle (MeV);
  //  teta1,teta2 - range of teta angle (radians);
  //  phi1,phi2   - range of phi  angle (radians);
  //  fTclev       - time of creation of level from which particle will be emitted (sec);
  //  fThlev       - level halflife (sec).
  //  tdlev       - time of decay of level (sec);
  //  tevst       - time of event's start (sec);
  fTclev=tclev;
  
  fNbPart=fNbPart+1;
  if (fNbPart>100) printf("PARTICLE: in event more than 100 particles: %d",fNbPart); 
 
  if (np<1 || np>50 || (np>32 && np<45)) printf("PARTICLE: unknown particle number: %d",np);
  
  npgeant[fNbPart]=np;
  pmass=datamass[np];
  phi=phi1+(phi2-phi1)*GetRandom();

  float ctet1=1.;
  float ctet2=-1.;
  float E, p, ctet;
  
  if(teta1!=0.) ctet1=cos(teta1);
  if(teta2!=pi) ctet2=cos(teta2);

  ctet=ctet1+(ctet2-ctet1)*GetRandom();
  float stet=sqrt(1.-ctet*ctet);
  E=E1;
  if(E1!=E2) E=E1+(E2-E1)*GetRandom();

  p=sqrt(E*(E+2.*pmass));
  fPmoment[0][fNbPart]=p*stet*cos(phi);
  fPmoment[1][fNbPart]=p*stet*sin(phi);
  fPmoment[2][fNbPart]=p*ctet;
  fTdlev=fTclev;

  if(thlev>0) fTdlev=fTclev-thlev/log(2.)*log(GetRandom());
  fPtime[fNbPart]=fTdlev;
  return; 
}

//---------------------------------------
void Decay::PbAtShell(int KLMenergy){
  // Subroutine describes in some approximation the deexcitation 
  // process in atomic shell of Pb after creation of electron 
  // vacation in K, L or M shell.
  // It is supposed electron binding energy on 
  // Pb K-shell = 88 keV, on L-shell = 15 keV, 
  // on M-shell = 3 keV. 
  // The following values of KLMenergy are allowed:
  // 88 (hole in K-shell), 15 (in L-shell) and 3 (in M-shell).
  // VIT, 7.07.1995, 22.10.1995.
  float p;
  int Lhole=0, Mhole=0;
  bool gonext=false;
  if (KLMenergy==88){
     p=100*GetRandom();
     if (p<=22.){
        particle(1,0.085,0.085,0.,pi,0.,twopi,0.,0.);
        Mhole=Mhole+1;
        for (int i=1; i<=Mhole;i++){
            particle(1,0.003,0.003,0.,pi,0.,twopi,0.,0.);
        }
      return;
     }
     else{
        p=100*GetRandom();
        if (p<=96.) particle(1,0.073,0.073,0.,pi,0.,twopi,0.,0.);
        else{
           particle(3,0.058,0.058,0.,pi,0.,twopi,0.,0.);
           Lhole=Lhole+1;
        }
     }
     Lhole=Lhole+1; 
     gonext=true; 
  }
  if (KLMenergy==15){
     Lhole=1;
     gonext=true; 
  }   
  if (gonext){
     for (int i=1; i<=Lhole;i++){
         p=100*GetRandom();
         if (p<=40.) particle(1,0.012,0.012,0.,pi,0.,twopi,0.,0.);
         else{
            particle(3,0.009,0.009,0.,pi,0.,twopi,0.,0.);
            Mhole=Mhole+1;
         }
         Mhole=Mhole+1;
     }  
     gonext=true;
  }
  if (KLMenergy==3){
     Mhole=1;
     gonext=true;
  } 
  if (gonext){
     for (int i=1; i<=Mhole;i++){
         particle(1,0.003,0.003,0.,pi,0.,twopi,0.,0.);
     }
     return;
  }
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////
float Decay::fe1_mod(float &e1)
{
  //probability distribution for energy of first e-/e+ for fModebb=1,2,3,7,10
  float fe1mod=0.;
  if(e1>e0) return fe1mod;

  float e2=e0-e1;
  float p1=sqrt(e1*(e1+2.*emass));
  float p2=sqrt(e2*(e2+2.*emass));
  if(fModebb==1) fe1mod=(e1+emass)*p1*fermi(fZdbb,e1)*(e2+emass)*p2*fermi(fZdbb,e2);
  if(fModebb==2) fe1mod=(e1+emass)*p1*fermi(fZdbb,e1)*(e2+emass)*p2*fermi(fZdbb,e2)
                        * (e0-2.*e1)*(e0-2.*e1);
  if(fModebb==3) fe1mod=p1*fermi(fZdbb,e1)*p2*fermi(fZdbb,e2)*(2.*p1*p1*p2*p2
                        +9.*((e1+emass)*(e2+emass)+emass*emass)*(p1*p1+p2*p2));

  if(fModebb==7) fe1mod=p1*fermi(fZdbb,e1)*p2*fermi(fZdbb,e2)*
                       ((e1+emass)*(e2+emass)+emass*emass)*(p1*p1+p2*p2);

  if(fModebb==10) fe1mod=(e1+emass)*p1*fermi(fZdbb,e1)*pow(e0-e1,5);

  return fe1mod;
}

/////////////////////////////////////////////////////////////////////////////////////////////
float Decay::fe2_mod(float &e2)
{
  //probability distribution for energy of second e-/e+ for fModebb=4,
  float fe2mod=0.;
  if(e2>e0-e1) return fe2mod;

  float p2=sqrt(e2*(e2+2.*emass));
  if(fModebb==4)  fe2mod=(e2+emass)*p2*fermi(fZdbb,e2)*pow(e0-e1-e2,5);
  if(fModebb==5)  fe2mod=(e2+emass)*p2*fermi(fZdbb,e2)*(e0-e1-e2);
  if(fModebb==6)  fe2mod=(e2+emass)*p2*fermi(fZdbb,e2)*pow(e0-e1-e2,3);
  if(fModebb==8)  fe2mod=(e2+emass)*p2*fermi(fZdbb,e2)*pow(e0-e1-e2,7)*pow(e1-e2,2);
  if(fModebb==13) fe2mod=(e2+emass)*p2*fermi(fZdbb,e2)*pow(e0-e1-e2,7);
  if(fModebb==14) fe2mod=(e2+emass)*p2*fermi(fZdbb,e2)*pow(e0-e1-e2,2);
  if(fModebb==15) fe2mod=(e2+emass)*p2*fermi(fZdbb,e2)*pow(e0-e1-e2,5)*
                        (9*pow(e0-e1-e2,2)+21*pow(e2-e1,2));
  if(fModebb==16) fe2mod=(e2+emass)*p2*fermi(fZdbb,e2)*pow(e0-e1-e2,5)*pow(e2-e1,2);

  return fe2mod;
}

/////////////////////////////////////////////////////////////////////////////////////////////
float Decay::fermi(const float &Z,const float &E)
{
 //Function fermi calculates the traditional function of Fermi 
 //in theory of beta decay to take into account the Coulomb correction
 //to the shape of electron/positron energy spectrum.
 // Input : Z - atomic number of daughter nuclei (>0 for e-, <0 for e+);
 //         E - kinetic energy of particle (MeV; E>50 eV).
 // Output: corr - value of correction factor (without normalization -
 //                constant factors are removed - for MC simulation).
  float dE=E;
  if (E<50e-06) dE=50e-06;

  float alfaz=Z/137.036;
  float w=dE/0.511+1.;
  double  p=sqrt(w*w-1.);
  float   y=alfaz*w/p;
  float   g=sqrt(1.-alfaz*alfaz);

  complex<double> carg;
  
  carg=complex<double>(g,y);
  float fermif=pow(p,(2.*g-2.))*exp(pi*y+2.*log(abs(cgamma(carg))));
  return fermif;
}

/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::beta(float Qbeta,float tcnuc,float thnuc)
{
  // Subroutine beta simulates the angles and energy of beta particles emitted 
  // in the beta decay of nucleus. The decay is considered as allowed.
  // Only Coulomb correction to the shape of energy spectrum is taken into consideration.
  // tcnuc - time of creation of nucleus (sec);
  // thnuc - nucleus halflife (sec).
  // tevst - time of event's start (sec);
  fQbeta=Qbeta;
  TF1 fbeta ("funbeta",this,&Decay::funbeta,0,fQbeta,3);
  float em=0.,fm=0., E,fe=0.;
  float f;
  int np;
  tgold(50.e-6,fQbeta,fbeta,0.001*fQbeta,2,em,fm);
  do{
     E=50.e-6+(fQbeta-50.e-6)*GetRandom();
     fe=fbeta.Eval(E);
     f=fm*GetRandom();
  }while(f>fe);
  if (fZdtr>=0.) np=3;
  if (fZdtr<0.)  np=2;
  particle(np,E,E,0.,pi,0.,twopi,tcnuc,thnuc);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::beta1f(float Qbeta,float tcnuc,float thnuc,float c1,float c2,float c3,float c4)
{
  // Subroutine beta1 simulates the angles and energy of 
  // beta particles emitted in beta decay of nucleus. 
  // The decay is considered as forbidden;
  // correction factor to the allowed spectrum shape has a form
  // typical for empirical corrections: 
  // cf(e)=(1+c1/w+c2*w+c3*w**2+c4*w**3), w=e/emass+1.
  // VIT, 30.07.1992; 15.10.1995; 31.03.2006

  fQbeta=Qbeta;
  fC1=c1;
  fC2=c2;
  fC3=c3;
  fC4=c4;
  TF1 fbeta1f ("funbeta1f",this,&Decay::funbeta1f,0,fQbeta,3);
  float em=0.,fm=0., E=0.,fe=0., f=0.;
  int np=0;
  tgold(50.e-6,fQbeta,fbeta1f,0.001*fQbeta,2,em,fm);
  do{
     E=50.e-6+(fQbeta-50.e-6)*GetRandom();
     fe=fbeta1f.Eval(E);
     f=fm*GetRandom();
  }while(f>fe);
  if (fZdtr>=0.) np=3;
  if (fZdtr<0.)  np=2; 
  particle(np,E,E,0.,pi,0.,twopi,tcnuc,thnuc);
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::beta1fu(float Qbeta,float tcnuc,float thnuc,float c1,float c2,float c3,float c4)
{
  /*
   beta1fu simulates the angles and energy of beta particles emitted in beta decay of nucleus. 
   The decay is considered as 1st-forbidden unique. Its shape is product of theoretical spectrum 
   shape for allowed decay and two correction factors:
   1. theoretical of BJ'1969     cf1(e)=pnu**2+lamda2*pel**2,
      where lambda2 is the Coulomb function calculated in BJ'1969,
      and pel and pnu are impulses of electron and neutrino:
      pel=sqrt(w**2-1), pnu=(Qbeta-e)/emass , w=e/emass+1;
   2. empirical correction   cf2(e)=(1+c1/w+c2*w+c3*w**2+c4*w**3).  
   
   Values of the "lambda2" Coulomb function for some Zdtr values from: 
   H.Behrens, J.Janecke, "Numerical tables for beta-decay and electron capture", Berlin, Springer-Verlag, 1969.
   Values are calculated as product between unscreened lambda2 (Table II of BJ'1969) and screened corrections (Table III), and are given for "standard" values of momentum (0.1-50.0 in units of m_e*c). 
   (Log values of these momenta are in array plog69.) 
  */
  int i;
  for (i=0 ;i<50;i++){
    fSl[i]=0;
  }
  fQbeta=Qbeta;
  fC1=c1;
  fC2=c2;
  fC3=c3;
  fC4=c4;
  if (fZdtr==19){ // 39Ar, fQbeta=0.565;
     double Sl[]={2.0929, 1.2337, 1.0747, 1.0234, 0.99977, 0.98728, 0.98024, 0.97624, 0.97445, 0.97377, 0.97406, 0.97549, 0.9757, 0.9754, 0.9754, 0.9756, 0.9760};
     fSlSize=17;
     for (i=0;i<fSlSize;i++){
         fSl[i]=Sl[i];
     }
  }   
  else if (fZdtr==20){ // 42Ar, Qbeta=0.600 + 42K, fQbeta=3.525;
     double Sl[48]={2.2248, 1.2634, 1.0851, 1.0275, 1.0008, 0.98693, 0.97884, 0.97426, 0.97213, 0.97128, 0.97138, 0.97276, 0.9731, 0.9728, 0.9728, 0.9731, 0.9735, 0.9740, 0.9745, 0.9750, 0.9756, 0.9762, 0.9768, 0.9774, 0.9780, 0.9794, 0.9808, 0.9821, 0.9834, 0.9846, 0.9859, 0.9870, 0.9882, 0.9903, 0.9924}; 
     fSlSize=35;
     for (i=0;i<fSlSize;i++){
         fSl[i]=Sl[i];
     }
  }
  else if (fZdtr==39){ // 90Sr, fQbeta=0.546;
    double Sl[48]={5.6836, 2.0435, 1.3704, 1.1386, 1.0327, 0.97761, 0.94571, 0.92621, 0.91383, 0.90577, 0.89708, 0.89379, 0.89354, 0.89479, 0.89695, 0.89953, 0.90229};
     fSlSize=17;
     for (i=0;i<fSlSize;i++){
         fSl[i]=Sl[i];
     }
  }
  else if (fZdtr==40){ //90Y, fQbeta=2.228;
     double Sl[]={5.8992, 2.0922, 1.3883, 1.1454, 1.0345, 0.97692, 0.94344, 0.92294, 0.90998, 0.90153, 0.89243, 0.88892, 0.88848, 0.88970, 0.89186, 0.89454, 0.89739, 0.90037, 0.90330, 0.90631, 0.90931, 0.91223, 0.91507, 0.9174, 0.9195, 0.9246, 0.9295, 0.9343, 0.9388, 0.9432};
     fSlSize=30;
     for (i=0;i<fSlSize;i++){
         fSl[i]=Sl[i];
     }
  }
  else if (fZdtr==56){ // 137-Cs, Qbeta=0.514 , to level 0.662
     double Sl[]={9.3262, 2.8592, 1.6650, 1.2481, 1.0580, 0.95794, 0.89948, 0.86350, 0.84043, 0.82535, 0.80875, 0.80209, 0.80046, 0.80152, 0.80409, 0.80752, 0.81167 };
     fSlSize=17;
     for (i=0;i<fSlSize;i++){
         fSl[i]=Sl[i];
     }
  }
  
  TF1 fbeta1fu ("funbeta1fu",this,&Decay::funbeta1fu,0,fQbeta,3);
  float em=0.,fm=0., E,fe=0.;
  float f;
  int np;
  tgold(50.e-6,fQbeta,fbeta1fu,0.001*fQbeta,2,em,fm);
  do{
     E=50.e-6+(fQbeta-50.e-6)*GetRandom();
     fe=fbeta1fu.Eval(E);
     f=fm*GetRandom();
  }while(f>fe);
  if(fZdtr>=0.) np=3;
  if(fZdtr<0.)  np=2;
  particle(np,E,E,0.,pi,0.,twopi,tcnuc,thnuc);
  return;
}  

/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::beta2f(float Qbeta,float tcnuc,float thnuc,int kf,float c1,float c2,float c3,float c4)
{
  // beta2 simulates the angles and energy of beta particles emitted 
  // in beta decay of nucleus. The decay is considered as forbidden;
  // correction factor to the allowed spectrum shape has one of a form,
  // typical for unique k-forbidden spectra:
  // k=1: cf(e)=pel**2+c1*       pnu**2,
  // k=2: cf(e)=pel**4+c1*pel**2*pnu**2+c2*       pnu**4,
  // k=3: cf(e)=pel**6+c1*pel**4*pnu**2+c2*pel**2*pnu**4+c3*       pnu**6,
  // k=4: cf(e)=pel**8+c1*pel**6*pnu**2+c2*pel**4*pnu**4+c3*pel**2*pnu**6+c4*pnu**8,
  // where pel and pnu are impulses of electron and neutrino:
  // pel=sqrt(w**2-1), pnu=(Qbeta-e)/emass , w=e/emass+1.
  // VIT, 30.07.1992; 15.10.1995; 31.03.2006.
  fQbeta=Qbeta;
  fKf=kf;
  fC1=c1;
  fC2=c2;
  fC3=c3;
  fC4=c4;
  TF1 fbeta2f ("funbeta2f",this,&Decay::funbeta2f,0,fQbeta,3);
  float em=0.,fm=0., E=0.,fe=0., f=0.;
  int np=0;
  tgold(50.e-6,fQbeta,fbeta2f,0.001*fQbeta,2,em,fm);
  do{
     E=50.e-6+(fQbeta-50.e-6)*GetRandom();
     fe=fbeta2f.Eval(E);
     f=fm*GetRandom();
  }while(f>fe);
  if (fZdtr>=0.) np=3;
  if (fZdtr<0.)  np=2;
  particle(np,E,E,0.,pi,0.,twopi,tcnuc,thnuc);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::nucltransK(float Egamma,float Ebinde,float conve,float convp)
{
  /* NucltransK choise one of the three concurent processes
     by which the transition from one nuclear state to another is
     occured: gamma-ray emission, internal conversion and internal
     pair creation. Conversion electrons are emitted only with one fixed energy
     (usually with Egamma-E(K)_binding_energy).
     Egamma - gamma-ray energy (MeV) [=difference in state energies];
     Ebinde - binding energy of electron (MeV);
     conve  - internal electron conversion coefficient [=Nelectron/Ngamma];
     convp  - pair conversion coefficient [=Npair/Ngamma];
     tevst  - time of event's start (sec);
     fNbPart - current number of last particle in event;
     npgeant[fNbPart]     - GEANT number for particle identification
                          (1 for gamma, 2 for e+ and 3 for e-);
     fPmoment[1-3][fNbPart] - x,y,z components of particle momentum (MeV);
     fPtime[fNbPart]       - time shift from previous time to calculate
                           when particle was emitted (sec).
  */
   float p=(1.+conve+convp)*GetRandom();
   if (p<=1) {
      particle(1,Egamma,Egamma,0,pi,0,twopi,fTclev,fThlev);//gamma
   }
   else if (p<=1.+conve){
      particle(3,Egamma-Ebinde,Egamma-Ebinde,0.,pi,0.,twopi,fTclev,fThlev);//electron
      particle(1,Ebinde,Ebinde,0.,pi,0.,twopi,0,0);//gamma
   }
   else pair(Egamma-2.*emass);//e+e- pair
   return; 
}
/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::nucltransKL(float Egamma,float EbindeK,float conveK,float EbindeL,float conveL,float convp)
{
  // Subroutine nucltransKL choise one of the three concurent processes
  // by which the transition from one nuclear state to another is
  // occured: gamma-ray emission, internal conversion and internal
  // pair creation. Conversion electrons are emitted with two 
  // fixed energies
  // (Egamma-E(K)_binding_energy and Egamma-E(L)_binding_energy).
  // VIT, 5.07.1995.
  float p=(1.+conveK+conveL+convp)*GetRandom();
  if (p<=1.)particle(1,Egamma,Egamma,0.,pi,0.,twopi,fTclev,fThlev);//gamma
  else if (p<=1.+conveK){
     particle(3,Egamma-EbindeK,Egamma-EbindeK,0.,pi,0.,twopi,fTclev,fThlev);//electron
     particle(1,EbindeK,EbindeK,0.,pi,0.,twopi,0,0);//gamma
  }
  else if (p<=1.+conveK+conveL){
     particle(3,Egamma-EbindeL,Egamma-EbindeL,0.,pi,0.,twopi,fTclev,fThlev);//electron
     particle(1,EbindeL,EbindeL,0.,pi,0.,twopi,0,0);//gamma
  }
  else pair(Egamma-2.*emass);
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::nucltransKLM(float Egamma,float EbindeK,float conveK,float EbindeL,float conveL,float EbindeM,float conveM, float convp)
{
  /*
  Subroutine nucltransKLM choises one of the three concurent processes
  by which the transition from one nuclear state to another is
  occured: gamma-ray emission, internal conversion and internal
  pair creation. Conversion electrons are emitted with three fixed energies:
  Egamma-E(K)_binding_energy, Egamma-E(L)_binding_energy and Egamma-E(M)_binding_energy.
  Egamma  - gamma-ray energy (MeV) [=difference in state energies];
  EbindeK - binding energy of electron (MeV) on K-shell;
  conveK  - internal conversion coefficient [=Nelectron/Ngamma] from K-shell;
  EbindeL - binding energy of electron (MeV) on L-shell;
  conveL  - internal conversion coefficient [=Nelectron/Ngamma] from L-shell;
  EbindeM - binding energy of electron (MeV) on M-shell;
  conveM  - internal conversion coefficient [=Nelectron/Ngamma] from M-shell;
  convp   - pair conversion coefficient [=Npair/Ngamma];
  tevst   - time of event's start (sec);
  fNbPart  - current number of last particle in event;
  npgeant[fNbPart]      - GEANT number for particle identification
                        (1 for gamma, 2 for e+ and 3 for e-);
  fPmoment[1-3][fNbPart] - x,y,z components of particle momentum (MeV);
  fPtime[fNbPart]        - time shift from previous time to calculate
                         when particle was emitted (sec).
  */
 
  float p=(1.+conveK+conveL+conveM+convp)*GetRandom(); 
  if (p<=1.) particle(1,Egamma,Egamma,0.,pi,0.,twopi,fTclev,fThlev);//gamma
  else if (p<=1.+conveK){
     particle(3,Egamma-EbindeK,Egamma-EbindeK,0.,pi,0.,twopi,fTclev,fThlev);//electron
     particle(1,EbindeK,EbindeK,0.,pi,0.,twopi,0,0);//gamma
  }
  else if (p<=1.+conveK+conveL){
     particle(3,Egamma-EbindeL,Egamma-EbindeL,0.,pi,0.,twopi,fTclev,fThlev);//electron
     particle(1,EbindeL,EbindeL,0.,pi,0.,twopi,0,0);//gamma
  }
  else if (p<=1.+conveK+conveL+conveM){
     particle(3,Egamma-EbindeM,Egamma-EbindeM,0.,pi,0.,twopi,fTclev,fThlev);//electron
     particle(1,EbindeM,EbindeM,0.,pi,0.,twopi,0,0);//gamma
  }
  else pair(Egamma-2.*emass);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::nucltransKLM_Pb(float Egamma,float EbindeK,float conveK,float EbindeL,float conveL,float EbindeM,float conveM, float convp)
{
  /*
  The same as nucltransKLM but two X rays are emitted after K conversion
  in deexcitation of 208-Pb in decay 208Tl->208Pb. 
  VIT, 4.02.2009.
 
  Subroutine nucltransKLM choises one of the three concurent processes
  by which the transition from one nuclear state to another is
  occured: gamma-ray emission, internal conversion and internal
  pair creation. Conversion electrons are emitted with three fixed energies:
  Egamma-E(K)_binding_energy, Egamma-E(L)_binding_energy and 
  Egamma-E(M)_binding_energy).
  VIT, 4.01.2007.
 */
  float p, p1;
  p=(1.+conveK+conveL+conveM+convp)*GetRandom();
  if (p<=1.) 
     particle(1,Egamma,Egamma,0.,pi,0.,twopi,fTclev,fThlev);//gamma
  else if (p<=1.+conveK){
     particle(3,Egamma-EbindeK,Egamma-EbindeK,0.,pi,0.,twopi,fTclev,fThlev);//electron
     p1=100.*GetRandom();
     if (p1<=73.9){
        particle(1,0.074,0.074,0.,pi,0.,twopi,0,0);//gamma
        particle(1,0.014,0.014,0.,pi,0.,twopi,0,0);//gamma
     }  
     else{
        particle(1,0.085,0.085,0.,pi,0.,twopi,0,0);//gamma
        particle(1,0.003,0.003,0.,pi,0.,twopi,0,0);//gamma
        // in 4.8% few low energy particles are emitted; they are neglected
     }
  }
  else if (p<=1.+conveK+conveL){
     particle(3,Egamma-EbindeL,Egamma-EbindeL,0.,pi,0.,twopi,fTclev,fThlev);//electron
     particle(1,EbindeL,EbindeL,0.,pi,0.,twopi,0,0);//gamma
  }
  else if (p<=1.+conveK+conveL+conveM){
     particle(3,Egamma-EbindeM,Egamma-EbindeM,0.,pi,0.,twopi,fTclev,fThlev);//electron
     particle(1,EbindeM,EbindeM,0.,pi,0.,twopi,0,0);//gamma
  }
  else pair(Egamma-2.*emass);

  return;
  
}
/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::pair(float Epair)
{
  // Generation of e+e- pair in zero-approximation for  INTERNAL pair creation:
  //  1) energy of e+ is equal to the energy of e-;
  //  2) e+ and e- are emitted in the same direction.

  phi=twopi*GetRandom();
  float ctet=-1.+2.*GetRandom();
  float teta=acos(ctet);
  float E=0.5*Epair;
  particle(2,E,E,teta,teta,phi,phi,fTclev,fThlev);
  particle(3,E,E,teta,teta,phi,phi,0.,0.);
  return; 
}


/////////////////////////////////////////////////////////////////////////////////////////////
void Decay::tgold(float a,float b,TF1 &fb,float eps,int minmax,float&xextr,float &fextr)
{
  // Subroutine tgold determines maximum or minimum of the function f(x) in 
  // the interval [a,b] by the gold section method. 
  float qc=0.61803395;
  float xR=a+(b-a)*qc;
  float xL=b-(b-a)*qc;
  fb.SetParameters(fQbeta,emass,fZdtr);
  float yR=fb.Eval(xR);
  float yL=fb.Eval(xL);
  do{
     if (minmax==1){
        if (yR<yL){
           a=xL;
           xL=xR;
           yL=yR;
           xR=a+(b-a)*qc;
           yR=fb.Eval(xR);
        }
     }
     else{ 
        if (yR>yL){
           a=xL;
           xL=xR;
      yL=yR;
      xR=a+(b-a)*qc;
      yR=fb.Eval(xR); 
        }
     }
     b=xR;
     xR=xL;
     yR=yL;
     xL=b-(b-a)*qc;
     yL=fb.Eval(xL);
  }while(b-a>eps);
  xextr=0.5*(a+b);
  fextr=fb.Eval(xextr);
}

/////////////////////////////////////////////////////////////////////////////////////////////
complex<double> Decay::cgamma(complex<double> z)
{
  complex<double> g,z0,z1;
  double x0=0.,q1=0.,q2=0.,x=0.,y=0.,th=0.,th1=0.,th2=0.;
  double g0=0.,gr=0.,gi=0.,gr1=0.,gi1=0.;
  double na=0.,t=0.,x1=0.,y1=0.,sr=0.,si=0.;
  int j=0,k=0;

  static double a[] = { 8.333333333333333e-02, -2.777777777777778e-03,
                        7.936507936507937e-04, -5.952380952380952e-04,
                        8.417508417508418e-04, -1.917526917526918e-03,
                        6.410256410256410e-03, -2.955065359477124e-02,
                        1.796443723688307e-01, -1.39243221690590};

  x = real(z);
  y = imag(z);
  if (x > 171) return complex<double>(1e308,0);
  if ((y == 0.0) && (x == (int)x) && (x <= 0.0))  return complex<double>(1e308,0);
  else if (x < 0.0) {
        x1 = x;
        y1 = y;
        x = -x;
        y = -y;
  }
  x0 = x;
  if (x <= 7.0) {
     na = (int)(7.0-x);
     x0 = x+na;
  }
  q1 = sqrt(x0*x0+y*y);
  th = atan(y/x0);
  gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(twopi);
  gi = th*(x0-0.5)+y*log(q1)-y;
  for (k=0;k<10;k++){
      t = pow(q1,-1.0-2.0*k);
      gr += (a[k]*t*cos((2.0*k+1.0)*th));
      gi -= (a[k]*t*sin((2.0*k+1.0)*th));
  }
  if (x <= 7.0) {
     gr1 = 0.0;
     gi1 = 0.0;
     for (j=0;j<na;j++) {
         gr1 += (0.5*log((x+j)*(x+j)+y*y));
         gi1 += atan(y/(x+j));
     }
     gr -= gr1;
     gi -= gi1;
  }
  if (x1 < 0.0) {
     q1 = sqrt(x*x+y*y);
     th1 = atan(y/x);
     sr = -sin(pi*x)*cosh(pi*y);
     si = -cos(pi*x)*sinh(pi*y);
     q2 = sqrt(sr*sr+si*si);
     th2 = atan(si/sr);
     if (sr < 0.0) th2 += pi;
     gr = log(pi/(q1*q2))-gr;
     gi = -th1-th2-gi;
     x = x1;
     y = y1;
  }
     g0 = exp(gr);
     gr = g0*cos(gi);
     gi = g0*sin(gi);
  g = complex<double>(gr,gi);
  return g;
}

/////////////////////////////////////////////////////////////////////////////////////////////
double Decay::divdif (double xtab[], double xval)
{
  //////////////
  /////////////
  int i=0,j=0 ;
  double difval=0.;
  int N=fSlSize;
  if(N<2 ) return difval;
  int arr[]={2,10,fSlSize-1};
  int M=TMath::MinElement(3,arr);
  int mplus=M+1;
  int ix=0, mid;
  int iy=N;
  if (xtab[0]>xtab[N-1]){
     do{
        mid=(ix+iy)/2;
       if(xval<=xtab[mid]) ix=mid;
       else iy=mid;
     }while ((iy-ix)>1);    
  }
  else{
     do{
       mid=(ix+iy)/2;
       if(xval>=xtab[mid]) ix=mid;
       else iy=mid;
     }while((iy-ix)>1);
  }
  int npts=M+2-M%2;
  int ip=0;
  int l=0, isub;
  double T[50],D[50];
  for (int ii=0;ii<N;ii++) {T[ii]=0; D[ii]=0;}
  bool dothat=false; 
  do {
     if (dothat){
        l=-l;
        if(l>=0) l=l+1;
     }
     isub=ix+l;
     if (0<=isub && isub<=N){
        
        T[ip]=xtab[isub];
        D[ip]=fSl[isub];
        ip++;
     } 
     else{
        npts=mplus;
     }
     dothat=true;
  }while (ip<npts);
  bool extra=(npts!=mplus);
  bool next=false;
  for (int l=1;l<=M;l++){
      if (extra){
         isub=mplus-l;
         D[M+1]=(D[M+1]-D[M-1])/(T[M+1]-T[isub-1]);
         next=true;
      }
      if (!extra||next){
         i=mplus-1;

         for (int j=l;j<=M;j++){
             isub=i-l;
             
             D[i]=(D[i]-D[i-1])/(T[i]-T[isub]);
             i--;
         }
      }
  }
  double sum=D[mplus-1];
  if(extra) sum=0.5*(sum+D[M+1]);
  j=M-1;
  for (int l=1;l<=M;l++){
      sum=D[j]+(xval-T[j])*sum;
      j--;
  }
  difval=sum;
  return difval;
}

/************************************************/
float Decay::GetRandom()
{
  return rnd->Rndm();
}

/************************************************/


int main()
{
  
  Decay decay;  
  decay.rnd = new TRandom1();
  decay.GENBBdia();
  for (int i=0;i<decay.nevents;i++){
      decay.icurrent=i;
      decay.GENBBsub();
      if (decay.ier!=0){
         printf("\033[31m incorrect input parameter in GENBBsub\033[m \n***********\n ");
         decay.GENBBdia();
      }
  }
  printf("all events: %d\n",(int)(round(decay.nevents*decay.toallevents)));
return 1;
}
