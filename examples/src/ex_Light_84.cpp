
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "../../include/DVR++.h"

using namespace std;

#define re 1.75  // in bohr
#define De 5.726  // in eV
#define alpha 1.22 // in bohr-1
#define mu 0.95    // in dalton
#define amu 1822.88863 // amu to au
#define hartree2eV  27.21138386 // hartree to eV


double pot(double r) {
  double PP;
  double fact;
  fact = exp(-alpha*(r-re));
  PP = De*fact*(fact-2.);
  return PP;
}

double analytic(int n) {
  double freq;
  double xe;
  double energy;
  double nn;

  freq = sqrt(2.*De/hartree2eV/mu/amu)*alpha*hartree2eV;  // in eV
  xe = 0.25*freq/De;
 
  nn = (double)n+0.5;
  energy = -De + freq*nn*(1. - xe*nn);
  
  return energy;
}


int main (int argc, char * argv[]) {

  DVR*  Dvr = new DVR();

  int    np = 12;
  
  double r0 = re;
  double freq0 = sqrt(2*De/mu)*alpha;

  double omega = freq0/sqrt(amu*hartree2eV);
  double delta = sqrt(1./mu/amu/omega);
  double fact = hartree2eV/mu/amu;

  Dvr->SetDimension(1);
  Dvr->AddHermBasis(0,np,fact,delta,r0);
  Dvr->PrepareDvr();
  
  for(int i=0;i<np;i++) {
    double QQ;
    Dvr->GetNodeCoord(i,&QQ);
    double Potential = pot(QQ);
    Dvr->LoadPotential(i,Potential);
  }
  
  Dvr->SolveHamiltonianLapack();

  for(int i=0;i<6;i++) {
    printf("%3d %16.4f %16.4f %16.4f\n",i,analytic(i),Dvr->GetEnergy(i),(Dvr->GetEnergy(i)-analytic(i)));
  }

  return EXIT_SUCCESS;
  
}
