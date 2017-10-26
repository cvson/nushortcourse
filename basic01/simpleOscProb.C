//
//  simpleOscProb.C
////////////////////////////////////////////////////
//
//  Simple Oscillation Prob.
//
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
#include "TMath.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

//In vacuum
//Two-flavor vacuum function
double model2flav_vac_MuToTau(double *x);
double model2flav_vac_MuSurvive(double *x);
double model2flav_vac_MuToElec(double *x);

//===================================
//Two-flavor vacuum function
//Mu2Tau
double model2flav_vac_MuToTau(double *x)
{
    double E = x[0]; //energy
    
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double sin_BL = sin(1.267*L*dmsq_23/E);

    double p1 = sinsq_2th23*sinsq_2th23*sin_BL*sin_BL;
    return p1;
}

//MuToMu

double model2flav_vac_MuSurvive(double *x)
{
    double p1 = model2flav_vac_MuToTau(x);
    return 1.-p1;
}



//MuToElec
double model2flav_vac_MuToElec(double *x)
{
    return 0.0;
}
