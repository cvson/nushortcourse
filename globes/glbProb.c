/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * To obtain the oscillation probabilities for T2K experiment
 * Compile with ``make glbProb''
 * Modified: cvson@ifirse.icise.vn or cvson@utexas.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <list>
#include <complex>
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include <globes/globes.h>   /* GLoBES library */


/* If filename given, write to file; for empty filename write to screen */
char MYFILE[]="glbProb_4vson.root";


/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{
    
    /* Define standard oscillation parameters (cf. NuFIT 5.2 http://www.nu-fit.org/?q=node/256) */
    double true_theta12 = asin(sqrt(0.303));
    double true_theta13 = asin(sqrt(0.02203));
    double true_theta23 = M_PI/4;//assume maximal mixing
    double true_deltacp = M_PI*(-1.)/2.;//assume maximal CP violation
    double true_sdm = 7.41e-5;
    double true_ldm = 2.511e-3;//assume normal mass ordering
    
    
    /* Initialize libglobes */
    glbInit(argv[0]);
    
    /* Initialize reactor experiment */
    glbInitExperiment("T2K.glb",&glb_experiment_list[0],&glb_num_of_exps);
    // glbInitExperiment("Reactor2.glb",&glb_experiment_list[0],&glb_num_of_exps);
    
    /* Initialize parameters */
    glb_params true_values = glbAllocParams();
    glbDefineParams(true_values,true_theta12,true_theta13,true_theta23,true_deltacp,true_sdm,true_ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    
    
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    /* Compute chi^2 without correlations */
    //InitOutput(MYFILE1,"Format: Log(10,s22th13)   sigma_E       chi^2 \n");
    TFile *poutfile = new TFile(MYFILE,"RECREATE");
    const int NPOINTENE = 600;
    double aprobmu2mu[NPOINTENE];
    double aprobmu2muanti[NPOINTENE];
    double aprobmu2e[NPOINTENE];
    double aprobmu2eanti[NPOINTENE];
    
    double aenergy[NPOINTENE];
    double energy = 0.2;//GeV minimum energy
    Int_t nbinprob = 0;
    double prob;
    //get probabilities as function of energy
    while(energy < 3.2){//only check up 3.2GeV
        aenergy[nbinprob] = energy;
        //numu --> numu prob.
        prob = glbProfileProbability(0, 2, 2, 1, energy);
        aprobmu2mu[nbinprob] = prob;
        //numu --> nue prob.
        prob = glbProfileProbability(0, 2, 1, 1, energy);
        aprobmu2e[nbinprob] = prob;
        //numubar -->numubar prob.
        prob = glbProfileProbability(0, 2, 2, -1, energy);
        aprobmu2muanti[nbinprob] = prob;
        //numubar -->nuebar prob.
        prob = glbProfileProbability(0, 2, 1, -1, energy);
        aprobmu2eanti[nbinprob] = prob;
        
        energy += 0.005;//energy step
        nbinprob +=1;
    }
    
    TGraph *pgrProbmu2mu = new TGraph(NPOINTENE,aenergy,aprobmu2mu);
    TGraph *pgrProbmu2e = new TGraph(NPOINTENE,aenergy,aprobmu2e);
    
    TGraph *pgrProbmu2muanti = new TGraph(NPOINTENE,aenergy,aprobmu2muanti);
    TGraph *pgrProbmu2eanti = new TGraph(NPOINTENE,aenergy,aprobmu2eanti);
    
    //for looking at probabilities in 2D of dcp vs. sinsqtheta23
    double dcpmin = -1.;
    int NPOINT = 100;
    double dcpstep = 0.02;
    double dcpmax = dcpmin+(NPOINT)*dcpstep;
    
    double sinsqth23min = 0.3;
    double sinsqth23max = 0.7;
    int NPOINTTH23 = 100;
    double sinsqth23step = (sinsqth23max-sinsqth23min)/(1.*NPOINTTH23);
    //numu-->numu  prob.
    TH2D *hdcpvssinsqth23_mu2mu = new TH2D("hdcpvssinsqth23_mu2mu","",NPOINTTH23,sinsqth23min,sinsqth23max,NPOINT,dcpmin*M_PI,dcpmax*M_PI);
    //numu-->nue prob.
    TH2D *hdcpvssinsqth23_mu2e = new TH2D("hdcpvssinsqth23_mu2e","",NPOINTTH23,sinsqth23min,sinsqth23max,NPOINT,dcpmin*M_PI,dcpmax*M_PI);
    int county =0;
    int countx;
    double x,y;
    double dcp_test, th23_test;
    for (y = dcpmin ; y < dcpmax ; y=y+dcpstep){
        dcp_test = (y+dcpstep*0.5)*M_PI;// get dcp value
        glbSetOscParams(true_values,dcp_test,GLB_DELTA_CP);//reset the dcp parameter
        county +=1;
        countx = 0;//reset the count for the following loop
        for (x = sinsqth23min ; x < sinsqth23max ; x=x+sinsqth23step) {
            th23_test = asin(sqrt(x+sinsqth23step/2.)) ;
            countx += 1;
            glbSetOscParams(true_values,th23_test,GLB_THETA_23);//reset the theta23 parameter
            glbSetOscillationParameters(true_values);//reset the oscillation prameter
            glbSetRates();
            //numu-->numu prob.
            double tmpprobmu2mu = glbProfileProbability(0, 2, 2, 1, 0.6);
            hdcpvssinsqth23_mu2mu->SetBinContent(countx,county,tmpprobmu2mu);
            // numu-->nue prob.
            double tmpprobmu2e = glbProfileProbability(0, 2, 1, 1, 0.6);
            hdcpvssinsqth23_mu2e->SetBinContent(countx,county,tmpprobmu2e);
            
        }
        
    }
    
    //open the output file and write the histogram
    poutfile->cd();
    pgrProbmu2mu->Write("probmu2mu");
    pgrProbmu2e->Write("probmu2e");
    
    pgrProbmu2muanti->Write("probmu2mu_anti");
    pgrProbmu2eanti->Write("probmu2e_anti");
    
    hdcpvssinsqth23_mu2mu->Write("hdcpvssinsqth23_mu2mu");
    hdcpvssinsqth23_mu2e->Write("hdcpvssinsqth23_mu2e");
    
    
    
    poutfile->Close();
    
    /* Destroy parameter and projection vector(s) */
    glbFreeParams(true_values);
    
    
    exit(0);
}


