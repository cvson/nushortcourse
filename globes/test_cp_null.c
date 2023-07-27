/* GLoBES -- General LOng Baseline Experiment Simulator */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "TFile.h"
#include "TH1.h"

#include <globes/globes.h>   /* GLoBES library */
//#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
//without smearing
/*char MYFILE[]="t2k2_final_nosmear_sensi_cp.root";
char AEDLFILE[]="t2k2_final_nosmear.glb";*/
char MYFILE[]="t2k2_final_wsmear_sensi_cp.root";
char AEDLFILE[]="t2k2_final_wsmear.glb";

int main(int argc, char *argv[])
{ 
    /* Initialize libglobes */
    glbInit(argv[0]);
    
    /* Initialize experiment NFstandard.glb */
    glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps);
    
    /* Intitialize output */
    TFile *poutfile = new TFile(MYFILE,"RECREATE");
    
    
    
    /* Define standard oscillation parameters */
    //http://www.nu-fit.org/?q=node/211 fit wo atmospheric
    double theta12 = asin(sqrt(0.846))/2.0;
    double theta13 = asin(sqrt(0.085))/2.0;
    double theta23 = M_PI/4;
    //double theta23 = M_PI*49.5/(4*45.0);
    //double theta23 = M_PI*43.5/(4*45.0);
    double deltacp = 0;
    double sdm = 7.53e-5;
    double ldm = 2.509e-3;
    
    double delta_true ;
    double delta_test ;
    
    //to define the mass hierarchy
    double signdm13 = -1. ;
    
    int count, i ;
    
    //locate 50 points of delta
    const int NDELTApoint = 50;//between [-M_PI, M_PI]
    double deltastep = 2./NDELTApoint;// in unit of M_PI
    float	delta[NDELTApoint] ;
    
    //chisquare array
    double chi0n_sys[NDELTApoint],chipn_sys[NDELTApoint],chi0i_sys[NDELTApoint],chipi_sys[NDELTApoint], chimin_sys ;
    double chi0n_cor[NDELTApoint],chipn_cor[NDELTApoint],chi0i_cor[NDELTApoint],chipi_cor[NDELTApoint], chimin_cor ;
    
    
    double thetheta13,x,y, res1,res2;
    
    double chimh_sys, chimh_cor, chimh_sys_min, chimh_cor_min;
    double chi_sys_min, chi_cor_min, chi_sys_min_mh, chi_cor_min_mh;
    
    chi_sys_min = 999. ;
    chi_cor_min = 999. ;
    chi_sys_min_mh = 999. ;
    chi_cor_min_mh = 999. ;
    
    /* Initialize parameter and projection vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params fit_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params minimum = 	glbAllocParams();
    glb_params central_values = glbAllocParams();
    
    //glb_projection theta13_projection = glbAllocProjection();
    
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbSetCentralValues(true_values);
    
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    
    /* Set starting values and input errors for all projections */
    //this to be understand
    glbDefineParams(input_errors,theta12*0.03,theta13*0.03,theta23*0.11,0.0,sdm*0.03,ldm*0.05);
    //	glbDefineParams(input_errors,theta12*0.0003,theta13*0.0003,theta23*0.0011,0.0,sdm*0.0003,ldm*0.0003);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);
    glbSetInputErrors(input_errors);
    
    /* Set true value delta = 0. , NH */
    delta_true = 0.*M_PI ;
    glbDefineParams(true_values,theta12,theta13,theta23,delta_true,sdm,ldm);
    glbSetOscParams(true_values,+1.*ldm,GLB_DM_31);
    glbSetOscParams(true_values,delta_true,GLB_DELTA_CP);
    glbSetCentralValues(true_values);
    
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    TH1D *hchisq_truedcp0_nh_sys = new TH1D("hchisq_truedcp0_nh_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcp0_nh_proj= new TH1D("hchisq_truedcp0_nh_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcp0_ih_sys = new TH1D("hchisq_truedcp0_ih_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcp0_ih_proj= new TH1D("hchisq_truedcp0_ih_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    // Loop over test values of delta for both normal and inverted hierarchy
    count = 0;
    for (y = -1. ; y < 1.0 ; y=y+deltastep){
        delta_test = (y+0.5*deltastep)*M_PI ;
        // for normal hierarchy and delta
        glbDefineParams(test_values,theta12,theta13,theta23,delta_test,sdm,1.*ldm);
        glbSetOscParams(test_values,1.*ldm,GLB_DM_31);
        glbSetOscParams(test_values,delta_test,GLB_DELTA_CP);
        glbSetOscParams(test_values,theta13,GLB_THETA_13);
        
        //chisqr with systematic only
        chi0n_sys[count]=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        if(chi0n_sys[count]<0.){chi0n_sys[count]=0. ;}
        hchisq_truedcp0_nh_sys->SetBinContent(count+1,chi0n_sys[count]);
        
        //to make project on dCP, need to put input error from disappearance?
        glbDefineParams(input_errors,theta12*0.03,theta13*0.03,theta23*0.11,0,sdm*0.03,ldm*0.05);
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        glbSetCentralValues(test_values);
        
        //projection on ChiDelta
        chi0n_cor[count]=glbChiDelta(test_values,minimum,GLB_ALL);
        if(chi0n_cor[count]<0.){chi0n_cor[count]=0. ;}
        hchisq_truedcp0_nh_proj->SetBinContent(count+1, chi0n_cor[count]);
        
        //for invered hierarchy; why need to invert delta as well?
        glbDefineParams(test_values,theta12,theta13,theta23,delta_test,sdm,-1.*ldm);
        glbSetOscParams(test_values,-1.*ldm,GLB_DM_31);
        glbSetOscParams(test_values,-1.*delta_test,GLB_DELTA_CP);
        glbSetOscParams(test_values,theta13,GLB_THETA_13);
        
        chi0i_sys[count]=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        if(chi0i_sys[count]<0.){chi0i_sys[count]=0. ;}
        hchisq_truedcp0_ih_sys->SetBinContent(count+1,chi0i_sys[count]);
        
        //projection on ChiDelta
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        glbSetCentralValues(test_values);
        chi0i_cor[count]=glbChiDelta(test_values,minimum,GLB_ALL);
        if(chi0i_cor[count]<0.){chi0i_cor[count]=0. ;}
        hchisq_truedcp0_ih_proj->SetBinContent(count+1, chi0i_cor[count]);
        
        delta[count]=y;
        
        printf("Computing w/ dcp=0 true: %d/%d  points \n",count, NDELTApoint) ;
        count++;
        	
    }
    
    
    /* Set true value delta = pi , NH */
    delta_true = 1.*M_PI ;
    glbDefineParams(true_values,theta12,theta13,theta23,delta_true,sdm,1.*ldm);
    glbSetOscParams(true_values,+1.*ldm,GLB_DM_31);
    
    glbSetOscParams(true_values,delta_true,GLB_DELTA_CP);
    glbSetCentralValues(true_values);
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    // Loop over test values of delta for both normal and inverted hierarchy
    count = 0;
    TH1D *hchisq_truedcppi_nh_sys = new TH1D("hchisq_truedcppi_nh_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcppi_nh_proj= new TH1D("hchisq_truedcppi_nh_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcppi_ih_sys = new TH1D("hchisq_truedcppi_ih_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcppi_ih_proj= new TH1D("hchisq_truedcppi_ih_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    
    for (y = -1. ; y < 1.0 ; y=y+deltastep){
        delta_test = (y+0.5*deltastep)*M_PI ;
        glbDefineParams(test_values,theta12,theta13,theta23,delta_test,sdm,ldm);
        glbSetOscParams(test_values,1.*ldm,GLB_DM_31);
        glbSetOscParams(test_values,delta_test,GLB_DELTA_CP);
        glbSetOscParams(test_values,theta13,GLB_THETA_13);
        
        chipn_sys[count]=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        if(chipn_sys[count]<0.){chipn_sys[count]=0. ;}
        hchisq_truedcppi_nh_sys->SetBinContent(count+1,chipn_sys[count]);
        
        //projection on ChiDelta
        glbDefineParams(input_errors,theta12*0.03,theta13*0.03,theta23*0.11,0,sdm*0.03,ldm*0.05);
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        glbSetCentralValues(test_values);
        
        chipn_cor[count]=glbChiDelta(test_values,minimum,GLB_ALL);
        if(chipn_cor[count]<0.){chipn_cor[count]=0. ;}
        hchisq_truedcppi_nh_proj->SetBinContent(count+1,chipn_cor[count]);
        
        //chisqr with systematic only
        glbDefineParams(test_values,theta12,theta13,theta23,delta_test,sdm,-1.*ldm);
        glbSetOscParams(test_values,-1.*ldm,GLB_DM_31);
        glbSetOscParams(test_values,-1.*delta_test,GLB_DELTA_CP);
        glbSetOscParams(test_values,theta13,GLB_THETA_13);
        
        chipi_sys[count]=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        if(chipi_sys[count]<0.){chipi_sys[count]=0. ;}
        hchisq_truedcppi_ih_sys->SetBinContent(count+1,chipi_sys[count]);
        
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        glbSetCentralValues(test_values);
        
        chipi_cor[count]=glbChiDelta(test_values,minimum,GLB_ALL);
        if(chipi_cor[count]<0.){chipi_cor[count]=0. ;}
        hchisq_truedcppi_ih_proj->SetBinContent(count+1,chipi_cor[count]);
        
        printf("Computing w/ dcp=pi true: %d/%d  points \n",count, NDELTApoint) ;
        count++;
    }
    
    //Test for minimum chi^2 for each value of delta
    TH1D *hchisq_min_glob_sys = new TH1D("hchisq_min_glob_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_min_nh_nh_sys= new TH1D("hchisq_min_nh_nh_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_min_glob_proj = new TH1D("hchisq_min_glob_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_min_nh_proj= new TH1D("hchisq_min_nh_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    
    for (i = 0; i < NDELTApoint ; i++){
        chi_sys_min = 999. ;
        chi_cor_min = 999. ;
        chi_sys_min_mh = 999. ;
        chi_cor_min_mh = 999. ;
        //global minimum
        if(chi0n_sys[i]<chi_sys_min){chi_sys_min=chi0n_sys[i] ;}
        if(chipn_sys[i]<chi_sys_min){chi_sys_min=chipn_sys[i] ;}
        if(chi0i_sys[i]<chi_sys_min){chi_sys_min=chi0i_sys[i] ;}
        if(chipi_sys[i]<chi_sys_min){chi_sys_min=chipi_sys[i] ;}
        hchisq_min_glob_sys->SetBinContent(i+1,sqrt(chi_sys_min));
        
        //if mh is known
        if(chi0n_sys[i]<chi_sys_min_mh){chi_sys_min_mh=chi0n_sys[i] ;}
        if(chipn_sys[i]<chi_sys_min_mh){chi_sys_min_mh=chipn_sys[i] ;}
        hchisq_min_nh_nh_sys->SetBinContent(i+1,sqrt(chi_sys_min_mh));
        
        //global mimum of projection
        if(chi0n_cor[i]<chi_cor_min){chi_cor_min=chi0n_cor[i] ;}
        if(chipn_cor[i]<chi_cor_min){chi_cor_min=chipn_cor[i] ;}
        if(chi0i_cor[i]<chi_cor_min){chi_cor_min=chi0i_cor[i] ;}
        if(chipi_cor[i]<chi_cor_min){chi_cor_min=chipi_cor[i] ;}
        hchisq_min_glob_proj->SetBinContent(i+1,sqrt(chi_cor_min));
        
        //if MH is known
        if(chi0n_cor[i]<chi_cor_min_mh){chi_cor_min_mh=chi0n_cor[i] ;}
        if(chipn_cor[i]<chi_cor_min_mh){chi_cor_min_mh=chipn_cor[i] ;}
        hchisq_min_nh_proj->SetBinContent(i+1,sqrt(chi_cor_min_mh));
        
        
    }
//output
    poutfile->cd();
    hchisq_truedcp0_nh_sys ->Write();
    hchisq_truedcp0_nh_proj ->Write();
    hchisq_truedcp0_ih_sys ->Write();
    hchisq_truedcp0_ih_proj ->Write();
    
    hchisq_truedcppi_nh_sys ->Write();
    hchisq_truedcppi_nh_proj ->Write();
    hchisq_truedcppi_ih_sys ->Write();
    hchisq_truedcppi_ih_proj ->Write();
    
    hchisq_min_glob_sys ->Write();
    hchisq_min_nh_nh_sys ->Write();
    hchisq_min_glob_proj ->Write();
    hchisq_min_nh_proj ->Write();
    
    
    poutfile->Close();
    glbFreeParams(true_values);
    glbFreeParams(test_values);
    glbFreeParams(input_errors);
    
    exit(0);
}
