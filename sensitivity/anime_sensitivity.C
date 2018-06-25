//
//  pred3_ccqe_numudisapearance.C
// using ouput from sensitivity_ccqe_numudisapearance.C
////////////////////////////////////////////////////
//
//  To make illustration for
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
double LogLikelihood(TH1* exp, TH1* obs)
{
    double llh_val = 0;
    // Hardcoded number of bins...
    for(int n = 1; n <= exp->GetNbinsX(); ++n){
        const double e_binval = exp->GetBinContent(n);
        const double o_binval = obs->GetBinContent(n);
       // if(e || o) ret += -2*log(TMath::Poisson(o,e));
        llh_val += 2*(e_binval-o_binval);
        if (o_binval!=0) {
            llh_val += -2*o_binval*log(e_binval/o_binval);
        }
    }
    return llh_val;
}


// formula for two-flavor oscillation in vacuum
//input theta_23 as sin^2 2theta 23 for convenience
double model2flav_vac_MuSurvive(double Ene, double baseL, double theta_23, double dmsq_23){
    double sin_BL = sin(1.267*baseL*dmsq_23/Ene);//L/E term
    double sinsq_2th23 = theta_23;//pow(sin(2*theta_23),2);//mixing term
    double probmu2mu = 1 - sinsq_2th23*sin_BL*sin_BL;//survival prob.
    
    return probmu2mu;
}

void anime_sensitivity(){
    
    
    TFile *pfile = new TFile("op_eventmc_4sensitivity.root","READ");
    const int NTHETA23BIN = 51;
    const double theta_MIN = 0.5;
    const double theta_MAX = 1;
    const double theta_BEST = 0.95;//45*Pi/180
    
    const int NDM23BIN = 51;
    const double dm23_MIN = 2e-3;
    const double dm23_MAX = 3e-3;
    const double dm23_BEST = 2.4e-3;//eV
    
    const double EnergyRange2Draw = 3.0;//use for drawing
    
    TH1D *hpred_osc[NDM23BIN][NTHETA23BIN];
    for (int idm=0; idm<NDM23BIN; ++idm) {
        for (int itheta=0; itheta<NTHETA23BIN; ++itheta) {
            hpred_osc[idm][itheta] = (TH1D*)pfile->Get(Form("hpred_osc_dmbin%d_thetabin%d",idm,itheta));
        }
    }
    
    TH1D *hpred_noosc = (TH1D*)pfile->Get("hpred_noosc");
    TH1D *hdata = (TH1D*)pfile->Get("hdata");
    
    gBenchmark->Start("hsen");
    gSystem->Unlink("hsensitivity_anime.gif"); // delete old file
    //condition to print
    const Int_t kUPDATE = 50;
    
    TCanvas *c1 = new TCanvas("c1","",0,0,1400,700);
    c1->Divide(2,1);
    c1->cd(1);
    gStyle->SetOptStat(0);
    hdata->SetMarkerStyle(8);
    hdata->SetMarkerColor(1);
    hdata->GetYaxis()->SetRangeUser(0,80);
    hdata->GetXaxis()->SetRangeUser(0,EnergyRange2Draw);
    hdata->GetYaxis()->SetTitle("Number of Events");
    hdata->GetXaxis()->SetTitle("Neutrino Energy [GeV]");

    hdata->Draw("ep");
    hdata->GetYaxis()->SetTitleOffset(1.5);
    hpred_noosc->SetLineWidth(2);
    hpred_noosc->SetLineColor(13);
    hpred_noosc->Draw("hist same");
    //to plot additional point
    int dmBIN_atBEST = (dm23_BEST-dm23_MIN)*NDM23BIN/(dm23_MAX-dm23_MIN)-0.5;
    int thBIN_atBEST = (theta_BEST-theta_MIN)*NTHETA23BIN/(theta_MAX-theta_MIN)-0.5;
    cout<<"BINS at best fit point for theta "<<thBIN_atBEST<<" for dm "<<dmBIN_atBEST<<endl;
    
    //emty histogram
    c1->cd(2);
     TGaxis::SetMaxDigits(3);
    TH2D *hsurface = new TH2D("hsurface","",NTHETA23BIN-1,theta_MIN,theta_MAX,NDM23BIN-1,dm23_MIN,dm23_MAX);
    hsurface->GetXaxis()->SetTitle("sin^{2}2#theta");
    hsurface->GetYaxis()->SetTitle("#Delta m^{2} [eV^{2}/c^{4}]");
    hsurface->GetYaxis()->SetTitleOffset(1.2);
    hsurface->Draw("AXIS");
    hsurface->GetZaxis()->SetRangeUser(0,25);
    Int_t nhist2plots=0;
    Int_t totalist = (NDM23BIN-1)* (NTHETA23BIN-1);
    for (int idm=1; idm<NDM23BIN; ++idm) {
        for (int itheta=1; itheta<NTHETA23BIN; ++itheta) {
            ++nhist2plots;
            c1->cd(1);
            int colorindex = (idm*itheta+1)%10;
            hpred_osc[idm-1][itheta-1]->SetLineColor(colorindex);
            hpred_osc[idm-1][itheta-1]->Draw("hist same");
            c1->cd(2);
            double likelihood = LogLikelihood(hpred_osc[idm-1][itheta-1],hdata);
            hsurface->SetBinContent(itheta,idm,likelihood);
            hsurface->Draw("colz");
            if(nhist2plots && (nhist2plots%kUPDATE==0)){
                c1->cd(1);
                hdata->Draw("ep same");
                c1->Update();
                if (gROOT->IsBatch()) {
                    c1->Print("hsensitivity_anime.gif+");
                    printf("i = %d/%d\n", nhist2plots,totalist);
                }
            }
            
        }
    }
    
    
 
    
    if (gROOT->IsBatch()) c1->Print("hsensitivity_anime.gif++");
    gBenchmark->Show("hsen");
    
    

    
  

    
    
}
