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

void pred3_ccqe_numudisapearance(){
    
    
    TFile *pfile = new TFile("op_eventmc_4sensitivity.root","READ");
    const int NTHETA23BIN = 50;
    const double theta_MIN = 0.5;
    const double theta_MAX = 1;
    const double theta_BEST = 0.95;//45*Pi/180
    
    const int NDM23BIN = 50;
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
    
    int ith_th23bin = 1;
    int ith_dm32bin = dmBIN_atBEST;
    double dm_val1 = dm23_MIN + (dm23_MAX-dm23_MIN)*(ith_dm32bin+0.5)/NDM23BIN;
    double theta_val1 = theta_MIN + (theta_MAX-theta_MIN)*(ith_th23bin+0.5)/NTHETA23BIN;
    cout <<"prediction at dm "<<dm_val1<<" theta "<<theta_val1<<endl;
    TH1D* hpred_osc1 = (TH1D*)hpred_osc[ith_dm32bin][ith_th23bin]->Clone("hpred_osc1");
    hpred_osc1->SetLineColor(2);
    hpred_osc1->SetLineWidth(2);
    hpred_osc1->Draw("hist same");
    
    ith_th23bin = thBIN_atBEST;
    double dm_val2 = dm23_MIN + (dm23_MAX-dm23_MIN)*(ith_dm32bin+0.5)/NDM23BIN;
    double theta_val2 = theta_MIN + (theta_MAX-theta_MIN)*(ith_th23bin+0.5)/NTHETA23BIN;
    cout <<"prediction at dm "<<dm_val2<<" theta "<<theta_val2<<endl;
    TH1D* hpred_osc2 = (TH1D*)hpred_osc[ith_dm32bin][ith_th23bin]->Clone("hpred_osc2");
    hpred_osc2->SetLineColor(4);
    hpred_osc2->SetLineWidth(2);
    hpred_osc2->Draw("hist same");
    
    ith_dm32bin = 49;
    double dm_val3 = dm23_MIN + (dm23_MAX-dm23_MIN)*(ith_dm32bin+0.5)/NDM23BIN;
    double theta_val3 = theta_MIN + (theta_MAX-theta_MIN)*(ith_th23bin+0.5)/NTHETA23BIN;
    cout <<"prediction at dm "<<dm_val3<<" theta "<<theta_val3<<endl;
    TH1D* hpred_osc3 = (TH1D*)hpred_osc[ith_dm32bin][ith_th23bin]->Clone("hpred_osc3");
    hpred_osc3->SetLineColor(8);
    hpred_osc3->SetLineWidth(2);
    hpred_osc3->Draw("hist same");
    
    TLegend* leg0 = new TLegend(.42, .58, 0.85, .82);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(18);
    leg0->SetTextFont(43);
    leg0->AddEntry(hdata, "Data", "ep");
    leg0->AddEntry(hpred_noosc, "Unoscillated Prediction", "l");
    leg0->AddEntry(hpred_osc1, Form("(#Delta m^{2}=%.5f,sin^{2}2#theta=%.2f)",dm_val1,theta_val1), "l");
    leg0->AddEntry(hpred_osc2, Form("(#Delta m^{2}=%.5f,sin^{2}2#theta=%.2f)",dm_val2,theta_val2), "l");
    leg0->AddEntry(hpred_osc3, Form("(#Delta m^{2}=%.5f,sin^{2}2#theta=%.2f)",dm_val3,theta_val3), "l");
    leg0->Draw();
    
    c1->cd(2);
    gStyle->SetOptStat(0);
    TH1D* hratio = (TH1D*)hdata->Clone("hratio");
    hratio->Sumw2();
    hratio->Divide(hpred_noosc);
    
    TH1D* hratio1 = (TH1D*)hpred_osc1->Clone("hratio1");
    hratio1->Sumw2();
    hratio1->Divide(hpred_noosc);
    
    TH1D* hratio2 = (TH1D*)hpred_osc2->Clone("hratio1");
    hratio2->Sumw2();
    hratio2->Divide(hpred_noosc);
    
    TH1D* hratio3 = (TH1D*)hpred_osc3->Clone("hratio1");
    hratio3->Sumw2();
    hratio3->Divide(hpred_noosc);
    
    
    for (Int_t ibin=1; ibin<=hratio->GetXaxis()->GetNbins(); ibin++) {
        Float_t err1;
        if (hpred_noosc->GetBinContent(ibin)!=0) {
            err1 = sqrt(hdata->GetBinContent(ibin))/hpred_noosc->GetBinContent(ibin);
            hratio->SetBinError(ibin,err1);
        }
      
    }
    hratio->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    hratio->GetYaxis()->SetTitle("Data / Unoscillated Pred.");
    
    hratio->GetYaxis()->SetRangeUser(0,2.);
    hratio->Draw("ep");
    hratio1->Draw("hist same");
    hratio2->Draw("hist same");
    hratio3->Draw("hist same");
    hratio->GetYaxis()->SetTitleOffset(1.5);
    TLine *poneline = new TLine(0,1.,EnergyRange2Draw,1.0);
    poneline->SetLineStyle(7);
    poneline->Draw();
    c1->Print("pred3_datavsnoosc.eps");
    
    
    
    
    //surface of chisq
    /*TH2D *hsurface = new TH2D("hsurface","",NTHETA23BIN,theta_MIN,theta_MAX,NDM23BIN,dm23_MIN,dm23_MAX);
    for (int idm=0; idm<NDM23BIN; ++idm) {
        for (int itheta=0; itheta<NTHETA23BIN; ++itheta) {
            double likelihood = LogLikelihood(hpred_osc[idm][itheta],hdata);
            hsurface->SetBinContent(itheta,idm,likelihood);
        }
    }*/
    
  

    
    
}
