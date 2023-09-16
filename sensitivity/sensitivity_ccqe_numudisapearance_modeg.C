//
//  sensitivity_ccqe_numudisapearance.C
////////////////////////////////////////////////////
//
//  To make sensitivity on atmospheric parameter with simple numu disappearance
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
void titleStyle2D(TH2* h1);
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
//implement correct model
double model2flav_vac_MuSurvive(double Ene, double baseL, double theta_23, double dmsq_23){
    double sin_BL = sin(1.267*baseL*dmsq_23/Ene);//L/E term
    double sinsq_2th23 = theta_23;//pow(sin(2*theta_23),2);//mixing term
    double probmu2mu = 1 - sinsq_2th23*sin_BL*sin_BL;//survival prob.
    
    return probmu2mu;
}

void sensitivity_ccqe_numudisapearance_modeg(){
    TFile *pfile = new TFile("../eventpred/op_eventrate.root","READ");
    
    const int NCHAN4SHOW = 6;
    const int NDET = 2;
    TH1D *heventrate4show[NDET][NCHAN4SHOW];
    for (Int_t ichan=0; ichan<NCHAN4SHOW; ++ichan) {
        heventrate4show[0][ichan] = (TH1D*)pfile->Get(Form("ND_eventrate_chan%d",ichan));
        heventrate4show[1][ichan] = (TH1D*)pfile->Get(Form("FD_eventrate_chan%d",ichan));
    }
    //Number of event for normalization
    double NEVENTPred_nom = heventrate4show[1][1] ->Integral(0,heventrate4show[1][1]->GetNbinsX()+1 );
    //here we use the FD event rate CCQE only
    const int NEVENTMC = 10000;
    //normalization factor
    double normfactor = NEVENTPred_nom/NEVENTMC*1.0;
    cout<<"normalization factor "<<normfactor<<endl;
    
    TFile *popfile = new TFile("op_eventmc_4sensitivity_modeg.root","RECREATE");
    TTree *ptree = new TTree("Event","Neutrino event");
    double Enu, ProbMu2mu;
    
    ptree->Branch("enu",&Enu,"enu/D");
    //ptree->Branch("probmu2mu",&ProbMu2mu,"probmu2mu/D");
    
    const double T2KBASELINE = 295;//km
    
    const int NTHETA23BIN = 201;//50+1
    const double theta_MIN = 0.3;
    const double theta_MAX = 0.7;
    const double theta_BEST = 0.95;//45*Pi/180
    
    /*const int NDM23BIN = 51;//50+1
    const double dm23_MIN = 2e-3;
    const double dm23_MAX = 3e-3;*/
    const int NDM23BIN = 201;//50+1
    const double dm23_MIN = -3e-3;
    const double dm23_MAX = 3e-3;
    
    const double dm23_BEST = 2.4e-3;//eV
    
    TH1D *hpred_osc[NDM23BIN][NTHETA23BIN];
    for (int idm=0; idm<NDM23BIN; ++idm) {
        for (int itheta=0; itheta<NTHETA23BIN; ++itheta) {
            hpred_osc[idm][itheta] = new TH1D(Form("hpred_osc_dmbin%d_thetabin%d",idm,itheta),"",100,0,5);
        }
    }
    TH1D *hpred_noosc = new TH1D("hpred_noosc","",100,0,5);
    TH1D *hdata = new TH1D("hdata","",100,0,5);
    
    for (int ievmc=0; ievmc<NEVENTMC; ++ievmc) {
        if(ievmc%1000==0) cout<<"processing "<<ievmc<<"/"<<NEVENTMC<<endl;
        Enu = heventrate4show[1][1] -> GetRandom();
        //ProbMu2mu = model2flav_vac_MuSurvive
        ptree->Fill();
        hpred_noosc->Fill(Enu);
        
        for (int idm=0; idm<NDM23BIN; ++idm) {
            double dm_val = dm23_MIN + (dm23_MAX-dm23_MIN)*(idm+0.5)/(NDM23BIN-1);
            for (int itheta=0; itheta<NTHETA23BIN; ++itheta) {
                double theta_val_tmp = theta_MIN + (theta_MAX-theta_MIN)*(itheta+0.5)/(NTHETA23BIN-1);//this is theta
                double theta_val = 4*theta_val_tmp*(1-theta_val_tmp);//this is sin^22theta
                ProbMu2mu = model2flav_vac_MuSurvive(Enu,T2KBASELINE,theta_val,dm_val);
                hpred_osc[idm][itheta]->Fill(Enu,ProbMu2mu);
            }
        }
        
        //best fit
        ProbMu2mu = model2flav_vac_MuSurvive(Enu,T2KBASELINE,theta_BEST,dm23_BEST);
        hdata->Fill(Enu, ProbMu2mu);
        
    }
    //take normalization
    hdata->Scale(normfactor);
    hpred_noosc->Scale(normfactor);
    for (int idm=0; idm<NDM23BIN; ++idm) {
        for (int itheta=0; itheta<NTHETA23BIN; ++itheta) {
            hpred_osc[idm][itheta]->Scale(normfactor);
            hpred_osc[idm][itheta]->Write(Form("hpred_osc_dmbin%d_thetabin%d",idm,itheta));
        }
    }
    
    hdata->Write("hdata");
    hpred_noosc->Write("hpred_noosc");
    
    //surface of chisq
    TH2D *hsurface = new TH2D("hsurface","",NTHETA23BIN-1,theta_MIN,theta_MAX,NDM23BIN-1,dm23_MIN,dm23_MAX);
    for (int idm=1; idm<NDM23BIN; ++idm) {
        for (int itheta=1; itheta<NTHETA23BIN; ++itheta) {
            double likelihood = LogLikelihood(hpred_osc[idm-1][itheta-1],hdata);
            hsurface->SetBinContent(itheta,idm,likelihood);
        }
    }
    hsurface->Write("hsurface");
    
    //plot the contour
    const int nlevels = 3;
    //http://www.reid.ai/2012/09/chi-squared-distribution-table-with.html
    double levels[nlevels] = {2.30, 4.61,6.18};//68%, 90% and 2 sigma
    //get minimum
    
    //double levels[nlevels] = {2.71, 6.63,9.0};//90%, 99%, 3 sigma

    new TCanvas;
    const Int_t NRGBs = 5;
      const Int_t NCont = 255;
      
      Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
      Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
      Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
      Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
      TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
      gStyle->SetNumberContours(NCont);
    //gStyle->SetPalette(kOcean);
    hsurface->GetXaxis()->SetTitle("sin^{2}#theta");
       hsurface->GetYaxis()->SetTitle("#Delta m^{2} [eV^{2}/c^{4}]");
       hsurface->GetYaxis()->SetTitleOffset(1.2);
       hsurface->GetZaxis()->SetRangeUser(0,25);
    gStyle->SetOptStat(0);
    titleStyle2D(hsurface);
    hsurface->Draw("colz");
    gPad->Print("sensitivesurface_modeg.pdf");
    //
    ptree->Write();
    popfile->Close();

}


void titleStyle2D(TH2* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.1);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.1);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.1);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.1);
    
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetXaxis()->SetTitleOffset(1.0);
}
