//
//  numu2num_allowedregion.C
////////////////////////////////////////////////////
//
//  To make allowed region for numu-> numu disappearance
//  root -b -q numu2num_allowedregion.C
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
double model2flav_vac_MuSurvive(double Ene, double baseL, double theta_23, double dmsq_23){
    double sin_BL = sin(1.267*baseL*dmsq_23/Ene);//L/E term
    double sinsq_2th23 = theta_23;//pow(sin(2*theta_23),2);//mixing term
    double probmu2mu = sinsq_2th23*sin_BL*sin_BL;//survival prob.
    
    return probmu2mu;
}
//should include the probability
//if near measured value, it's likelihood

void numu2numu_allowedregion(){
    gROOT->ProcessLine(".x ../rootlogon.C");
    const Bool_t isSmoothGraph = true;
    //p = 1-sin^22\theta_23 sin^2
    //calculate from T2K best fit values
    double sinsq_th23_best = 0.536;
    double sinsq_th23_err = 0.005;
    double dmsq_23_best = 2.43e-3;
    double dmsq_23_err = 0.06e-3;
    double sinsq_2th23_best = 4*sinsq_th23_best*(1-sinsq_th23_best);
    
    const double T2KBASELINE = 295;//km
    const double T2KEnergy = 0.6;//GeV
    
    
    double prob_measured = 0.1;//model2flav_vac_MuSurvive(T2KEnergy,T2KBASELINE,sinsq_2th23_best,dmsq_23_best);
    double proberr_measured = 0.05;//0.05*prob_measured;
    cout<<"prob measured "<<prob_measured<<" error "<<proberr_measured<<endl;
    double prob_uplimit = prob_measured+proberr_measured;
    double prob_lolimit = prob_measured-proberr_measured;
    
    
    
    const int NTHETA23BIN = 101;
    const double theta_MIN = 1.e-3;
    const double theta_MAX = 1.e0;
    
    const int NDM23BIN = 101;
    const double dm23_MIN = 1.e-5;
    const double dm23_MAX = 1.e+3;
    
    Double_t *ybins = new Double_t[NDM23BIN+1];
    double dybin = 8.0/NDM23BIN;
    double l10 = TMath::Log(10);
    for (int i=0;i<=NDM23BIN;i++) {
        ybins[i] = TMath::Exp(l10*(i-NDM23BIN*5/8.)*dybin);
    }
    
    Double_t *xbins = new Double_t[NTHETA23BIN+1];
    double dxbin = 3.0/NTHETA23BIN;
    for (int i=0; i<=NTHETA23BIN; i++) {
        xbins[i] = TMath::Exp(l10*(i-NTHETA23BIN)*dxbin);
    }
    
    
    TH2D *hsurface = new TH2D("hsurface","",NTHETA23BIN,xbins,NDM23BIN,ybins);
    
    const int NEVENTMC = 10000000;
    
    for (int idm=1; idm<=NDM23BIN; ++idm) {
        double dm_val = (ybins[idm]+ybins[idm-1])/2.;//dm23_MIN + (dm23_MAX-dm23_MIN)*(idm+0.5)/NDM23BIN;
        for (int itheta=1; itheta<=NTHETA23BIN; ++itheta) {
            double theta_val = (xbins[itheta]+xbins[itheta-1])/2.;//theta_MIN + (theta_MAX-theta_MIN)*(itheta+0.5)/NTHETA23BIN;
            double ProbMu2mu = model2flav_vac_MuSurvive(T2KEnergy,T2KBASELINE,theta_val,dm_val);
            hsurface->SetBinContent(itheta,idm,ProbMu2mu);
           
        }
    }
    
   
    
    
    
    TFile *popfile = new TFile("op_numu2numu_allowedregion.root","RECREATE");
    hsurface->Write();

    
    new TCanvas;
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gPad->SetLogx();
    hsurface->Draw("colz");
    hsurface->GetXaxis()->SetTitle("sin^{2}2#theta");
    hsurface->GetYaxis()->SetTitle("#Delta m^{2} (eV^{2}/c^{4})");
    gPad->Print("numu2nux_prob.eps");
    
    
    //narrow range
    new TCanvas;
    gPad->SetLogy();
    gPad->SetLogx();
    
    Int_t colors[]={2,4};
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    //allowed region
    Double_t levels[]={0.01,0.2};
    hsurface->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    hsurface->SetTitle("");
    hsurface->Draw("cont1");
    gPad->Update();
   //smooth
    if (isSmoothGraph) {
        Int_t nmaxloop=10;
        for (Int_t nloop=0; nloop<nmaxloop; ++nloop) {
            hsurface->Smooth(1,"k5b");
        }
    }
    gPad->Print(isSmoothGraph?"numu2numu_allowedregion_wcontour.eps":"numu2numu_allowedregion_wcontour_nosmooth.eps");
    
    //excluded region only
    Double_t levels2[]={0.01};
    hsurface->SetContour((sizeof(levels2)/sizeof(Double_t)), levels2);
    hsurface->SetTitle("");
    hsurface->Draw("cont1");
    gPad->Update();
    if (isSmoothGraph) {
    Int_t nmaxloop=10;
        for (Int_t nloop=0; nloop<nmaxloop; ++nloop) {
            hsurface->Smooth(1,"k5b");
        }
    }
    gPad->Print(isSmoothGraph?"numu2numu_excludedregion_wcontour.eps":"numu2numu_excludedregion_wcontour_nosmooth.eps");
    
    
   
}
