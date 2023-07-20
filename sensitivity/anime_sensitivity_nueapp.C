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
//input dcp_23 as sin^2 2theta 23 for convenience

void anime_sensitivity_nueapp(){
    
    
    TFile *pfile = new TFile("op_eventmc_4sensitivity_app.root","READ");
    const int NDCPBIN = 51;//50+1
    const double dcp_MIN = -TMath::Pi();
    const double dcp_MAX = TMath::Pi();
    const double dcp_BEST = -TMath::Pi()/2.;//45*Pi/180
    
    const int NTH13BIN = 51;//50+1
    const double th13_MIN = 0.020;
    const double th13_MAX = 0.024;
    const double th13_BEST = 0.022;
    
    
    
    const double EnergyRange2Draw = 1.5;//use for drawing
    
    TH1D *hpred_osc[NTH13BIN][NDCPBIN];
    for (int ith13=0; ith13<NTH13BIN; ++ith13) {
        for (int idcp=0; idcp<NDCPBIN; ++idcp) {
            hpred_osc[ith13][idcp] = (TH1D*)pfile->Get(Form("hpred_osc_th13bin%d_dcpbin%d",ith13,idcp));
        }
    }
    
    TH1D *hpred_noosc = (TH1D*)pfile->Get("hpred_noosc");
    TH1D *hdata = (TH1D*)pfile->Get("hdata");
    
    gBenchmark->Start("hsen");
    gSystem->Unlink("hsensitivity_anime_nueapp.gif"); // delete old file
    //condition to print
    const Int_t kUPDATE = 50;
    
    TCanvas *c1 = new TCanvas("c1","",0,0,1400,700);
    c1->Divide(2,1);
    c1->cd(1);
    gStyle->SetOptStat(0);
    hdata->SetMarkerStyle(8);
    hdata->SetMarkerColor(1);
    hdata->GetYaxis()->SetRangeUser(0,40);
    hdata->GetXaxis()->SetRangeUser(0,EnergyRange2Draw);
    hdata->GetYaxis()->SetTitle("Number of Events");
    hdata->GetXaxis()->SetTitle("Neutrino Energy [GeV]");

    hdata->Draw("ep");
    hdata->GetYaxis()->SetTitleOffset(1.5);
    hpred_noosc->SetLineWidth(2);
    hpred_noosc->SetLineColor(13);
    hpred_noosc->Draw("hist same");
    //to plot additional point
    int th13BIN_atBEST = (th13_BEST-th13_MIN)*NTH13BIN/(th13_MAX-th13_MIN)-0.5;
    int dcpBIN_atBEST = (dcp_BEST-dcp_MIN)*NDCPBIN/(dcp_MAX-dcp_MIN)-0.5;
    cout<<"BINS at best fit point for dcp "<<dcpBIN_atBEST<<" for th13 "<<th13BIN_atBEST<<endl;
    
    //emty histogram
    c1->cd(2);
     TGaxis::SetMaxDigits(3);
    TH2D *hsurface = new TH2D("hsurface","",NDCPBIN-1,dcp_MIN,dcp_MAX,NTH13BIN-1,th13_MIN,th13_MAX);
    hsurface->GetXaxis()->SetTitle("#delta_{CP}");
    hsurface->GetYaxis()->SetTitle("#theta_{13}");
    hsurface->GetYaxis()->SetTitleOffset(1.2);
    hsurface->Draw("AXIS");
    hsurface->GetZaxis()->SetRangeUser(0,25);
    Int_t nhist2plots=0;
    Int_t totalist = (NTH13BIN-1)* (NDCPBIN-1);
    for (int ith13=1; ith13<NTH13BIN; ++ith13) {
        for (int idcp=1; idcp<NDCPBIN; ++idcp) {
            ++nhist2plots;
            c1->cd(1);
            int colorindex = (ith13*idcp+1)%10;
            hpred_osc[ith13-1][idcp-1]->SetLineColor(colorindex);
            hpred_osc[ith13-1][idcp-1]->Draw("hist same");
            c1->cd(2);
            double likelihood = LogLikelihood(hpred_osc[ith13-1][idcp-1],hdata);
            hsurface->SetBinContent(idcp,ith13,likelihood);
            hsurface->Draw("colz");
            if(nhist2plots && (nhist2plots%kUPDATE==0)){
                c1->cd(1);
                hdata->Draw("ep same");
                c1->Update();
                if (gROOT->IsBatch()) {
                    c1->Print("hsensitivity_anime_nueapp.gif+");
                    printf("i = %d/%d\n", nhist2plots,totalist);
                }
            }
            
        }
    }
    
    
 
    
    if (gROOT->IsBatch()) c1->Print("hsensitivity_anime_nueapp.gif++");
    gBenchmark->Show("hsen");
    
    

    
  

    
    
}
