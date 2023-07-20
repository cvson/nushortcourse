//
//  pred_ccqe_numudisapearance.C
// using ouput from sensitivity_ccqe_numudisapearance.C
////////////////////////////////////////////////////
//
//  To make illustration for
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu

void pred_ccqe_nueappearance(){
    
    
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
    for (int idm=0; idm<NTH13BIN; ++idm) {
        for (int itheta=0; itheta<NDCPBIN; ++itheta) {
            hpred_osc[idm][itheta] = (TH1D*)pfile->Get(Form("hpred_osc_th13bin%d_dcpbin%d",idm,itheta));
        }
    }
    TH1D *hpred_noosc = (TH1D*)pfile->Get("hpred_noosc");
    TH1D *hdata = (TH1D*)pfile->Get("hdata");
    
    TCanvas *c1 = new TCanvas("c1","",0,0,800,600);
    
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
    TLegend* leg0 = new TLegend(.52, .58, 0.85, .82);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(18);
    leg0->SetTextFont(43);
    leg0->AddEntry(hdata, "Data", "ep");
    leg0->AddEntry(hpred_noosc, "Unoscillated Prediction", "l");
    leg0->Draw();
    
    /*c1->cd(2);
    gStyle->SetOptStat(0);
    TH1D* hratio = (TH1D*)hdata->Clone("hratio");
    hratio->Sumw2();
    hratio->Divide(hpred_noosc);
    for (Int_t ibin=1; ibin<=hratio->GetXaxis()->GetNbins(); ibin++) {
        Float_t err1;
        if (hpred_noosc->GetBinContent(ibin)!=0) {
            err1 = sqrt(hdata->GetBinContent(ibin))/hpred_noosc->GetBinContent(ibin);
            hratio->SetBinError(ibin,err1);
        }
      
    }
    hratio->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    hratio->GetYaxis()->SetTitle("Data / Prediction");
    
    hratio->GetYaxis()->SetRangeUser(0,2.);
    hratio->Draw("ep");
    hratio->GetYaxis()->SetTitleOffset(1.5);
    TLine *poneline = new TLine(0,1.,EnergyRange2Draw,1.0);
    poneline->SetLineStyle(7);
    poneline->Draw();
     */
    c1->Print("nueapp_pred_datavsnoosc.eps");
    
    
    
    

  

    
    
}
