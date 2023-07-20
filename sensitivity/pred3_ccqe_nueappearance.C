//
//  pred2_ccqe_numudisapearance.C
// using ouput from sensitivity_ccqe_numudisapearance.C
////////////////////////////////////////////////////
//
//  To make illustration for
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu

void pred3_ccqe_nueappearance(){
    
    
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
    //to plot additional point
    int th13BIN_atBEST = (th13_BEST-th13_MIN)*(NTH13BIN-1)/(th13_MAX-th13_MIN)-0.5;
    int dcpBIN_atBEST = (dcp_BEST-dcp_MIN)*(NDCPBIN-1)/(dcp_MAX-dcp_MIN)-0.5;
    cout<<"BINS at best fit point for dcp "<<dcpBIN_atBEST<<" for th13 "<<th13BIN_atBEST<<endl;
    
    int ith_dcpbin = 1;
    int ith_th13bin = th13BIN_atBEST;
    double th13_val1 = th13_MIN + (th13_MAX-th13_MIN)*(ith_th13bin+0.5)/(NTH13BIN-1);
    double dcp_val1 = dcp_MIN + (dcp_MAX-dcp_MIN)*(ith_dcpbin+0.5)/(NDCPBIN-1);
    cout <<"prediction at th13 "<<th13_val1<<" dcp "<<dcp_val1<<endl;
    
    TH1D* hpred_osc1 = (TH1D*)hpred_osc[ith_th13bin][ith_dcpbin]->Clone("hpred_osc1");
    hpred_osc1->SetLineColor(2);
    hpred_osc1->SetLineWidth(2);
    hpred_osc1->Draw("hist same");
    
    ith_dcpbin = dcpBIN_atBEST;
    double th13_val2 = th13_MIN + (th13_MAX-th13_MIN)*(ith_th13bin+0.5)/(NTH13BIN-1);
    double dcp_val2 = dcp_MIN + (dcp_MAX-dcp_MIN)*(ith_dcpbin+0.5)/(NDCPBIN-1);
    
    cout <<"prediction at dm "<<th13_val2<<" theta "<<dcp_val2<<endl;
    TH1D* hpred_osc2 = (TH1D*)hpred_osc[ith_th13bin][ith_dcpbin]->Clone("hpred_osc2");
    hpred_osc2->SetLineColor(4);
    hpred_osc2->SetLineWidth(2);
    hpred_osc2->Draw("hist same");
    
    ith_th13bin = 49;
    double th13_val3 = th13_MIN + (th13_MAX-th13_MIN)*(ith_th13bin+0.5)/(NTH13BIN-1);
    double dcp_val3 = dcp_MIN + (dcp_MAX-dcp_MIN)*(ith_dcpbin+0.5)/(NDCPBIN-1);
    cout <<"prediction at dm "<<th13_val3<<" theta "<<dcp_val3<<endl;
    TH1D* hpred_osc3 = (TH1D*)hpred_osc[ith_th13bin][ith_dcpbin]->Clone("hpred_osc3");
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
    leg0->AddEntry(hpred_osc1, Form("(#theta_{13}=%.4f,#delta_{CP}=%.2f)",th13_val1,dcp_val1), "l");
    leg0->AddEntry(hpred_osc2, Form("(#theta_{13}=%.4f,#delta_{CP}=%.2f)",th13_val2,dcp_val2), "l");
    
    leg0->AddEntry(hpred_osc3, Form("(#theta_{13}=%.4f,#delta_{CP}=%.2f)",th13_val3,dcp_val3), "l");
    leg0->Draw();
    
  
    c1->Print("nuapp_pred3_datavsnoosc.eps");
    
    
    
    
  
    
  

    
    
}
