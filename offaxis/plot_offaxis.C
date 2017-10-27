void titleStyle(TH1* h1){
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.4);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.4);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleOffset(0.85);
}

void plot_offaxis(){
    TF1 *fnu_energy = new TF1("fnu_energy","0.427*x/(1+pow([0]*x,2)/pow(0.139570,2))",0.,20.);
    fnu_energy->SetParameter(0,0);
    new TCanvas;
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gPad->SetBottomMargin(gPad->GetBottomMargin()*1.4);
    gPad->SetLeftMargin(gPad->GetLeftMargin()*1.3);
    Int_t ci;
    TH1 *h1 = new TH1F("h1","h1",100,0,20);
    h1->GetXaxis()->SetTitle("True pion energy (GeV)");
    h1->GetYaxis()->SetTitle("True neutrino energy (GeV)");
    h1->GetYaxis()->SetRangeUser(0.,10.);
    h1->SetTitle("");
    titleStyle(h1);
    h1->DrawClone();
    Double_t xgr[4]={0,20,20,0};
    Double_t ygr[4]={2.5,2.5,1.5,1.5};
    TGraph *graph = new TGraph(4,xgr,ygr);
    graph->SetFillColor(16);
    graph->Draw("f");
     //h1->DrawClone("");
    fnu_energy->SetLineColor(kBlack);
    fnu_energy->SetTitle("");
    fnu_energy->SetLineWidth(3);
    fnu_energy->DrawClone("same");
    
    fnu_energy->SetParameter(0,7e-3);
    //UT orange
    ci = TColor::GetColor("#B45F04");
    fnu_energy->SetLineColor(ci);
    fnu_energy->DrawClone("same");
    
    fnu_energy->SetParameter(0,14e-3);

    //dodger blue
    ci = TColor::GetColor("#336699");
    fnu_energy->SetLineColor(ci);
    fnu_energy->DrawClone("same");
    
    fnu_energy->SetParameter(0,21e-3);
    //navy blue
    ci = TColor::GetColor("#669933");
    fnu_energy->SetLineColor(ci);
    fnu_energy->DrawClone("same");

    
    TString legForm="E_{#nu}=#frac{0.427E_{#pi}}{1+#gamma^{2}#theta^{2}}";
    TLatex* tlx=new TLatex(0.2, 0.75,legForm);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.055);
    tlx->SetTextFont(42);
    tlx->Draw();
    
    TString legOn="On axis";
    TLatex* tlx0=new TLatex(0.48, 0.45,legOn);
    tlx0->SetNDC(kTRUE); // <- use NDC coordinate
    tlx0->SetTextSize(0.055);
    tlx0->SetTextAngle(32);
    tlx0->Draw();
    
    TString leg7mrad="7 mrad off-axis";
    TLatex* tlx1=new TLatex(0.62, 0.46,leg7mrad);
    tlx1->SetNDC(kTRUE); // <- use NDC coordinate
    tlx1->SetTextSize(0.055);
    tlx1->SetTextAngle(8);
    //UT orange
    ci = TColor::GetColor("#B45F04");
    tlx1->SetTextColor(ci);
    tlx1->Draw();
    
    TString leg14mrad="14 mrad off-axis";
    TLatex* tlx2=new TLatex(0.55, 0.32,leg14mrad);
    tlx2->SetNDC(kTRUE); // <- use NDC coordinate
    tlx2->SetTextSize(0.055);
    tlx2->SetTextAngle(-3);
    //dodger blue
    ci = TColor::GetColor("#336699");
    tlx2->SetTextColor(ci);
    tlx2->Draw();
    
    TString leg21mrad="21 mrad off-axis";
    TLatex* tlx3=new TLatex(0.35, 0.2,leg21mrad);
    tlx3->SetNDC(kTRUE); // <- use NDC coordinate
    tlx3->SetTextSize(0.055);
    tlx3->SetTextAngle(-3);
    //navy blue
    ci = TColor::GetColor("#669933");
    tlx3->SetTextColor(ci);
    tlx3->Draw();
        gPad->Print("offaxis_energy.eps");

}