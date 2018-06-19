{
    gROOT->ProcessLine(".x ../rootlogon.C");
    gROOT->ProcessLine(".L ../basicPlotUtil.C");
    
    int Nexperiments = 100000;
    int NtruthValue = 6;
    int NObservation = 10;
    double pvalue = 0;
    TRandom3 prandom;
    TH1F *hExp = new TH1F( "hExp", "hExp",  25, 0, 25);
    for (int i=0; i < Nexperiments; i++) {
        double nToy = prandom.Poisson(NtruthValue);
        hExp->Fill(nToy);
        if(nToy>=NObservation) pvalue += 1;
    }
    cout<<"P value is"<<pvalue/Nexperiments<<endl;
    
    new TCanvas;
    
    hExp->GetYaxis()->SetTitle("Number of toy experiments");
    hExp->GetXaxis()->SetTitle("Value of toy experiments");
    titleStyle(hExp);
    double ymax = hExp->GetMaximum();
    hExp->GetYaxis()->SetRangeUser(0, ymax*1.2);
    TGaxis::SetMaxDigits(3);
    hExp->Draw("hist");
    TLine *pline = new TLine(10,0,10,ymax*1.2);
    pline->SetLineWidth(2);
    pline->SetLineColor(2);
    pline->Draw();
    gPad->Print("plots_poison_test.eps");
    
}
