//
//  poison_test_nova.C
////////////////////////////////////////////////////
//
//  Simple calculation P value for NOvA nuebar appearance
//  result presented at Neutrino2018
//  run: root -b -q poison_test_nova.C
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
{
    gROOT->ProcessLine(".x ../rootlogon.C");
    gROOT->ProcessLine(".L ../basicPlotUtil.C");
    
    int Nexperiments = 100000000;//10000000
    int NBackground = 5.3;//4.92+0.55
    int NObservation = 18;
    double pvalue = 0;
    TRandom3 prandom;
    TH1F *hExp = new TH1F( "hExp", "hExp",  50, 0, 50);
    for (int i=0; i < Nexperiments; i++) {
        if (i%10000==0) cout<<"doing "<<i<<"/"<<Nexperiments<<endl;
        double nToy = prandom.Poisson(NBackground);
        hExp->Fill(nToy);
        if(nToy>=NObservation) pvalue += 1;
    }
    pvalue /= Nexperiments;
    cout<<"P value is"<<std::setprecision(6)<<pvalue<<" sigma "<<TMath::NormQuantile(1-pvalue)<<endl;
    
    //another check with build in function
    double pvalue_cdf = 1-ROOT::Math::poisson_cdf(NObservation-1,NBackground);
    double sigma_cdf = TMath::NormQuantile(1-pvalue_cdf);
    cout<<"Pvalue CDF "<<pvalue_cdf<<" sigma "<<sigma_cdf<<endl;
    
    //make plots
    new TCanvas;
    hExp->GetYaxis()->SetTitle("Number of toy experiments");
    hExp->GetXaxis()->SetTitle("Value of toy experiments");
    titleStyle(hExp);
    double ymax = hExp->GetMaximum();
    hExp->GetYaxis()->SetRangeUser(0.1, ymax*10);
    TGaxis::SetMaxDigits(3);
    gPad->SetLogy();
    hExp->Draw("hist");
    TLine *pline = new TLine(NObservation,0.1,NObservation,ymax*10);
    pline->SetLineWidth(2);
    pline->SetLineColor(2);
    pline->Draw();
    gPad->Print("plots_nova_discover_nuebar.eps");
    
}
