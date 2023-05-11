{
    
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TString.h"
#include "OscCalculatorPMNS.h"
    gROOT->ProcessLine(".L OscCalculatorPMNS.C+");
    
    double x=0.0; //Energy
    
    
    
    //input parameters
    OscCalculatorPMNS myProb;
    double myDeltam31 = 2.514;
    double myDeltam21=7.42;//change from 7.54
    double myDeltam2=myDeltam31-myDeltam21*1e-2;//change from 2.42
    double myTheta12=33.44*TMath::Pi()/180.;//0.60126;// this is radia   TMath::ASin(TMath::Sqrt(0.8704))/2
    double myTheta23=49.0*TMath::Pi()/180.;//0.7854;// TMath::ASin(TMath::Sqrt(1))/2
    double myTheta13=8.57*TMath::Pi()/180.;//0.1479;//TMath::ASin(TMath::Sqrt(0.085))/2
    double myL=1.6;
    
    double mydcp = 0.0;
    
    
    
    //Set all parameters here
    myProb.SetL(myL);
    
    myProb.SetRho(0.001);
    myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
    myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
    
    myProb.SetTh12(myTheta12);
    myProb.SetTh23(myTheta23);
    myProb.SetTh13(myTheta13);
    myProb.SetdCP(0.0001);
    
    //Flag to trigger antineutrino probabilities +1=neutrinos, -1=antineutrinos
    int isAnti=-1;
    
    //Note, probability flavors are defined by pdg code:
    //nue  =12, antinue  =-12
    //numu =14, antinue  =-14
    //nutau=16, antinutau=-16
    //For instance, to compute P(numu->nue) you would call:
    // myProb.P(14,12,E)
    
    //Read back parameters to check
    std::cout << myProb.GetL() << endl
    << myProb.GetRho() << endl
    << myProb.GetDmsq21() << endl
    << myProb.GetDmsq32() << endl
    << myProb.GetTh12() << endl
    << myProb.GetTh23() << endl
    << myProb.GetTh13() << endl
    << myProb.GetdCP() << endl
    << endl;
    
    //This step is done automatically when computing probabilities,
    //but useful to do it here to check code works
    myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
    myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
    
    //Print mixing matrix and delta_m_sqrs.
    myProb.PrintMix();
    myProb.PrintDeltaMsqrs();
    
    //Keep this here for reference
    /*
     double par[9]={0.61,                       //Theta12
     myTheta23,                  //Theta23
     0.086,                       //Theta13
     myDeltam2*TMath::Power(10,-3.), //DeltaM32^2
     7.59*TMath::Power(10,-5.),  //DeltaM21^2
     0.0001,                       //DeltaCP (3Pi/2)
     2.71,                        //Rock Density
     myL,                        //L
     1};                        //NuAntiNu
     */
    double epsilon=0.00001;
    
    
    TGraph *e2eDcp0NH=0;
    TGraph *e2eDcp0IH =0;
    TGraph *Antie2eDcp0dm21only=0;
    
    
    const Int_t npoints=1000;
    Double_t loveremax =0.01;//energy range GeV
    
    Double_t musx[npoints];
    Double_t musynueDcp0[npoints];
    Double_t musynueDcp270[npoints];
    Double_t antimusynueDcp0[npoints];
    
    //for dcp = 0
    for (Int_t i=1; i < npoints+1; i++){
        
        x = i/(npoints*1.0/loveremax) + epsilon;//vary energy range
        musx[i-1]=x*1000;
        
        myProb.Reset();
        myProb.SetL(myL);
        myProb.SetRho(0.001);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        musynueDcp0[i-1]=myProb.P(isAnti*12,isAnti*12,x);
        
    }
    
    //for dcp = 0, IH
    for (Int_t i=1; i < npoints+1; i++){
        
        x = i/(npoints*1.0/loveremax) + epsilon;//vary energy range
        
        myProb.Reset();
        myProb.SetL(myL);
        myProb.SetRho(0.001);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(-1.0*myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(3*2.0*TMath::Pi()/4 + epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        musynueDcp270[i-1]=myProb.P(isAnti*12,isAnti*12,x);
        
    }
    
    //for dcp=0, dm32 = 0
    
    for (Int_t i=1; i < npoints+1; i++){
        
        x = i/(npoints*1.0/loveremax) + epsilon;//vary energy range
        
        myProb.Reset();
        myProb.SetL(myL);
        myProb.SetRho(0.001);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(0.0*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        antimusynueDcp0[i-1]=myProb.P(isAnti*12,isAnti*12,x);
        
    }
    
   
    
    
    
    
    e2eDcp0NH = new TGraph(npoints,musx,musynueDcp0);
    e2eDcp0IH = new TGraph(npoints,musx,musynueDcp270);
    Antie2eDcp0dm21only = new TGraph(npoints,musx,antimusynueDcp0);
    
    
    
    //Plot everything
    
    TH2F *histo1=0;
    histo1 = new TH2F("NewOscProb","NewOscProb",npoints, 0., loveremax*1000, 500, 0, 1.);
    
    TCanvas *canvas1=0;
    canvas1 = new TCanvas("Can1","Can1",1000,800);
    
    canvas1->SetLeftMargin(canvas1->GetLeftMargin()*1.2);
    canvas1->SetBottomMargin(canvas1->GetBottomMargin()*1.2);
    
    
    canvas1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    //gPad->SetLeftMargin(gPad->GetLeftMargin()*1.4);
    double xcoord = 0.62;
    double yshift = -0.4;
    TString delMStr32=Form("%2.2f #times 10^{-3} eV^{2}",myDeltam2);
    TString legStr32="|#Deltam^{2}_{32}|=";
    TLatex* tlx=new TLatex(xcoord, 0.85+yshift,legStr32+delMStr32);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    //tlx->SetTextSize(0.03);
    tlx->SetTextSize(22);
    tlx->SetTextFont(43);
    tlx->SetTextAlign(12);
    
    TString delMStr21=Form("%2.2f #times 10^{-5} eV^{2}",myDeltam21);
    TString legStr21="#Deltam^{2}_{21}=";
    TLatex* tlx1=new TLatex(xcoord, 0.80+yshift,legStr21+delMStr21);
    tlx1->SetNDC(kTRUE); // <- use NDC coordinate
    //tlx1->SetTextSize(0.03);
    tlx1->SetTextSize(22);
    tlx1->SetTextFont(43);
    tlx1->SetTextAlign(12);
    
    
    TString Th23Str=Form("sin^{2}2#theta_{23}=%2.2f",TMath::Sin(2*myTheta23)*TMath::Sin(2*myTheta23));
    TLatex* tlx2=new TLatex(xcoord, 0.75+yshift,Th23Str);
    tlx2->SetNDC(kTRUE); // <- use NDC coordinate
    //tlx2->SetTextSize(0.03);
    tlx2->SetTextSize(22);
    tlx2->SetTextFont(43);
    tlx2->SetTextAlign(12);
    
    TString Th12Str=Form("sin^{2}2#theta_{12}=%2.2f",TMath::Sin(2*myTheta12)*TMath::Sin(2*myTheta12));
    TLatex* tlx3=new TLatex(xcoord, 0.70+yshift,Th12Str);
    tlx3->SetNDC(kTRUE); // <- use NDC coordinate
    //tlx3->SetTextSize(0.03);
    tlx3->SetTextSize(22);
    tlx3->SetTextFont(43);
    tlx3->SetTextAlign(12);
    
    TString Th13Str=Form("sin^{2}2#theta_{13}=%2.2f",TMath::Sin(2*myTheta13)*TMath::Sin(2*myTheta13));
    TLatex* tlx4=new TLatex(xcoord, 0.65+yshift,Th13Str);
    tlx4->SetNDC(kTRUE); // <- use NDC coordinate
    //tlx4->SetTextSize(0.03);
    tlx4->SetTextSize(22);
    tlx4->SetTextFont(43);
    tlx4->SetTextAlign(12);
    
    
    TString LStr=Form("L=%2.1f km",myL);
    TLatex* tlx5=new TLatex(xcoord, 0.60+yshift,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    //tlx5->SetTextSize(0.03);
    tlx5->SetTextSize(22);
    tlx5->SetTextFont(43);
    tlx5->SetTextAlign(12);
    
    
    histo1->SetTitle("");
    
    histo1->GetXaxis()->SetTitle("Neutrino Energy (MeV)");
    histo1->GetYaxis()->SetTitle("Prob. (#bar{#nu}_{e}#rightarrow#bar{#nu}_{e}) ");
    histo1->GetYaxis()->CenterTitle();
    histo1->GetXaxis()->CenterTitle();
    histo1->GetXaxis()->SetLabelSize(histo1->GetXaxis()->GetTitleSize()*1.2);
    histo1->GetYaxis()->SetLabelSize(histo1->GetYaxis()->GetTitleSize()*1.2);
    histo1->GetXaxis()->SetTitleSize(histo1->GetXaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleSize(histo1->GetYaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    histo1->GetXaxis()->SetRangeUser(0.2,loveremax);
    histo1->Draw();
    e2eDcp0NH->SetLineWidth(2);
    e2eDcp0IH->SetLineWidth(2);
    Antie2eDcp0dm21only->SetLineWidth(2);
    
    Int_t ci;
    
    ci = TColor::GetColor("#0072B2");
    e2eDcp0NH->SetLineColor(ci);
    e2eDcp0NH->Draw("C same");
    
    ci = TColor::GetColor("#D55E00");
    e2eDcp0IH->SetLineColor(ci);
    e2eDcp0IH->Draw("C same");
    
    
    ci = TColor::GetColor("#000000");
    Antie2eDcp0dm21only->SetLineColor(ci);
    Antie2eDcp0dm21only->SetLineStyle(3);
    Antie2eDcp0dm21only->Draw("C same");
    

    
    tlx->Draw("same");
    tlx1->Draw("same");
    tlx2->Draw("same");
    tlx3->Draw("same");
    tlx4->Draw("same");
    tlx5->Draw("same");
    Double_t xlegmin = 0.35;
    Double_t ylegmin = 0.5;
    TLegend* leg = new TLegend(xlegmin,ylegmin,xlegmin+0.3,ylegmin+0.23);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(22);
    leg->SetTextFont(43);
    leg->AddEntry(e2eDcp0NH, "Normal Mass Ordering","l");
    leg->AddEntry(e2eDcp0IH, "Inverted Mass Ordering","l");
    leg->AddEntry(Antie2eDcp0dm21only, "(#theta_{12},#Deltam^{2}_{21}) oscillation only","l");
    leg->Draw();
 canvas1->SaveAs("plots/dayabay_threeflav_exact_probability_matter_nue2nue.eps");
    canvas1->SaveAs("plots/daybay_threeflav_exact_probability_matter_nue2nue.pdf");
    
    
    
    
    
    
    
}
