{
    
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TString.h"
#include "OscCalculatorPMNS.h"
    gROOT->ProcessLine(".L OscCalculatorPMNS.C+");
    TString savename="numu2numu";
    double x=0.0; //Energy
    Bool_t isTheta23Random = false;
    if(isTheta23Random) savename +="_th23rand";
    
    
    //input parameters
    OscCalculatorPMNS myProb;
    double myDeltam31 = 2.514;
    double myDeltam21=7.42;//change from 7.54
    double myDeltam2=myDeltam31-myDeltam21*1e-2;//change from 2.42
    double myTheta12=33.44*TMath::Pi()/180.;//0.60126;// this is radia   TMath::ASin(TMath::Sqrt(0.8704))/2
    double myTheta23=49.0*TMath::Pi()/180.;//0.7854;// TMath::ASin(TMath::Sqrt(1))/2
    double myTheta13=8.57*TMath::Pi()/180.;//0.1479;//TMath::ASin(TMath::Sqrt(0.085))/2
    double myL_t2k=295;
    double myL_nova = 810;
    
    double mydcp = 0.0;
    
    
    
    //Set all parameters here
    myProb.SetL(myL_t2k);
    myProb.SetRho(2.71);
    myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
    myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
    
    myProb.SetTh12(myTheta12);
    myProb.SetTh23(myTheta23);
    myProb.SetTh13(myTheta13);
    myProb.SetdCP(0.0001);
    
    //Flag to trigger antineutrino probabilities +1=neutrinos, -1=antineutrinos
    int isAnti=1;
    
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
    
    
    TGraph *grProb_t2k_nh=0;
    TGraph *grProb_t2k_ih =0;
    TGraph *grProb_nova_nh=0;
    TGraph *grProb_nova_ih=0;
    
    TGraph *grProb_t2k_nh_max=0;
    TGraph *grProb_t2k_ih_max =0;
    TGraph *grProb_nova_nh_max=0;
    TGraph *grProb_nova_ih_max=0;
    
    TGraph *grProb_t2k_nh_min=0;
    TGraph *grProb_t2k_ih_min =0;
    TGraph *grProb_nova_nh_min=0;
    TGraph *grProb_nova_ih_min=0;
    
    const Int_t npoints=1000;
    Double_t loveremax =5;//energy range GeV
    
    Double_t musx[npoints];
    Double_t mu2mut2k_nh[npoints];
    Double_t mu2mut2k_ih[npoints];
    Double_t mu2munova_nh[npoints];
    Double_t mu2munova_ih[npoints];
    
    Double_t mu2mut2k_nh_max[npoints];
    Double_t mu2mut2k_ih_max[npoints];
    Double_t mu2munova_nh_max[npoints];
    Double_t mu2munova_ih_max[npoints];
    
    Double_t mu2mut2k_nh_min[npoints];
    Double_t mu2mut2k_ih_min[npoints];
    Double_t mu2munova_nh_min[npoints];
    Double_t mu2munova_ih_min[npoints];
    
    Double_t mu2mut2k_nh_rand;
    Double_t mu2mut2k_ih_rand;
    Double_t mu2munova_nh_rand;
    Double_t mu2munova_ih_rand;
    
    
    
    Double_t theta23_rand;
    const Int_t nstepth23 = 100;
    double mintheta23 = 39.6*TMath::Pi()/180.;
    double maxtheta23 = 51.8*TMath::Pi()/180.;
    
    //for dcp = 0
    for (Int_t i=1; i < npoints+1; i++){
        
        x = i/(npoints*1.0/loveremax) + epsilon;//vary energy range
        musx[i-1]=x;
        
        myProb.Reset();
        myProb.SetL(myL_t2k);
        myProb.SetRho(0.001);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        mu2mut2k_nh[i-1]=myProb.P(isAnti*14,isAnti*14,x);
        
        //for dcp = 0, IH
        myProb.Reset();
        myProb.SetL(myL_t2k);
        myProb.SetRho(0.001);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(-1.0*myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(3*2.0*TMath::Pi()/4 + epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        mu2mut2k_ih[i-1]=myProb.P(isAnti*14,isAnti*14,x);
        
        //nova + NH
        myProb.Reset();
        myProb.SetL(myL_nova);
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        mu2munova_nh[i-1]=myProb.P(isAnti*14,isAnti*14,x);
        
        //nova + NH
        myProb.Reset();
        myProb.SetL(myL_nova);
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(-myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        mu2munova_ih[i-1]=myProb.P(isAnti*14,isAnti*14,x);
        
        mu2mut2k_nh_max[i-1] = mu2munova_nh[i-1];
        mu2mut2k_ih_max[i-1] = mu2mut2k_ih[i-1];
        mu2munova_nh_max[i-1] = mu2munova_nh[i-1];
        mu2munova_ih_max[i-1] = mu2munova_ih[i-1];
        
        mu2mut2k_nh_min[i-1] = mu2munova_nh[i-1];
        mu2mut2k_ih_min[i-1] = mu2mut2k_ih[i-1];
        mu2munova_nh_min[i-1] = mu2munova_nh[i-1];
        mu2munova_ih_min[i-1] = mu2munova_ih[i-1];
        
   
        
        
        //loop over theta_23
        for (Int_t ith23=0; ith23 <= nstepth23; ith23++){
            theta23_rand = mintheta23 + (maxtheta23-mintheta23)*ith23/nstepth23;
            
            myProb.Reset();
            myProb.SetL(myL_t2k);
            myProb.SetRho(0.001);
            myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
            myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
            myProb.SetTh12(myTheta12);
            myProb.SetTh23(theta23_rand);
            myProb.SetTh13(myTheta13);
            myProb.SetdCP(epsilon);
            myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
            myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
            mu2mut2k_nh_rand=myProb.P(isAnti*14,isAnti*14,x);
            
            //for dcp = 0, IH
            myProb.Reset();
            myProb.SetL(myL_t2k);
            myProb.SetRho(0.001);
            myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
            myProb.SetDmsq32(-1.0*myDeltam2*TMath::Power(10,-3.));
            myProb.SetTh12(myTheta12);
            myProb.SetTh23(theta23_rand);
            myProb.SetTh13(myTheta13);
            myProb.SetdCP(3*2.0*TMath::Pi()/4 + epsilon);
            myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
            myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
            mu2mut2k_ih_rand=myProb.P(isAnti*14,isAnti*14,x);
            
            //nova + NH
            myProb.Reset();
            myProb.SetL(myL_nova);
            myProb.SetRho(2.71);
            myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
            myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
            myProb.SetTh12(myTheta12);
            myProb.SetTh23(theta23_rand);
            myProb.SetTh13(myTheta13);
            myProb.SetdCP(epsilon);
            myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
            myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
            mu2munova_nh_rand=myProb.P(isAnti*14,isAnti*14,x);
            
            //nova + NH
            myProb.Reset();
            myProb.SetL(myL_nova);
            myProb.SetRho(2.71);
            myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
            myProb.SetDmsq32(-myDeltam2*TMath::Power(10,-3.));
            myProb.SetTh12(myTheta12);
            myProb.SetTh23(theta23_rand);
            myProb.SetTh13(myTheta13);
            myProb.SetdCP(epsilon);
            myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
            myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
            mu2munova_ih_rand=myProb.P(isAnti*14,isAnti*14,x);
            
            //assign max and min
        
            
            if(mu2mut2k_nh_rand>mu2mut2k_nh_max[i-1]) mu2mut2k_nh_max[i-1] = mu2mut2k_nh_rand;
            if(mu2mut2k_ih_rand>mu2mut2k_ih_max[i-1]) mu2mut2k_ih_max[i-1] = mu2mut2k_ih_rand;
            if(mu2munova_nh_rand>mu2munova_nh_max[i-1]) mu2munova_nh_max[i-1] = mu2munova_nh_rand;
            if(mu2munova_ih_rand>mu2munova_ih_max[i-1]) mu2munova_ih_max[i-1] = mu2munova_ih_rand;
            
            if(mu2mut2k_nh_rand<mu2mut2k_nh_min[i-1]) mu2mut2k_nh_min[i-1] = mu2mut2k_nh_rand;
            if(mu2mut2k_ih_rand<mu2mut2k_ih_min[i-1]) mu2mut2k_ih_min[i-1] = mu2mut2k_ih_rand;
            if(mu2munova_nh_rand<mu2munova_nh_min[i-1]) mu2munova_nh_min[i-1] = mu2munova_nh_rand;
            if(mu2munova_ih_rand<mu2munova_ih_min[i-1]) mu2munova_ih_min[i-1] = mu2munova_ih_rand;
            
            
        }
    }
    
    
    
    
    grProb_t2k_nh = new TGraph(npoints,musx,mu2mut2k_nh);
    grProb_t2k_ih = new TGraph(npoints,musx,mu2mut2k_ih);
    grProb_nova_nh = new TGraph(npoints,musx,mu2munova_nh);
    grProb_nova_ih = new TGraph(npoints,musx,mu2munova_ih);
    
    grProb_t2k_nh_max = new TGraph(npoints,musx,mu2mut2k_nh_max);
    grProb_t2k_ih_max = new TGraph(npoints,musx,mu2mut2k_ih_max);
    grProb_nova_nh_max = new TGraph(npoints,musx,mu2munova_nh_max);
    grProb_nova_ih_max = new TGraph(npoints,musx,mu2munova_ih_max);
    
    grProb_t2k_nh_min = new TGraph(npoints,musx,mu2mut2k_nh_min);
    grProb_t2k_ih_min = new TGraph(npoints,musx,mu2mut2k_ih_min);
    grProb_nova_nh_min = new TGraph(npoints,musx,mu2munova_nh_min);
    grProb_nova_ih_min = new TGraph(npoints,musx,mu2munova_ih_min);
    
    //Plot everything
    
    TH2F *histo1=0;
    histo1 = new TH2F("NewOscProb","NewOscProb",npoints, 0., loveremax, 500, 0, 1.);
    
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
    
    
    TString LStr=Form("#delta_{CP}=%2.1f",epsilon);
    TLatex* tlx5=new TLatex(xcoord, 0.60+yshift,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    //tlx5->SetTextSize(0.03);
    tlx5->SetTextSize(22);
    tlx5->SetTextFont(43);
    tlx5->SetTextAlign(12);
    
    
    histo1->SetTitle("");
    
    histo1->GetXaxis()->SetTitle("Neutrino Energy (MeV)");
    histo1->GetYaxis()->SetTitle("Prob. (#nu_{#mu}#rightarrow #nu_{#mu}) ");
    histo1->GetYaxis()->CenterTitle();
    histo1->GetXaxis()->CenterTitle();
    histo1->GetXaxis()->SetLabelSize(histo1->GetXaxis()->GetTitleSize()*1.2);
    histo1->GetYaxis()->SetLabelSize(histo1->GetYaxis()->GetTitleSize()*1.2);
    histo1->GetXaxis()->SetTitleSize(histo1->GetXaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleSize(histo1->GetYaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    histo1->GetXaxis()->SetRangeUser(0.2,loveremax);
    histo1->Draw();
    grProb_t2k_nh->SetLineWidth(2);
    grProb_t2k_ih->SetLineWidth(2);
    grProb_nova_nh->SetLineWidth(2);
    grProb_nova_ih->SetLineWidth(2);
    
    Int_t ci;
    ci = TColor::GetColor("#0072B2");
    grProb_t2k_nh->SetLineColor(ci);
    grProb_t2k_ih->SetLineColor(ci);
    grProb_t2k_ih->SetLineStyle(3);
    grProb_t2k_nh->Draw("C same");
    grProb_t2k_ih->Draw("C same");
    
    if(isTheta23Random)
    {
        Double_t* exxcp2 =grProb_t2k_nh_max->GetX();
    Double_t* eyycp2 =grProb_t2k_nh_max->GetY();
    Int_t ngraphpointcp2=grProb_t2k_nh_max->GetN();
    double xmincp2 = TMath::MinElement(grProb_t2k_nh_max->GetN(), exxcp2);
    double xmaxcp2 = TMath::MaxElement(grProb_t2k_nh_max->GetN(), exxcp2);

    cout<<"xmin "<<xmincp2<<" xmax "<<xmaxcp2<<endl;
    TGraph *gabovecp2 = new TGraph(ngraphpointcp2+3);
    for (Int_t i=0;i<ngraphpointcp2;i++) gabovecp2->SetPoint(i,exxcp2[i],eyycp2[i]);
    gabovecp2->SetPoint(ngraphpointcp2,xmaxcp2,grProb_t2k_nh_max->GetMaximum());
    gabovecp2->SetPoint(ngraphpointcp2+1,xmincp2,grProb_t2k_nh_max->GetMaximum());
    gabovecp2->SetPoint(ngraphpointcp2+2,xmincp2,grProb_t2k_nh_max->GetMinimum());
    gabovecp2->SetFillColor(ci);
    gabovecp2->SetLineColor(kWhite);
    gabovecp2->Draw("f");
    
    Double_t* exxcpmin =grProb_t2k_nh_min->GetX();
    Double_t* eyycpmin =grProb_t2k_nh_min->GetY();
    Int_t ngraphpointcpmin=grProb_t2k_nh_min->GetN();
    double xmincpmin = TMath::MinElement(grProb_t2k_nh_min->GetN(), exxcpmin);
    double xmaxcpmin = TMath::MaxElement(grProb_t2k_nh_min->GetN(), exxcpmin);

    cout<<"xmin "<<xmincpmin<<" xmax "<<xmaxcpmin<<endl;
    TGraph *gabovecpmin = new TGraph(ngraphpointcpmin+3);
    for (Int_t i=0;i<ngraphpointcpmin;i++) gabovecpmin->SetPoint(i,exxcpmin[i],eyycpmin[i]);
    gabovecpmin->SetPoint(ngraphpointcpmin,xmaxcpmin,grProb_t2k_nh_min->GetMaximum());
    gabovecpmin->SetPoint(ngraphpointcpmin+1,xmincpmin,grProb_t2k_nh_min->GetMaximum());
    gabovecpmin->SetPoint(ngraphpointcpmin+2,xmincpmin,grProb_t2k_nh_min->GetMinimum());
    gabovecpmin->SetFillColor(kWhite);
    gabovecpmin->Draw("f");
    }
    
    
    ci = TColor::GetColor("#D55E00");
    grProb_nova_nh->SetLineColor(ci);
    grProb_nova_ih->SetLineColor(ci);
    grProb_nova_ih->SetLineStyle(3);
    
    grProb_nova_nh->Draw("C same");
    
    
    if (isTheta23Random){
    Double_t* exxcp2novamax =grProb_nova_nh_max->GetX();
    Double_t* eyycp2novamax =grProb_nova_nh_max->GetY();
    Int_t ngraphpointcp2novamax=grProb_nova_nh_max->GetN();
    double xmincp2novamax = TMath::MinElement(grProb_nova_nh_max->GetN(), exxcp2novamax);
    double xmaxcp2novamax = TMath::MaxElement(grProb_nova_nh_max->GetN(), exxcp2novamax);

    cout<<"xmin "<<xmincp2novamax<<" xmax "<<xmaxcp2novamax<<endl;
    TGraph *gabovecp2novamax = new TGraph(ngraphpointcp2novamax+3);
    for (Int_t i=0;i<ngraphpointcp2novamax;i++) gabovecp2novamax->SetPoint(i,exxcp2novamax[i],eyycp2novamax[i]);
    gabovecp2novamax->SetPoint(ngraphpointcp2novamax,xmaxcp2novamax,grProb_nova_nh_max->GetMaximum());
    gabovecp2novamax->SetPoint(ngraphpointcp2novamax+1,xmincp2novamax,grProb_nova_nh_max->GetMaximum());
    gabovecp2novamax->SetPoint(ngraphpointcp2novamax+2,xmincp2novamax,grProb_nova_nh_max->GetMinimum());
    gabovecp2novamax->SetFillColor(ci);
    gabovecp2novamax->SetLineColor(kWhite);
    gabovecp2novamax->Draw("f");
    
    Double_t* exxcp2novamin =grProb_nova_nh_min->GetX();
    Double_t* eyycp2novamin =grProb_nova_nh_min->GetY();
    Int_t ngraphpointcp2novamin=grProb_nova_nh_min->GetN();
    double xmincp2novamin = TMath::MinElement(grProb_nova_nh_min->GetN(), exxcp2novamin);
    double xmaxcp2novamin = TMath::MaxElement(grProb_nova_nh_min->GetN(), exxcp2novamin);

    cout<<"xmin "<<xmincp2novamin<<" xmax "<<xmaxcp2novamin<<endl;
    TGraph *gabovecp2novamin = new TGraph(ngraphpointcp2novamin+3);
    for (Int_t i=0;i<ngraphpointcp2novamin;i++) gabovecp2novamin->SetPoint(i,exxcp2novamin[i],eyycp2novamin[i]);
    gabovecp2novamin->SetPoint(ngraphpointcp2novamin,xmaxcp2novamin,grProb_nova_nh_min->GetMaximum());
    gabovecp2novamin->SetPoint(ngraphpointcp2novamin+1,xmincp2novamin,grProb_nova_nh_min->GetMaximum());
    gabovecp2novamin->SetPoint(ngraphpointcp2novamin+2,xmincp2novamin,grProb_nova_nh_min->GetMinimum());
    gabovecp2novamin->SetFillColor(kWhite);
    gabovecp2novamin->Draw("f");
    }
    

    
    
    

    grProb_nova_ih->Draw("C same");
    
    if (isTheta23Random){
    Double_t* exxcp2novaihmax =grProb_nova_ih_max->GetX();
    Double_t* eyycp2novaihmax =grProb_nova_ih_max->GetY();
    Int_t ngraphpointcp2novaihmax=grProb_nova_ih_max->GetN();
    double xmincp2novaihmax = TMath::MinElement(grProb_nova_ih_max->GetN(), exxcp2novaihmax);
    double xmaxcp2novaihmax = TMath::MaxElement(grProb_nova_ih_max->GetN(), exxcp2novaihmax);

    cout<<"xmin "<<xmincp2novaihmax<<" xmax "<<xmaxcp2novaihmax<<endl;
    TGraph *gabovecp2novaihmax = new TGraph(ngraphpointcp2novaihmax+3);
    for (Int_t i=0;i<ngraphpointcp2novaihmax;i++) gabovecp2novaihmax->SetPoint(i,exxcp2novaihmax[i],eyycp2novaihmax[i]);
    gabovecp2novaihmax->SetPoint(ngraphpointcp2novaihmax,xmaxcp2novaihmax,grProb_nova_ih_max->GetMaximum());
    gabovecp2novaihmax->SetPoint(ngraphpointcp2novaihmax+1,xmincp2novaihmax,grProb_nova_ih_max->GetMaximum());
    gabovecp2novaihmax->SetPoint(ngraphpointcp2novaihmax+2,xmincp2novaihmax,grProb_nova_ih_max->GetMinimum());
    gabovecp2novaihmax->SetFillColor(ci);
    gabovecp2novaihmax->SetLineColor(kWhite);
    gabovecp2novaihmax->Draw("f");
    
    Double_t* exxcp2novaihmin =grProb_nova_ih_min->GetX();
    Double_t* eyycp2novaihmin =grProb_nova_ih_min->GetY();
    Int_t ngraphpointcp2novaihmin=grProb_nova_ih_min->GetN();
    double xmincp2novaihmin = TMath::MinElement(grProb_nova_ih_min->GetN(), exxcp2novaihmin);
    double xmaxcp2novaihmin = TMath::MaxElement(grProb_nova_ih_min->GetN(), exxcp2novaihmin);

    cout<<"xmin "<<xmincp2novaihmin<<" xmax "<<xmaxcp2novaihmin<<endl;
    TGraph *gabovecp2novaihmin = new TGraph(ngraphpointcp2novaihmin+3);
    for (Int_t i=0;i<ngraphpointcp2novaihmin;i++) gabovecp2novaihmin->SetPoint(i,exxcp2novaihmin[i],eyycp2novaihmin[i]);
    gabovecp2novaihmin->SetPoint(ngraphpointcp2novaihmin,xmaxcp2novaihmin,grProb_nova_ih_min->GetMaximum());
    gabovecp2novaihmin->SetPoint(ngraphpointcp2novaihmin+1,xmincp2novaihmin,grProb_nova_ih_min->GetMaximum());
    gabovecp2novaihmin->SetPoint(ngraphpointcp2novaihmin+2,xmincp2novaihmin,grProb_nova_ih_min->GetMinimum());
    gabovecp2novaihmin->SetFillColor(kWhite);
    gabovecp2novaihmin->Draw("f");
    }
    
 
    

    
    tlx->Draw("same");
    tlx1->Draw("same");
    tlx2->Draw("same");
    tlx3->Draw("same");
    tlx4->Draw("same");
    tlx5->Draw("same");
    Double_t xlegmin = 0.35;
    Double_t ylegmin = 0.5;
    TLegend* leg = new TLegend(xlegmin,ylegmin,xlegmin+0.3,ylegmin+0.2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(22);
    leg->SetTextFont(43);
    leg->SetMargin(0.15);
    leg->AddEntry(grProb_t2k_nh, "T2K, Normal MO","l");
    leg->AddEntry(grProb_t2k_ih, "T2K, Inverted MO","l");
    leg->AddEntry(grProb_nova_nh, "NOvA, Normal MO","l");
    leg->AddEntry(grProb_nova_ih, "NOvA, Inverted MO","l");
    leg->Draw();
   canvas1->SaveAs(Form("plots/threeflav_exact_probability_matter_%s.eps",savename.Data()));
    
    canvas1->SaveAs(Form("plots/threeflav_exact_probability_matter_%s.pdf",savename.Data()));
    
    
    
    
    
}
