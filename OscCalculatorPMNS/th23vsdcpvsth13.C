{
    
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TString.h"
#include "TMath.h"
#include "OscCalculatorPMNS.h"
    gROOT->ProcessLine(".L OscCalculatorPMNS.C+");
    TString savename = "th23vsdcpvsth13";
    double x=0.0; //Energy
    Bool_t isTheta23Random = false;
    if(isTheta23Random) savename +="_th23rand";
    //Bool_t isNOVA=false;
    //Bool_t isLBNE=false;//T2K if false
    
    
    //input parameters
    TString nufitver="nufit5d2_2022nov";
    OscCalculatorPMNS myProb;
    double myDeltam31=2.511;
    double myDeltam21=7.41;//change from 7.54
    double myDeltam2=myDeltam31-myDeltam21*1e-2;//change from 2.42
    double myDeltam2Inv = -1.0*myDeltam31;
    double myTheta12=33.41*TMath::Pi()/180.;//0.60126;// this is radia   TMath::ASin(TMath::Sqrt(0.8704))/2
    double myTheta23=45.0*TMath::Pi()/180.;//49.0*TMath::Pi()/180.;//0.7854;// TMath::ASin(TMath::Sqrt(1))/2
   
    double mintheta23 = 39.6*TMath::Pi()/180.;
    double maxtheta23 = 51.8*TMath::Pi()/180.;
    
    //double myTheta13=8.57*TMath::Pi()/180.;//0.1479;//TMath::ASin(TMath::Sqrt(0.085))/2
    //T2K present at Neutrino 2022
    //root [0] TMath::ASin(TMath::Sqrt(0.0861))/2
    //(double) 0.14890537
    double myTheta13=8.54*TMath::Pi()/180.;//0.14890537;
    
    
    
    double myLT2K = 295.;
    double ent2k = 0.6;
    TString strEnergy = "ene0d6";
    /*double ent2k = 0.7;
     TString strEnergy = "ene0d7";*/
    /*double ent2k = 0.5;
    TString strEnergy = "ene0d5";*/
    
    double myLNOVA = 810.;
    double ennova = 2.0;
    
    double myLDUNE = 1250.;
    
    double myL=myLT2K;//isNOVA? 810:295.;
    double myEne = ent2k;
    
    double mydcp = -1.0*TMath::Pi()/2.;//0.0;
    
    
    //Set all parameters here
    myProb.SetL(myL);
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
    Double_t epsilon=0.00001;
    const Int_t nstepdcp=201;
    double dcpmin = -1.0*TMath::Pi();
    double dcpmax = 1.0*TMath::Pi();
    double BGdcpProj =0;
    
    
    const Int_t nstepth23 = 201;
    double sinsqth23min = 0.4;
    double sinsqth23max = 0.6;
    double BGsinsqth23Proj = 0.5;
    
    const Int_t nstepth13 = 201;
    double sinsqth13min = 0.019;
    double sinsqth13max = 0.025;//0.022+0.0006*5
    double BGsinsqth13Proj = 0.022;
    
    TH3D *hdcpvssinsqth23vssinsqth13_mu2mu = new TH3D("hdcpvssinsqth23vssinsqth13_mu2mu","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax,nstepth13,sinsqth13min,sinsqth13max);
    TH3D *hdcpvssinsqth23vssinsqth13_mu2e = new TH3D("hdcpvssinsqth23vssinsqth13_mu2e","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax,nstepth13,sinsqth13min,sinsqth13max);
    TH3D *hdcpvssinsqth23vssinsqth13_mu2mubar = new TH3D("hdcpvssinsqth23vssinsqth13_mu2mubar","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax,nstepth13,sinsqth13min,sinsqth13max);
    TH3D *hdcpvssinsqth23vssinsqth13_mu2ebar = new TH3D("hdcpvssinsqth23vssinsqth13_mu2ebar","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax,nstepth13,sinsqth13min,sinsqth13max);
    Int_t binBGdcpProj = hdcpvssinsqth23vssinsqth13_mu2mu->GetYaxis()->FindBin(BGdcpProj);
    Int_t binBGsinsqth23Proj = hdcpvssinsqth23vssinsqth13_mu2mu->GetXaxis()->FindBin(BGsinsqth23Proj);
    Int_t binBGsinsqth13Proj = hdcpvssinsqth23vssinsqth13_mu2mu->GetXaxis()->FindBin(BGsinsqth13Proj);
    
    TH2D *hdcpvssinsqth23_mu2mu = new TH2D("hdcpvssinsqth23_mu2mu","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth23_mu2e = new TH2D("hdcpvssinsqth23_mu2e","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth23_mu2mubar = new TH2D("hdcpvssinsqth23_mu2mubar","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth23_mu2ebar = new TH2D("hdcpvssinsqth23_mu2ebar","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax);
    
    TH2D *hdcpvssinsqth13_mu2mu = new TH2D("hdcpvssinsqth13_mu2mu","",nstepth13,sinsqth13min,sinsqth13max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth13_mu2e = new TH2D("hdcpvssinsqth13_mu2e","",nstepth13,sinsqth13min,sinsqth13max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth13_mu2mubar = new TH2D("hdcpvssinsqth13_mu2mubar","",nstepth13,sinsqth13min,sinsqth13max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth13_mu2ebar = new TH2D("hdcpvssinsqth13_mu2ebar","",nstepth13,sinsqth13min,sinsqth13max,nstepdcp,dcpmin,dcpmax);
    
    TH2D *hsinsqth13vssinsqth23_mu2mu = new TH2D("hsinsqth13vssinsqth23_mu2mu","",nstepth23,sinsqth23min,sinsqth23max,nstepth13,sinsqth13min,sinsqth13max);
    TH2D *hsinsqth13vssinsqth23_mu2e = new TH2D("hsinsqth13vssinsqth23_mu2e","",nstepth23,sinsqth23min,sinsqth23max,nstepth13,sinsqth13min,sinsqth13max);
    TH2D *hsinsqth13vssinsqth23_mu2mubar = new TH2D("hsinsqth13vssinsqth23_mu2mubar","",nstepth23,sinsqth23min,sinsqth23max,nstepth13,sinsqth13min,sinsqth13max);
    TH2D *hsinsqth13vssinsqth23_mu2ebar = new TH2D("hsinsqth13vssinsqth23_mu2ebar","",nstepth23,sinsqth23min,sinsqth23max,nstepth13,sinsqth13min,sinsqth13max);
    
    TH2D *hdcpvssinsqth23_mu2mu_ih = new TH2D("hdcpvssinsqth23_mu2mu_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth23_mu2e_ih = new TH2D("hdcpvssinsqth23_mu2e_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth23_mu2mubar_ih = new TH2D("hdcpvssinsqth23_mu2mubar_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth23_mu2ebar_ih = new TH2D("hdcpvssinsqth23_mu2ebar_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax);
    
    TH2D *hdcpvssinsqth13_mu2mu_ih = new TH2D("hdcpvssinsqth13_mu2mu_ih","",nstepth13,sinsqth13min,sinsqth13max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth13_mu2e_ih = new TH2D("hdcpvssinsqth13_mu2e_ih","",nstepth13,sinsqth13min,sinsqth13max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth13_mu2mubar_ih = new TH2D("hdcpvssinsqth13_mu2mubar_ih","",nstepth13,sinsqth13min,sinsqth13max,nstepdcp,dcpmin,dcpmax);
    TH2D *hdcpvssinsqth13_mu2ebar_ih = new TH2D("hdcpvssinsqth13_mu2ebar_ih","",nstepth13,sinsqth13min,sinsqth13max,nstepdcp,dcpmin,dcpmax);
    
    TH2D *hsinsqth13vssinsqth23_mu2mu_ih = new TH2D("hsinsqth13vssinsqth23_mu2mu_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepth13,sinsqth13min,sinsqth13max);
    TH2D *hsinsqth13vssinsqth23_mu2e_ih = new TH2D("hsinsqth13vssinsqth23_mu2e_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepth13,sinsqth13min,sinsqth13max);
    TH2D *hsinsqth13vssinsqth23_mu2mubar_ih = new TH2D("hsinsqth13vssinsqth23_mu2mubar_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepth13,sinsqth13min,sinsqth13max);
    TH2D *hsinsqth13vssinsqth23_mu2ebar_ih = new TH2D("hsinsqth13vssinsqth23_mu2ebar_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepth13,sinsqth13min,sinsqth13max);
    
    TH3D *hdcpvssinsqth23vssinsqth13_mu2mu_ih = new TH3D("hdcpvssinsqth23vssinsqth13_mu2mu_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax,nstepth13,sinsqth13min,sinsqth13max);
    TH3D *hdcpvssinsqth23vssinsqth13_mu2e_ih = new TH3D("hdcpvssinsqth23vssinsqth13_mu2e_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax,nstepth13,sinsqth13min,sinsqth13max);
    TH3D *hdcpvssinsqth23vssinsqth13_mu2mubar_ih = new TH3D("hdcpvssinsqth23vssinsqth13_mu2mubar_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax,nstepth13,sinsqth13min,sinsqth13max);
    TH3D *hdcpvssinsqth23vssinsqth13_mu2ebar_ih = new TH3D("hdcpvssinsqth23vssinsqth13_mu2ebar_ih","",nstepth23,sinsqth23min,sinsqth23max,nstepdcp,dcpmin,dcpmax,nstepth13,sinsqth13min,sinsqth13max);
    
    double probvaltmp;
    for (Int_t iz =0; iz<nstepth13; iz++) {
        double thsintheta13 = sinsqth13min+ (iz+0.5)*(sinsqth13max-sinsqth13min)/nstepth13;
        for (Int_t ix =0; ix<nstepth23; ix++) {
            double thsintheta23 = sinsqth23min+ (ix+0.5)*(sinsqth23max-sinsqth23min)/nstepth23;
            for (Int_t iy=0; iy<nstepdcp; iy++) {
                double thdcp = dcpmin + (iy+0.5)*(dcpmax-dcpmin)/nstepdcp;
                
                
                //for T2K, NH, dcp = 0
                myProb.Reset();
                myProb.SetL(myLT2K);
                myProb.SetRho(2.71);
                myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
                myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
                myProb.SetTh12(myTheta12);
                
                myTheta23 = TMath::ASin(TMath::Sqrt(thsintheta23));
                myProb.SetTh23(myTheta23);
                myTheta13 = TMath::ASin(TMath::Sqrt(thsintheta13));
                myProb.SetTh13(myTheta13);
                
                myProb.SetdCP(thdcp);
                myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
                myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
                
                isAnti = 1;
                probvaltmp = myProb.P(isAnti*14,isAnti*14,myEne);
                hdcpvssinsqth23vssinsqth13_mu2mu->SetBinContent(ix+1,iy+1,iz+1,probvaltmp);
                if(iz==binBGsinsqth13Proj)hdcpvssinsqth23_mu2mu->SetBinContent(ix+1,iy+1,probvaltmp);
                if(ix==binBGsinsqth23Proj)hdcpvssinsqth13_mu2mu->SetBinContent(iz+1,iy+1,probvaltmp);
                if(iy==binBGdcpProj)hsinsqth13vssinsqth23_mu2mu->SetBinContent(ix+1,iz+1,probvaltmp);
                
                probvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
                hdcpvssinsqth23vssinsqth13_mu2e->SetBinContent(ix+1,iy+1,iz+1,probvaltmp);
                if(iz==binBGsinsqth13Proj)hdcpvssinsqth23_mu2e->SetBinContent(ix+1,iy+1,probvaltmp);
                if(ix==binBGsinsqth23Proj)hdcpvssinsqth13_mu2e->SetBinContent(iz+1,iy+1,probvaltmp);
                if(iy==binBGdcpProj)hsinsqth13vssinsqth23_mu2e->SetBinContent(ix+1,iz+1,probvaltmp);
                
                isAnti = -1;
                probvaltmp = myProb.P(isAnti*14,isAnti*14,myEne);
                hdcpvssinsqth23vssinsqth13_mu2mubar->SetBinContent(ix+1,iy+1,iz+1,probvaltmp);
                if(iz==binBGsinsqth13Proj)hdcpvssinsqth23_mu2mubar->SetBinContent(ix+1,iy+1,probvaltmp);
                if(ix==binBGsinsqth23Proj)hdcpvssinsqth13_mu2mubar->SetBinContent(iz+1,iy+1,probvaltmp);
                if(iy==binBGdcpProj)hsinsqth13vssinsqth23_mu2mubar->SetBinContent(ix+1,iz+1,probvaltmp);
                
                probvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
                hdcpvssinsqth23vssinsqth13_mu2ebar->SetBinContent(ix+1,iy+1,iz+1,probvaltmp);
                if(iz==binBGsinsqth13Proj)hdcpvssinsqth23_mu2ebar->SetBinContent(ix+1,iy+1,probvaltmp);
                if(ix==binBGsinsqth23Proj)hdcpvssinsqth13_mu2ebar->SetBinContent(iz+1,iy+1,probvaltmp);
                if(iy==binBGdcpProj)hsinsqth13vssinsqth23_mu2ebar->SetBinContent(ix+1,iz+1,probvaltmp);
                //take inverted mass ordering
                myProb.Reset();
                myProb.SetL(myLT2K);
                myProb.SetRho(2.71);
                myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
                myProb.SetDmsq32(myDeltam2Inv*TMath::Power(10,-3.));//inverted
                myProb.SetTh12(myTheta12);
                
                myTheta23 = TMath::ASin(TMath::Sqrt(thsintheta23));
                myProb.SetTh23(myTheta23);
                myProb.SetTh13(myTheta13);
                
                myProb.SetdCP(thdcp);
                myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
                myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
                
                isAnti = 1;
                probvaltmp = myProb.P(isAnti*14,isAnti*14,myEne);
                hdcpvssinsqth23vssinsqth13_mu2mu_ih->SetBinContent(ix+1,iy+1,iz+1,probvaltmp);
                if(iz==binBGsinsqth13Proj)hdcpvssinsqth23_mu2mu_ih->SetBinContent(ix+1,iy+1,probvaltmp);
                if(ix==binBGsinsqth23Proj)hdcpvssinsqth13_mu2mu_ih->SetBinContent(iz+1,iy+1,probvaltmp);
                if(iy==binBGdcpProj)hsinsqth13vssinsqth23_mu2mu_ih->SetBinContent(ix+1,iz+1,probvaltmp);
                
                probvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
                hdcpvssinsqth23vssinsqth13_mu2e_ih->SetBinContent(ix+1,iy+1,iz+1,probvaltmp);
                if(iz==binBGsinsqth13Proj)hdcpvssinsqth23_mu2e_ih->SetBinContent(ix+1,iy+1,probvaltmp);
                if(ix==binBGsinsqth23Proj)hdcpvssinsqth13_mu2e_ih->SetBinContent(iz+1,iy+1,probvaltmp);
                if(iy==binBGdcpProj)hsinsqth13vssinsqth23_mu2e_ih->SetBinContent(ix+1,iz+1,probvaltmp);
                
                isAnti = -1;
                probvaltmp = myProb.P(isAnti*14,isAnti*14,myEne);
                hdcpvssinsqth23vssinsqth13_mu2mubar_ih->SetBinContent(ix+1,iy+1,iz+1,probvaltmp);
                if(iz==binBGsinsqth13Proj)hdcpvssinsqth23_mu2mubar_ih->SetBinContent(ix+1,iy+1,probvaltmp);
                if(ix==binBGsinsqth23Proj)hdcpvssinsqth13_mu2mubar_ih->SetBinContent(iz+1,iy+1,probvaltmp);
                if(iy==binBGdcpProj)hsinsqth13vssinsqth23_mu2mubar_ih->SetBinContent(ix+1,iz+1,probvaltmp);
                
                probvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
                hdcpvssinsqth23vssinsqth13_mu2ebar_ih->SetBinContent(ix+1,iy+1,iz+1,probvaltmp);
                if(iz==binBGsinsqth13Proj)hdcpvssinsqth23_mu2ebar_ih->SetBinContent(ix+1,iy+1,probvaltmp);
                if(ix==binBGsinsqth23Proj)hdcpvssinsqth13_mu2ebar_ih->SetBinContent(iz+1,iy+1,probvaltmp);
                if(iy==binBGdcpProj)hsinsqth13vssinsqth23_mu2ebar_ih->SetBinContent(ix+1,iz+1,probvaltmp);
                
            }
        }
    }
    
    TFile *poutfile = new TFile(Form("oscprob_%s_%s_%s.root",savename.Data(),nufitver.Data(),strEnergy.Data()),"RECREATE");
    hdcpvssinsqth23_mu2mu->Write();
    hdcpvssinsqth23_mu2e->Write();
    hdcpvssinsqth23_mu2mubar->Write();
    hdcpvssinsqth23_mu2ebar->Write();
    hdcpvssinsqth23_mu2mu_ih->Write();
    hdcpvssinsqth23_mu2e_ih->Write();
    hdcpvssinsqth23_mu2mubar_ih->Write();
    hdcpvssinsqth23_mu2ebar_ih->Write();
    
    hdcpvssinsqth13_mu2mu->Write();
    hdcpvssinsqth13_mu2e->Write();
    hdcpvssinsqth13_mu2mubar->Write();
    hdcpvssinsqth13_mu2ebar->Write();
    hdcpvssinsqth13_mu2mu_ih->Write();
    hdcpvssinsqth13_mu2e_ih->Write();
    hdcpvssinsqth13_mu2mubar_ih->Write();
    hdcpvssinsqth13_mu2ebar_ih->Write();
    
    hsinsqth13vssinsqth23_mu2mu->Write();
    hsinsqth13vssinsqth23_mu2e->Write();
    hsinsqth13vssinsqth23_mu2mubar->Write();
    hsinsqth13vssinsqth23_mu2ebar->Write();
    hsinsqth13vssinsqth23_mu2mu_ih->Write();
    hsinsqth13vssinsqth23_mu2e_ih->Write();
    hsinsqth13vssinsqth23_mu2mubar_ih->Write();
    hsinsqth13vssinsqth23_mu2ebar_ih->Write();
    
    hdcpvssinsqth23vssinsqth13_mu2mu->Write();
    hdcpvssinsqth23vssinsqth13_mu2e->Write();
    hdcpvssinsqth23vssinsqth13_mu2mubar->Write();
    hdcpvssinsqth23vssinsqth13_mu2ebar->Write();
    hdcpvssinsqth23vssinsqth13_mu2mu_ih->Write();
    hdcpvssinsqth23vssinsqth13_mu2e_ih->Write();
    hdcpvssinsqth23vssinsqth13_mu2mubar_ih->Write();
    hdcpvssinsqth23vssinsqth13_mu2ebar_ih->Write();
    
    poutfile->Close();
    
   
}
