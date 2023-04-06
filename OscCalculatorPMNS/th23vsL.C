{
    
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TString.h"
#include "TMath.h"
#include "OscCalculatorPMNS.h"
    gROOT->ProcessLine(".L OscCalculatorPMNS.C+");
    TString savename = "th23vsL";
    double x=0.0; //Energy
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
    double myth23max = 45.0*TMath::Pi()/180.;
    
    //double myTheta13=8.57*TMath::Pi()/180.;//0.1479;//TMath::ASin(TMath::Sqrt(0.085))/2
    //T2K present at Neutrino 2022
    //root [0] TMath::ASin(TMath::Sqrt(0.0861))/2
    //(double) 0.14890537
    double myTheta13=8.54*TMath::Pi()/180.;//0.14890537;
    double fixdcp = 0.0001;
    
    
    
    double myLT2K = 295.;
    double ent2k = 0.6;
    TString strEnergy = "ene0d6";
    /*double ent2k = 0.7;
     TString strEnergy = "ene0d7";*/
    /*double ent2k = 0.5;
     TString strEnergy = "ene0d5";*/
    
    double myLNOVA = 810.;
    double ennova = 2.0;
    //TString strEnergy = "novaL810en2d0gev";
    
    double myLDUNE = 1285.; //1285km
    double endune = 2.5;//
    
    /*double myL=myLDUNE;//isNOVA? 810:295.;
     double myEne = endune;
     TString strEnergy = "duneL1285en2d5gev";*/
    
    double myLESSnuSB = 360;
    double enESSnuSB = 0.23;
    //TString strEnergy = "essnusb_ene0d23";
    
    double myL=myLT2K;//isNOVA? 810:295.;
    double myEne = ent2k;
    
    
    
    
    //Set all parameters here
    myProb.SetL(myL);
    myProb.SetRho(2.71);
    myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
    myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
    myProb.SetTh12(myTheta12);
    myProb.SetTh23(myTheta23);
    myProb.SetTh13(myTheta13);
    myProb.SetdCP(fixdcp);
    
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
    const Int_t nstepdcp=200;
    double dcpmin = -1.0*TMath::Pi();
    double dcpmax = 1.0*TMath::Pi();
    double dcpmin4hist = dcpmin-0.5*(dcpmax-dcpmin)/nstepdcp;
    double dcpmax4hist = dcpmax-0.5*(dcpmax-dcpmin)/nstepdcp;
    double BGdcpProj =0; //fix value of dcp in 2D th23 vs th13
    
    
    const Int_t nstepth23 = 200;
    
    double sinsqth23min = 0.45;//0.4
    double sinsqth23max = 0.55;//0.6
    double sinsqth23min4hist = sinsqth23min-0.5*(sinsqth23max-sinsqth23min)/nstepth23;
    double sinsqth23max4hist = sinsqth23max-0.5*(sinsqth23max-sinsqth23min)/nstepth23;
    double BGsinsqth23Proj = 0.5;//fix value of th23 in 2d dcp vs th13
    
    const Int_t nstepth13 = 200;
    double sinsqth13min = 0.019;
    double sinsqth13max = 0.025;//0.022+0.0006*5
    double sinsqth13min4hist = sinsqth13min - 0.5*(sinsqth13max-sinsqth13min)/nstepth13;
    double sinsqth13max4hist = sinsqth13max - 0.5*(sinsqth13max-sinsqth13min)/nstepth13;
    double BGsinsqth13Proj = 0.022;//fix value of th13 in 2D dcp vs th23
    
    const Int_t nstepL = 250;
    double Lmin = 10;//km
    double Lmax = 510;//km
    double Lmin4hist = Lmin - 0.5*(Lmax-Lmin)/nstepL;
    double Lmax4hist = Lmax - 0.5*(Lmax-Lmin)/nstepL;
    
    
    
    TH2D *hth23vsL_mu2e = new TH2D("hth23vsL_mu2e","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    Int_t BinIndexTh23Max = hth23vsL_mu2e->GetYaxis()->FindBin(0.5);
    
    TH2D *hth23vsL_mu2ebar = new TH2D("hth23vsL_mu2ebar","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    
    TH2D *hth23vsL_mu2mu = new TH2D("hth23vsL_mu2mu","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    TH2D *hth23vsL_mu2mubar = new TH2D("hth23vsL_mu2mubar","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    
    //chisq2 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.053015
    TH2D *hth23vsL_chisqs1nue = new TH2D("hth23vsL_chisqs1nue","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    TH2D *hth23vsL_chisqs1nuebar = new TH2D("hth23vsL_chisqs1nuebar","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    TH2D *hth23vsL_chisqs1nuesum = new TH2D("hth23vsL_chisqs1nuesum","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    
    TH2D *hth23vsL_chisqs1ORnue = new TH2D("hth23vsL_chisqs1ORnue","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    TH2D *hth23vsL_chisqs1ORnuebar = new TH2D("hth23vsL_chisqs1ORnuebar","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    TH2D *hth23vsL_chisqs1ORnuesum = new TH2D("hth23vsL_chisqs1ORnuesum","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    
    
    TH2D *hth23vsL_chisqs1numu = new TH2D("hth23vsL_chisqs1numu","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    TH2D *hth23vsL_chisqs1numubar = new TH2D("hth23vsL_chisqs1numubar","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    TH2D *hth23vsL_chisqs1numusum = new TH2D("hth23vsL_chisqs1numusum","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
   
    TH2D *hth23vsL_mu2e_ih = new TH2D("hth23vsL_mu2e_ih","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    TH2D *hth23vsL_mu2ebar_ih = new TH2D("hth23vsL_mu2ebar_ih","",nstepL,Lmin4hist,Lmax4hist,nstepth23,sinsqth23min4hist,sinsqth23max4hist);
    
    
    
    
    
    double probvaltmp;
    double probvaltmpAnti;
    
    double probvalnumutmp;
    double probvalnumutmpAnti;
   
    double th23maxprobvaltmp;
    double th23maxprobvaltmpAnti;
    
    double th23maxprobvalnumutmp;
    double th23maxprobvalnumutmpAnti;
    
    double tmpchisq1ORsubnu;
    double tmpchisq1ORsubnubar;
    
    double tmpchisq1subnu;
    double tmpchisq1subnubar;
    
    double tmpchisq1subnumu;
    double tmpchisq1subnumubar;
    


    
    double tmpMinchisq;
    
    for (Int_t ix =0; ix<nstepL; ix++) {
        double thL = Lmin4hist + (ix+0.5)*(Lmax-Lmin)/nstepL;
        //get probability for 0 and pi
        //for T2K, normal ordering
        myProb.Reset();
        myProb.SetL(thL);
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myth23max);
        myProb.SetTh13(myTheta13);
        
        myProb.SetdCP(fixdcp);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        
        isAnti = 1;
        th23maxprobvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
        th23maxprobvalnumutmp = myProb.P(isAnti*14,isAnti*14,myEne);
        
        isAnti = -1;
        th23maxprobvaltmpAnti = myProb.P(isAnti*14,isAnti*12,myEne);
        th23maxprobvalnumutmpAnti = myProb.P(isAnti*14,isAnti*14,myEne);
        
        for (Int_t iy=0; iy<nstepth23; iy++) {
            double tmpsinsqth23 = sinsqth23min4hist + (iy+0.5)*(sinsqth23max-sinsqth23min)/nstepth23;
            myTheta23 = TMath::ASin(TMath::Sqrt(tmpsinsqth23));
            //for T2K, normal ordering
            myProb.Reset();
            myProb.SetL(thL);
            myProb.SetRho(2.71);
            myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
            myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
            myProb.SetTh12(myTheta12);
            myProb.SetTh23(myTheta23);
            myProb.SetTh13(myTheta13);
            
            myProb.SetdCP(fixdcp);
            myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
            myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
            
            isAnti = 1;
            probvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
            hth23vsL_mu2e->SetBinContent(ix+1,iy+1,probvaltmp);
            
            probvalnumutmp = myProb.P(isAnti*14,isAnti*14,myEne);
            hth23vsL_mu2mu->SetBinContent(ix+1,iy+1,probvalnumutmp);
          
            
            isAnti = -1;
            probvaltmpAnti = myProb.P(isAnti*14,isAnti*12,myEne);
            hth23vsL_mu2ebar->SetBinContent(ix+1,iy+1,probvaltmpAnti);
            
            probvalnumutmpAnti = myProb.P(isAnti*14,isAnti*14,myEne);
            hth23vsL_mu2mubar->SetBinContent(ix+1,iy+1,probvalnumutmpAnti);
            
           
            
            
            tmpchisq1subnu =(1./pow(thL,2.))*(pow(probvaltmp-th23maxprobvaltmp,2)/probvaltmp);
            hth23vsL_chisqs1nue->SetBinContent(ix+1,iy+1,tmpchisq1subnu);
            
            tmpchisq1subnubar = (1./pow(thL,2.))*(pow(probvaltmpAnti-th23maxprobvaltmpAnti,2)/probvaltmpAnti);
            hth23vsL_chisqs1nuebar->SetBinContent(ix+1,iy+1,tmpchisq1subnubar);
            
            hth23vsL_chisqs1nuesum->SetBinContent(ix+1,iy+1,tmpchisq1subnu+tmpchisq1subnubar);
            
            tmpchisq1subnumu =(1./pow(thL,2.))*(pow(probvalnumutmp-th23maxprobvalnumutmp,2)/probvalnumutmp);
            hth23vsL_chisqs1numu->SetBinContent(ix+1,iy+1,tmpchisq1subnumu);
            
            tmpchisq1subnumubar = (1./pow(thL,2.))*(pow(probvalnumutmpAnti-th23maxprobvalnumutmpAnti,2)/probvalnumutmpAnti);
            hth23vsL_chisqs1numubar->SetBinContent(ix+1,iy+1,tmpchisq1subnumubar);
            
            hth23vsL_chisqs1nuesum->SetBinContent(ix+1,iy+1,tmpchisq1subnu+tmpchisq1subnubar);
            hth23vsL_chisqs1numusum->SetBinContent(ix+1,iy+1,tmpchisq1subnumu+tmpchisq1subnumubar);
            //octant resolving
            
            
            //take inverted mass ordering
            myProb.Reset();
            myProb.SetL(thL);
            myProb.SetRho(2.71);
            myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
            myProb.SetDmsq32(myDeltam2Inv*TMath::Power(10,-3.));//inverted
            myProb.SetTh12(myTheta12);
            myProb.SetTh23(myTheta23);
            myProb.SetTh13(myTheta13);
            
            myProb.SetdCP(fixdcp);
            myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
            myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
            
            isAnti = 1;
            probvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
            hth23vsL_mu2e_ih->SetBinContent(ix+1,iy+1,probvaltmp);
                       
            isAnti = -1;
            probvaltmpAnti = myProb.P(isAnti*14,isAnti*12,myEne);
            hth23vsL_mu2ebar_ih->SetBinContent(ix+1,iy+1,probvaltmpAnti);
            
            //octant resolving
            /*if((iy+1)==BinIndexTh23Max){
                hth23vsL_chisqs1ORnue->SetBinContent(ix+1, iy+1,0.0);
                hth23vsL_chisqs1ORnuebar->SetBinContent(ix+1, iy+1,0.0);
                hth23vsL_chisqs1ORnuesum->SetBinContent(ix+1, iy+1,0.0);
            }
            else if ((iy+1)>BinIndexTh23Max){
                
                
            }
            else {
                
            }*/
            
        }
    }
    
    TFile *poutfile = new TFile(Form("oscprob_%s_%s_%s.root",savename.Data(),nufitver.Data(),strEnergy.Data()),"RECREATE");
    hth23vsL_mu2e->Write();
    hth23vsL_mu2ebar->Write();
    hth23vsL_mu2mu->Write();
    hth23vsL_mu2mubar->Write();
    hth23vsL_mu2e_ih->Write();
    hth23vsL_mu2ebar_ih->Write();
   
   
    hth23vsL_chisqs1nue->Write();
    hth23vsL_chisqs1nuebar->Write();
    hth23vsL_chisqs1nuesum->Write();
    
    hth23vsL_chisqs1numu->Write();
    hth23vsL_chisqs1numubar->Write();
    hth23vsL_chisqs1numusum->Write();
    
    poutfile->Close();
    
    
}
