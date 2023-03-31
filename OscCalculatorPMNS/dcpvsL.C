{
    
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TString.h"
#include "TMath.h"
#include "OscCalculatorPMNS.h"
    gROOT->ProcessLine(".L OscCalculatorPMNS.C+");
    TString savename = "dcpvsL";
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
    //TString strEnergy = "novaL810en2d0gev";
    
    double myLDUNE = 1285.; //1285km
    double endune = 2.5;//
    
    /*double myL=myLDUNE;//isNOVA? 810:295.;
     double myEne = endune;
     TString strEnergy = "duneL1285en2d5gev";*/
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
    const Int_t nstepdcp=200;
    double dcpmin = -1.0*TMath::Pi();
    double dcpmax = 1.0*TMath::Pi();
    double dcpmin4hist = dcpmin-0.5*(dcpmax-dcpmin)/nstepdcp;
    double dcpmax4hist = dcpmax-0.5*(dcpmax-dcpmin)/nstepdcp;
    double BGdcpProj =0; //fix value of dcp in 2D th23 vs th13
    
    
    const Int_t nstepth23 = 200;
    
    double sinsqth23min = 0.4;
    double sinsqth23max = 0.6;
    double sinsqth23min4hist = sinsqth23min-0.5*(sinsqth23max-sinsqth23min)/nstepth23;
    double sinsqth23max4hist = sinsqth23max-0.5*(sinsqth23max-sinsqth23min)/nstepth23;
    double BGsinsqth23Proj = 0.5;//fix value of th23 in 2d dcp vs th13
    
    const Int_t nstepth13 = 200;
    double sinsqth13min = 0.019;
    double sinsqth13max = 0.025;//0.022+0.0006*5
    double sinsqth13min4hist = sinsqth13min - 0.5*(sinsqth13max-sinsqth13min)/nstepth13;
    double sinsqth13max4hist = sinsqth13max - 0.5*(sinsqth13max-sinsqth13min)/nstepth13;
    double BGsinsqth13Proj = 0.022;//fix value of th13 in 2D dcp vs th23
    
    const Int_t nstepL = 200;
    double Lmin = 10;//km
    double Lmax = 410;//km
    double Lmin4hist = Lmin - 0.5*(Lmax-Lmin)/nstepL;
    double Lmax4hist = Lmax - 0.5*(Lmax-Lmin)/nstepL;
    
    
    
    TH2D *hdcpvsL_mu2e = new TH2D("hdcpvsL_mu2e","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    TH2D *hdcpvsL_mu2ebar = new TH2D("hdcpvsL_mu2ebar","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    
    TH2D *hdcpvsL_Acp = new TH2D("hdcpvsL_Acp","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    TH2D *hdcpvsL_Acp2invL2 = new TH2D("hdcpvsL_Acp2invL2","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    TH2D *hdcpvsL_Acp2invL = new TH2D("hdcpvsL_Acp2invL","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    TH2D *hdcpvsL_Acp2invLPnue = new TH2D("hdcpvsL_Acp2invLPnue","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    //acp/stat. eror
    TH2D *hdcpvsL_chisq0 = new TH2D("hdcpvsL_chisq0","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    //absolute Acp
    TH2D *hdcpvsL_chisq0sub = new TH2D("hdcpvsL_chisq0sub","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    
    //chisq2 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.053015
    TH2D *hdcpvsL_chisqs1sub = new TH2D("hdcpvsL_chisq1sub","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    TH2D *hdcpvsL_chisqs2sub = new TH2D("hdcpvsL_chisq2sub","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    TH2D *hdcpvsL_chisqs3sub = new TH2D("hdcpvsL_chisq3sub","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    
    TH2D *hdcpvsL_mu2e_ih = new TH2D("hdcpvsL_mu2e_ih","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    TH2D *hdcpvsL_mu2ebar_ih = new TH2D("hdcpvsL_mu2ebar_ih","",nstepL,Lmin4hist,Lmax4hist,nstepdcp,dcpmin4hist,dcpmax4hist);
    
    
    
    
    
    double probvaltmp;
    double probvaltmpAnti;
    double tmpAcp;
    double tmpAcp2invL2;
    double tmpAcp2invL;
    double tmpAcp2invLPnue;
    double tmpchisq0;
    
    double tmpchisq1subdcp0;
    double tmpchisq1subdcppi;
    
    double tmpchisq2subdcp0;
    double tmpchisq2subdcppi;
    
    double dcp0tmpchisq0;
    double dcppitmpchisq0;
    
    double dcp0probtmp;
    double dcp0probtmpAnti;
    double dcppiprobtmp;
    double dcppiprobtmpAnti;
    
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
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        
        myProb.SetdCP(0.0);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        
        isAnti = 1;
        dcp0probtmp = myProb.P(isAnti*14,isAnti*12,myEne);
        isAnti = -1;
        dcp0probtmpAnti = myProb.P(isAnti*14,isAnti*12,myEne);
        
        
        
        //for pi
        myProb.Reset();
        myProb.SetL(thL);
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        
        myProb.SetdCP(TMath::Pi());
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        
        isAnti = 1;
        dcppiprobtmp = myProb.P(isAnti*14,isAnti*12,myEne);
        isAnti = -1;
        dcppiprobtmpAnti = myProb.P(isAnti*14,isAnti*12,myEne);
        
        
        
        
        for (Int_t iy=0; iy<nstepdcp; iy++) {
            double thdcp = dcpmin4hist + (iy+0.5)*(dcpmax-dcpmin)/nstepdcp;
            
            //for T2K, normal ordering
            myProb.Reset();
            myProb.SetL(thL);
            myProb.SetRho(2.71);
            myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
            myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
            myProb.SetTh12(myTheta12);
            myProb.SetTh23(myTheta23);
            myProb.SetTh13(myTheta13);
            
            myProb.SetdCP(thdcp);
            myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
            myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
            
            isAnti = 1;
            probvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
            hdcpvsL_mu2e->SetBinContent(ix+1,iy+1,probvaltmp);
          
            
            isAnti = -1;
            probvaltmpAnti = myProb.P(isAnti*14,isAnti*12,myEne);
            hdcpvsL_mu2ebar->SetBinContent(ix+1,iy+1,probvaltmpAnti);
            
            tmpAcp = (probvaltmp-probvaltmpAnti)/(probvaltmp+probvaltmpAnti);
            hdcpvsL_Acp->SetBinContent(ix+1,iy+1,tmpAcp);
            
            tmpAcp2invL = tmpAcp/thL;
            hdcpvsL_Acp2invL->SetBinContent(ix+1,iy+1,tmpAcp2invL);
            
            tmpAcp2invL2 = tmpAcp/(thL*thL);
            hdcpvsL_Acp2invL2->SetBinContent(ix+1,iy+1,tmpAcp2invL2);
            
            tmpAcp2invLPnue = tmpAcp*sqrt(probvaltmp+probvaltmpAnti)/thL;
            hdcpvsL_Acp2invLPnue->SetBinContent(ix+1,iy+1,tmpAcp2invLPnue);
            
            tmpchisq0= (1./pow(thL,2.))*pow(probvaltmp-probvaltmpAnti,2)/(probvaltmp+probvaltmpAnti);
            hdcpvsL_chisq0->SetBinContent(ix+1,iy+1,tmpchisq0);
            
                        
            dcp0tmpchisq0= pow(probvaltmp-probvaltmpAnti-dcp0probtmp+dcp0probtmpAnti,2);
            tmpchisq2subdcp0 = (1./pow(thL,2.))*dcp0tmpchisq0/(probvaltmp+probvaltmpAnti);

            dcppitmpchisq0= pow(probvaltmp-probvaltmpAnti-dcppiprobtmp+dcppiprobtmpAnti,2);
            tmpchisq2subdcppi = (1./pow(thL,2.))*dcppitmpchisq0/(probvaltmp+probvaltmpAnti);
            
            tmpMinchisq = TMath::Min(tmpchisq2subdcp0,tmpchisq2subdcppi);
            hdcpvsL_chisqs2sub->SetBinContent(ix+1,iy+1,tmpMinchisq);
            
            
            tmpchisq1subdcp0 =(1./pow(thL,2.))*(pow(probvaltmp-dcp0probtmp,2)/probvaltmp+pow(probvaltmpAnti-dcp0probtmpAnti,2)/probvaltmpAnti);
            tmpchisq1subdcppi = (1./pow(thL,2.))*(pow(probvaltmp-dcppiprobtmp,2)/probvaltmp+pow(probvaltmpAnti-dcppiprobtmpAnti,2)/probvaltmpAnti);
            tmpMinchisq = TMath::Min(tmpchisq1subdcp0,tmpchisq1subdcppi);
            hdcpvsL_chisqs1sub->SetBinContent(ix+1,iy+1,tmpMinchisq);
            
            
            tmpchisq1subdcp0 =(1./pow(thL,2.))*pow(probvaltmp*dcp0probtmpAnti-probvaltmpAnti*dcp0probtmp,2)/(probvaltmp*pow(dcp0probtmpAnti,2)+probvaltmpAnti*pow(dcp0probtmp,2));
            tmpchisq1subdcppi = (1./pow(thL,2.))*pow(probvaltmp*dcppiprobtmpAnti-probvaltmpAnti*dcppiprobtmp,2)/(probvaltmp*pow(dcppiprobtmpAnti,2)+probvaltmpAnti*pow(dcppiprobtmp,2));
            
            tmpMinchisq = TMath::Min(tmpchisq1subdcp0,tmpchisq1subdcppi);
            hdcpvsL_chisqs3sub->SetBinContent(ix+1,iy+1,tmpMinchisq);
            
            //take inverted mass ordering
            myProb.Reset();
            myProb.SetL(thL);
            myProb.SetRho(2.71);
            myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
            myProb.SetDmsq32(myDeltam2Inv*TMath::Power(10,-3.));//inverted
            myProb.SetTh12(myTheta12);
            myProb.SetTh23(myTheta23);
            myProb.SetTh13(myTheta13);
            
            myProb.SetdCP(thdcp);
            myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
            myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
            
            isAnti = 1;
           
            
            probvaltmp = myProb.P(isAnti*14,isAnti*12,myEne);
            hdcpvsL_mu2e_ih->SetBinContent(ix+1,iy+1,probvaltmp);
           
            
            isAnti = -1;

            
            probvaltmpAnti = myProb.P(isAnti*14,isAnti*12,myEne);
            hdcpvsL_mu2ebar_ih->SetBinContent(ix+1,iy+1,probvaltmpAnti);
           
            
        }
    }
    
    TFile *poutfile = new TFile(Form("oscprob_%s_%s_%s.root",savename.Data(),nufitver.Data(),strEnergy.Data()),"RECREATE");
    hdcpvsL_mu2e->Write();
    hdcpvsL_mu2ebar->Write();
    hdcpvsL_mu2e_ih->Write();
    hdcpvsL_mu2ebar_ih->Write();
   
    hdcpvsL_Acp->Write();
    hdcpvsL_Acp2invL2->Write();
    hdcpvsL_Acp2invL->Write();
    hdcpvsL_Acp2invLPnue->Write();
    hdcpvsL_chisq0->Write();
    hdcpvsL_chisqs1sub->Write();
    hdcpvsL_chisqs2sub->Write();
    hdcpvsL_chisqs3sub->Write();
    
    poutfile->Close();
    
    
}
