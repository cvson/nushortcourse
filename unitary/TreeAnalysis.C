#define TreeAnalysis_cxx
#include "TreeAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TComplex.h>

void TreeAnalysis::Loop()
{
bool isNormalHierarchy=false;
    TString filename = isNormalHierarchy?"HistToy_normal.root":"HistToy_invert.root";
    TFile ffile(filename,"RECREATE");
    TH1D *hs12 = new TH1D("hs12","hs12",100,0,1);
    TH1D *hc12 = new TH1D("hc12","hc12",100,0,1);
    
    TH1D *hs23 = new TH1D("hs23","hs23",100,0,1);
    TH1D *hc23 = new TH1D("hc23","hc23",100,0,1);
    
    TH1D *hs13 = new TH1D("hs13","hs13",100,0,1);
    TH1D *hc13 = new TH1D("hc13","hc13",100,0,1);
    
    TH1D *hscp = new TH1D("hscp","hscp",200,-1,1);
    TH1D *hccp = new TH1D("hccp","hccp",200,-1,1);
    
    TH1D *hUe1 = new TH1D("hUe1","hUe1",200,-1,1);
    TH1D *hUe2 = new TH1D("hUe2","hUe2",200,-1,1);
    TH1D *hUe3Re = new TH1D("hUe3Re","hUe3Re",200,-1,1);
    TH1D *hUe3Im = new TH1D("hUe3Im","hUe3Im",200,-1,1);
    
    TH1D *hUmu1Re = new TH1D("hUmu1Re","hUmu1Re",200,-1,1);
    TH1D *hUmu1Im = new TH1D("hUmu1Im","hUmu1Im",200,-1,1);
    TH1D *hUmu2Re = new TH1D("hUmu2Re","hUmu2Re",200,-1,1);
    TH1D *hUmu2Im = new TH1D("hUmu2Im","hUmu2Im",200,-1,1);
    TH1D *hUmu3 = new TH1D("hUmu3","hUmu3",200,-1,1);
    
    TH1D *hUtau1Re = new TH1D("hUtau1Re","hUtau1Re",200,-1,1);
    TH1D *hUtau1Im = new TH1D("hUtau1Im","hUtau1Im",200,-1,1);
    TH1D *hUtau2Re = new TH1D("hUtau2Re","hUtau2Re",200,-1,1);
    TH1D *hUtau2Im = new TH1D("hUtau2Im","hUtau2Im",200,-1,1);
    TH1D *hUtau3 = new TH1D("hUtau3","hUtau3",200,-1,1);
    
    TH1D *hUe1Umu1Re = new TH1D("hUe1Umu1Re","hUe1Umu1Re",200,-1,1);
    TH1D *hUe1Umu1Im = new TH1D("hUe1Umu1Im","hUe1Umu1Im",200,-1,1);
    
    TH1D *hUe2Umu2Re = new TH1D("hUe2Umu2Re","hUe2Umu2Re",200,-1,1);
    TH1D *hUe2Umu2Im = new TH1D("hUe2Umu2Im","hUe2Umu2Im",200,-1,1);
    
    TH1D *hUe3Umu3Re = new TH1D("hUe3Umu3Re","hUe3Umu3Re",200,-1,1);
    TH1D *hUe3Umu3Im = new TH1D("hUe3Umu3Im","hUe3Umu3Im",200,-1,1);
    
    TH1D *hsideARe = new TH1D("hsideARe","hsideARe",200,-1,1);
    TH1D *hsideAIm = new TH1D("hsideAIm","hsideAIm",200,-1,1);
    
    TH1D *hsideBRe = new TH1D("hsideBRe","hsideBRe",200,-1,1);
    TH1D *hsideBIm = new TH1D("hsideBIm","hsideBIm",200,-1,1);
    
    TH2D *hCorrelationA = new TH2D("hCorrelationA","hCorrelationA",200,-1,1,200,-1,1);
    TH2D *hCorrelationB = new TH2D("hCorrelationB","hCorrelationB",200,-1,1,200,-1,1);
    
    TH1D *hrhoA = new TH1D("hrhoA","hrhoA",200,-1,1);
    TH1D *hetaA = new TH1D("hetaA","hetaA",200,-1,1);
    
    TH1D *hrhoB = new TH1D("hrhoB","hrhoB",200,-1,1);
    TH1D *hetaB = new TH1D("hetaB","hetaB",200,-1,1);
    
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       if ((sinsq12>1) ||(sinsq23>1)||(sinsq13>1)) {
           continue;
       }
       //write in sin and cos
       double s12=TMath::Sqrt(sinsq12);
       hs12->Fill(s12);
       double c12=TMath::Sqrt(1-sinsq12);
       hc12->Fill(c12);
       
       double s23=TMath::Sqrt(sinsq23);
       hs23->Fill(s23);
       double c23=TMath::Sqrt(1-sinsq23);
       hc23->Fill(c23);
       
       double s13=TMath::Sqrt(sinsq13);
       hs13->Fill(s13);
       double c13=TMath::Sqrt(1-sinsq13);
       hc13->Fill(c13);
       
       //using TComplex for CP phase
       double scp=TMath::Sin(cpphase*TMath::Pi());
       hscp->Fill(scp);
       double ccp=TMath::Cos(cpphase*TMath::Pi());
       hccp->Fill(ccp);
       TComplex *cpPhasePos=new TComplex(ccp,scp);
       TComplex *cpPhaseNeg=new TComplex(ccp,-scp);
       
       //PMNS element
       //electron
       double Ue1=c12*c13;
       hUe1->Fill(Ue1);
       double Ue2=s12*c13;
       hUe2->Fill(Ue2);
       TComplex *Ue3= new TComplex(cpPhaseNeg->Re()*s13,cpPhaseNeg->Im()*s13);
       hUe3Re->Fill(Ue3->Re());
       hUe3Im->Fill(Ue3->Im());
       
       //muon
       TComplex *Umu1= new TComplex(-s12*c23-c12*s23*s13*(cpPhasePos->Re()), -c12*s23*s13*(cpPhasePos->Im()));
       hUmu1Re->Fill(Umu1->Re());
       hUmu1Im->Fill(Umu1->Im());
       TComplex *Umu2= new TComplex(c12*c23-s12*s23*s13*(cpPhasePos->Re()),-s12*s23*s13*(cpPhasePos->Im()));
       hUmu2Re->Fill(Umu2->Re());
       hUmu2Im->Fill(Umu2->Im());
       double Umu3=s23*c13;
       hUmu3->Fill(Umu3);
       
       //tau
       TComplex *Utau1=new TComplex(s12*s23-c12*c23*s13*(cpPhasePos->Re()),-c12*c23*s13*(cpPhasePos->Im()));
       hUtau1Re->Fill(Utau1->Re());
       hUtau1Im->Fill(Utau1->Im());
       TComplex *Utau2=new TComplex(-c12*s23-s12*c23*s13*(cpPhasePos->Re()),-s12*c23*s13*(cpPhasePos->Im()));
       hUtau2Re->Fill(Utau2->Re());
       hUtau2Im->Fill(Utau2->Im());
       double Utau3=c23*c13;
       hUtau3->Fill(Utau3);
       cout<<"Utau3  "<<Utau3<<endl;
       
       //calculate the triangle
       TComplex *Ue1Umu1=new TComplex(Ue1*Umu1->Re(),Ue1*Umu1->Im());//Ue1 is real
       hUe1Umu1Re->Fill(Ue1Umu1->Re());
       hUe1Umu1Im->Fill(Ue1Umu1->Im());
       
       TComplex *Ue2Umu2=new TComplex(Ue2*Umu2->Re(),Ue2*Umu2->Im());//Ue2 is real
       hUe2Umu2Re->Fill(Ue2Umu2->Re());
       hUe2Umu2Im->Fill(Ue2Umu2->Im());
       
       TComplex *Ue3Umu3=new TComplex((Ue3->Re())*Umu3,-(Ue3->Im())*Umu3);//conjugate of e3
       hUe3Umu3Re->Fill(Ue3Umu3->Re());
       hUe3Umu3Im->Fill(Ue3Umu3->Im());
       
       // sides
       //TComplex sideA=Ue1Umu1/Ue2Umu2;
       //TComplex sideB=Ue3Umu3/Ue2Umu2;
       TComplex *sideA=new TComplex((Ue1Umu1->Re()*Ue2Umu2->Re()+Ue1Umu1->Im()*Ue2Umu2->Im())/Ue2Umu2->Rho2(),(-Ue1Umu1->Re()*Ue2Umu2->Im()+Ue1Umu1->Im()*Ue2Umu2->Re())/Ue2Umu2->Rho2());
       hsideARe->Fill(sideA->Re());
       hsideAIm->Fill(sideA->Im());
       
        TComplex *sideB=new TComplex((Ue3Umu3->Re()*Ue2Umu2->Re()+Ue3Umu3->Im()*Ue2Umu2->Im())/Ue2Umu2->Rho2(),(-Ue3Umu3->Re()*Ue2Umu2->Im()+Ue3Umu3->Im()*Ue2Umu2->Re())/Ue2Umu2->Rho2());
       hsideBRe->Fill(sideB->Re());
       hsideBIm->Fill(sideB->Im());
       
       hCorrelationA->Fill(1+sideA->Re(),sideA->Im());
       hCorrelationB->Fill(-sideB->Re(),-sideB->Im());
       hrhoA->Fill(1+sideA->Re());
       hetaA->Fill(sideA->Im());
       hrhoB->Fill(-sideB->Re());
       hetaB->Fill(-sideB->Im());
       
       
       
      // if (Cut(ientry) < 0) continue;
       //cout<<s12<<endl;
   }
    hs12->Write();
    hc12->Write();
    
    hs23->Write();
    hc23->Write();
    
    hs13->Write();
    hc13->Write();
    
    hscp->Write();
    hccp->Write();
    
    hUe1->Write();
    hUe2->Write();
    hUe3Re->Write();
    hUe3Im->Write();
    
    hUmu1Re->Write();
    hUmu1Im->Write();
    hUmu2Re->Write();
    hUmu2Im->Write();
    hUmu3->Write();
    
    hUtau1Re->Write();
    hUtau1Im->Write();
    hUtau2Re->Write();
    hUtau2Im->Write();
    hUtau3->Write();
    
    hUe1Umu1Re->Write();
    hUe1Umu1Im->Write();
    
    hUe2Umu2Re->Write();
    hUe2Umu2Im->Write();
    
    hUe3Umu3Re->Write();
    hUe3Umu3Im->Write();
    
    hsideARe->Write();
    hsideAIm->Write();
    hsideBRe->Write();
    hsideBIm->Write();
    hCorrelationA->Write();
    hCorrelationB->Write();
    hrhoA->Write();
    hetaA->Write();
    hrhoB->Write();
    hetaB->Write();
    
    ffile.Close();
}
