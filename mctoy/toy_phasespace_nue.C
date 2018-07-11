//
//  toy_phasespace.C
////////////////////////////////////////////////////
//
//  Simple simulation of numu + n --> p + mu-
//  using TGenPhaseSpace and store basic information
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
{
    //this is to load style for plot.
    gROOT->ProcessLine(".x rootlogon.C");
    
    gSystem.Load("libPhysics");
    double electronmass = 0.5109989461*1e-3; // electron mass in GeV
    double protonmass = 0.93827208;
    double neutronmass = 0.93956542;
    double nuE = 0.6;//Neutrino energy fixed
    //define the 4-momentum vector for neutrino as beam and neutron as target
    TLorentzVector target(0.0, 0.0, 0.0, neutronmass);//neutron, assumed at rest
    TLorentzVector beam(0.0, 0.0, nuE, nuE);//neutrino
    //invariant mass
    TLorentzVector W = beam + target;
    //mass of two produced particle
    Double_t masses[2] = { protonmass, electronmass} ;
    //define a TGenPhaseSpace object
    TGenPhaseSpace event;
    event.SetDecay(W,2,masses);
    //histogram to store basic info. for lepton
    TH1D *helectron_momentum = new TH1D("helectron_momentum", "", 100, 0, 1);
    TH1D *helectron_angle = new TH1D("helectron_angle", "", 180, 0, 180);
    TH2D *helectron_momvsangle = new TH2D("helectron_momvsangle", "", 100, 0, 1,180,0,180);
    
    //store basic info. in a tree
    TFile *pfile_ouput = new TFile("toy_phasespace_ccqe_nue.root","RECREATE");
    TTree *ptree = new TTree("event","event");
    double Enu;
    double Emu,Pmu,Thetamu;
    double Epro,Ppro, Thetapro;
    
    ptree->Branch("enu",&Enu,"enu/D");
    ptree->Branch("emu",&Emu,"emu/D");
    ptree->Branch("pmu",&Pmu,"pmu/D");
    ptree->Branch("thetamu",&Thetamu,"thetamu/D");
    
    ptree->Branch("epro",&Epro,"epro/D");
    ptree->Branch("ppro",&Ppro,"ppro/D");
    ptree->Branch("thetapro",&Thetapro,"thetapro/D");
    //define No. of events to generate
    int NEVENTGEN = 100000;
    
    for (Int_t n=0;n<NEVENTGEN;n++) {
        if(n%10000==0) cout<<n<<"/"<<NEVENTGEN<<" generated"<<endl;
        event.Generate();
        TLorentzVector *pProton = event.GetDecay(0);
        TLorentzVector *pelectron = event.GetDecay(1);
        
        //store neutrino energy
        Enu = nuE;
        
        //store basic info. for electron
        Emu = pelectron->E();
        Pmu = pelectron->P();
        Thetamu = pelectron->Theta()*180/TMath::Pi();
        
        //store basic info. for electron
        Epro = pProton->E();
        Ppro = pProton->P();
        Thetapro = pProton->Theta()*180/TMath::Pi();
        
        //store basic info. for proton
        helectron_momentum->Fill(pelectron->P());
        helectron_angle->Fill(pelectron->Theta()*180/TMath::Pi());
        helectron_momvsangle->Fill(pelectron->P(),pelectron->Theta()*180/TMath::Pi());
        ptree->Fill();
    }
    helectron_momentum->Write("helectron_momentum");
    helectron_angle->Write("helectron_angle");
    helectron_momvsangle->Write("helectron_momvsangle");
    ptree->Write();
    pfile_ouput->Close();
    
    //you can make simple plot here
    new TCanvas;
    helectron_angle->Draw();
    helectron_angle->GetXaxis()->SetTitle("Scattering angle (^{#circ}) of induced electrons");
    helectron_angle->GetYaxis()->SetTitle("Number of event generated");
    gPad->Print("electron_angle.png");
    gPad->Print("electron_angle.eps");
    
    new TCanvas;
    helectron_momentum->Draw();
    helectron_momentum->GetXaxis()->SetTitle("Momentum [GeV] of induced electrons");
    helectron_momentum->GetYaxis()->SetTitle("Number of event generated");
    gPad->Print("electron_momentum.png");
    gPad->Print("electron_momentum.eps");
    
    new TCanvas;
    //increase right margin of Pad
    gPad->SetRightMargin(gPad->GetRightMargin()*1.2);
    helectron_momvsangle->Draw("colz");
    helectron_momvsangle->GetXaxis()->SetTitle("Momentum [GeV] of induced electrons");
    helectron_momvsangle->GetYaxis()->SetTitle("Scattering angle (^{#circ}) of induced electrons");
    gPad->Print("electron_momvsangle.png");
    gPad->Print("electron_momvsangle.eps");
    
}
