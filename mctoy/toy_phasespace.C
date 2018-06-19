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
    double muonmass = 0.105658357; // muon mass in GeV
    double protonmass = 0.93827208;
    double neutronmass = 0.93956542;
    double nuE = 0.6;//Neutrino energy fixed
    //define the 4-momentum vector for neutrino as beam and neutron as target
    TLorentzVector target(0.0, 0.0, 0.0, neutronmass);//neutron, assumed at rest
    TLorentzVector beam(0.0, 0.0, nuE, nuE);//neutrino
    //invariant mass
    TLorentzVector W = beam + target;
    //mass of two produced particle
    Double_t masses[2] = { protonmass, muonmass} ;
    //define a TGenPhaseSpace object
    TGenPhaseSpace event;
    event.SetDecay(W,2,masses);
    //histogram to store basic info. for lepton
    TH1D *hmuon_momentum = new TH1D("hmuon_momentum", "", 100, 0, 1);
    TH1D *hmuon_angle = new TH1D("hmuon_angle", "", 180, 0, 180);
    TH2D *hmuon_momvsangle = new TH2D("hmuon_momvsangle", "", 100, 0, 1,180,0,180);
    
    //store basic info. in a tree
    TFile *pfile_ouput = new TFile("toy_phasespace_ccqe.root","RECREATE");
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
        TLorentzVector *pMuon = event.GetDecay(1);
        
        //store neutrino energy
        Enu = nuE;
        
        //store basic info. for muon
        Emu = pMuon->E();
        Pmu = pMuon->P();
        Thetamu = pMuon->Theta()*180/TMath::Pi();
        
        //store basic info. for muon
        Epro = pProton->E();
        Ppro = pProton->P();
        Thetapro = pProton->Theta()*180/TMath::Pi();
        
        //store basic info. for proton
        hmuon_momentum->Fill(pMuon->P());
        hmuon_angle->Fill(pMuon->Theta()*180/TMath::Pi());
        hmuon_momvsangle->Fill(pMuon->P(),pMuon->Theta()*180/TMath::Pi());
        ptree->Fill();
    }
    hmuon_momentum->Write("hmuon_momentum");
    hmuon_angle->Write("hmuon_angle");
    hmuon_momvsangle->Write("hmuon_momvsangle");
    ptree->Write();
    pfile_ouput->Close();
    
    //you can make simple plot here
    new TCanvas;
    hmuon_angle->Draw();
    hmuon_angle->GetXaxis()->SetTitle("Scattering angle (^{#circ}) of induced muons");
    hmuon_angle->GetYaxis()->SetTitle("Number of event generated");
    gPad->Print("muon_angle.png");
    gPad->Print("muon_angle.eps");
    
    new TCanvas;
    hmuon_momentum->Draw();
    hmuon_momentum->GetXaxis()->SetTitle("Momentum [GeV] of induced muons");
    hmuon_momentum->GetYaxis()->SetTitle("Number of event generated");
    gPad->Print("muon_momentum.png");
    gPad->Print("muon_momentum.eps");
    
    new TCanvas;
    //increase right margin of Pad
    gPad->SetRightMargin(gPad->GetRightMargin()*1.2);
    hmuon_momvsangle->Draw("colz");
    hmuon_momvsangle->GetXaxis()->SetTitle("Momentum [GeV] of induced muons");
    hmuon_momvsangle->GetYaxis()->SetTitle("Scattering angle (^{#circ}) of induced muons");
    gPad->Print("muon_momvsangle.png");
    gPad->Print("muon_momvsangle.eps");
    
}
