//
//  event_generator_frflux.C
////////////////////////////////////////////////////
//
//  To get event rate from production of flux and xsec
//   extension: include detection efficiency
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
void event_generator_frflux(){
    //get flux file
    TFile *pfile_flux = new TFile("../expflux/T2Kflux2013/t2kflux_2013_horn250kA.root","READ");
    const int NDET = 2;//for near detector and far detector
    char *namedet[NDET]={"nd280","sk"};
    
    const int NNUTYPE = 4;//numu, numubar, nue, nuebar
    char *namenutype[NNUTYPE]={"numu","numub","nue","nueb"};
    char *namenutype4leg[NNUTYPE]={"#nu_{#mu} flux","#bar{#nu}_{#mu} flux","#nu_e flux","#bar{#nu}_e flux"};
    //to retrive histogram
    TH1D* hflux[NDET][NNUTYPE];
    
    Int_t ci;
    for (Int_t idet=0; idet<NDET; ++idet) {
        for (Int_t itype=0; itype<NNUTYPE; ++itype) {
            hflux[idet][itype] = (TH1D*)pfile_flux->Get(Form("enu_%s_%s",namedet[idet],namenutype[itype]));
        }//end itype
    }//end idet
    
    TFile *pfile_op = new TFile("EVENT ","RECREATE");
    
    TTree *ptree = new TTree("Event","Neutrino event");
    
    tree->Branch("enu",&Enu,"enu/D");
    tree->Branch("pm",&pm,"pm/D");
    tree->Branch("cos",&xcos,"cos/D");
    tree->Branch("r",&r,"r/D");
    tree->Branch("dp",&dp,"dp/D");
    tree->Branch("tp",&tp,"tp/D");
    tree->Branch("crosx",&crosx,"crosx/D");
    tree->Branch("ntries",&nt,"ntries/I");
    tree->Branch("tpmin",&tpmin,"tpmin/D");
    
    tree->Branch("Pnu",Pnu,"Pnu[4]/D");
    tree->Branch("Pl",Pl,"Pl[4]/D");
    tree->Branch("Pt",Pt,"Pt[4]/D");
    tree->Branch("Pp",Pp,"Pp[4]/D");
    
    tree->Branch("Nuc",&Nuc,"Nuc/I");
    double electronmass = 0.000510998910;  // lepton mass in GeV
    double muonmass = 0.105658357; // muon mass in GeV
    double taumass = 1.77682; // tau mass in GeV
    
    double protonmass = 0.93827208;
    double neutronmass = 0.93956542;
    
    
    int NEVENTGEN = 1000;
    double Enu;
     double Pnu[4],Pt[4],Pl[4],Pp[4];
    for (int ievent=0; ievent<NEVENTGEN; ++ievent) {
        Enu = hflux[0][0]->GetRandom();
        double pnu[4],p[4][4];
        //int idpart[4],parent[4];
        //double R;
        
        pnu[1] = pnu[2] = 0;
        pnu[0] = pnu[3] = Enu;
        //  n1p1h.GenerateVectors(14,Nuc,pnu,p,idpart,parent,R,crosx);
          double Enu2 = Enu*Enu;
        double muonmass2 = muonmass*muonmass;
        
        //generate CCQE for simple
        //proton is at rest
        //numu +n = mu- + p
        //lepton kinetic energy
        double Emumax = Enu;
        double Emumin = muonmass;
        
        double Tlepmax = Emumax-muonmass;
        double Tlepmin = Emumin-muonmass;
        
        //randomize lepton momentum
        double Tlep = Tlepmin+(Tlepmax-Tlepmin)*Random();
        
        //then randomize angle
        
        ptree->Fill();
    }
    ptree->Write();
    pfile_op->Close();
    
}
