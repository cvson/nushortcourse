//To generate Tree for 3 mixing angle

{
 bool isNormalHierarchy=false;
    TString filename = isNormalHierarchy?"EventToy_normal.root":"EventToy_invert.root";
    TFile ffile(filename,"RECREATE");
    TTree* tree = new TTree("tree","tree") ;
    Double_t* sinsq12 = new Double_t ;
    Double_t* sinsq23 = new Double_t ;
    Double_t* sinsq13 = new Double_t ;
    Double_t* cpphase = new Double_t ;
    tree->Branch("sinsq12",sinsq12,"sinsq12/D") ;
    tree->Branch("sinsq23",sinsq23,"sinsq23/D") ;
    tree->Branch("sinsq13",sinsq13,"sinsq13/D") ;
    tree->Branch("cpphase",cpphase,"cpphase/D") ;
    double gsinsq12=0.307;
    double gsinsq12va=0.0170294;
    
    double gsinsq23= isNormalHierarchy?0.386:0.392;
    double gsinsq23va= isNormalHierarchy?0.0225499:0.0316623;
    
    double gsinsq13= isNormalHierarchy?0.0241:0.0244;
    double gsinsq13va= isNormalHierarchy?0.0025:0.00240208;

    double gcp= isNormalHierarchy?1.08:1.09;
    double gcpva= isNormalHierarchy?0.295381:0.325576;
    
    for (int i=0 ; i<100000 ; i++) {
        *sinsq12 = gRandom->Gaus(gsinsq12,gsinsq12va) ;
        *sinsq23 = gRandom->Gaus(gsinsq23,gsinsq23va) ; ;
        *sinsq13 = gRandom->Gaus(gsinsq13,gsinsq13va) ;
        *cpphase = gRandom->Gaus(gcp,gcpva) ;
        tree->Fill() ;
    }
    tree->Write();
    ffile.Close();
    
}