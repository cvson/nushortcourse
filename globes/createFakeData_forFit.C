void createFakeData_forFit(){
    Bool_t isApplyingPoissonFluc = true;
    gRandom->SetSeed(0);
    TFile *pfile = new TFile("eventrate_all_t2k2_final_wsmear.root","READ");
    const int NRULES = 4;
    TH1D *hrate_all_rule[NRULES];
    int Nobs;
    for (Int_t irule=0; irule<NRULES; ++irule) {
        hrate_all_rule[irule] = (TH1D*)pfile->Get(Form("hrate_all_rule%d",irule));
        
        cout<<"for rule "<<irule<<endl;
        for (Int_t ibin=0; ibin<hrate_all_rule[irule]->GetNbinsX(); ++ibin) {
            Nobs = hrate_all_rule[irule]->GetBinContent(ibin+1);
            if (isApplyingPoissonFluc) {
                Nobs = gRandom->Poisson(hrate_all_rule[irule]->GetBinContent(ibin+1));
            }
            cout<<Nobs<<", ";
        }
        cout<<endl;
    }
    
}
