
void titleStyle(TH1* h1);
void checkWithSmearing_all_final(){
    //this event rate got with GloBES before applying any post-smearing-effiency
    //should be produced with the same oscillation prob. as the publication materials which we use to modify
    Bool_t isComparetoExperimentalPred = true;
    TFile *pfile = new TFile("eventrate_all_t2k2_final_wsmear.root","READ");
    
    //app_e, app_ebar, dis_mu, dis_mubar
    const int NRULES = 4;
    TH1D *hrate_all_rule[NRULES];
    TH1D *hrate_sig_rule[NRULES];
    TH1D *hrate_bkg_rule[NRULES];
    const int NCHANNEL = 27;//defined in your T2HKall.glb file
    //7 for app_e, 8 for app_ebar, 8 for dis_mu, 8 for dis_mubar
    const char *channelname[NCHANNEL] = {
        /*For appearance nue*/
        "#nu_{#mu}#rightarrow #nu_{e} CCQE,#nu-mode",
        "#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e} CCQE, #nu-mode",
        "#nu_{#mu}#rightarrow #nu_{#mu} CC,#nu-mode",
        "#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu} CC, #nu-mode",
        "beam #nu_{e}#rightarrow #nu_{e} CC, #nu-mode",
        "beam #bar{#nu}_{e}#rightarrow #bar{#nu}_{e} CC, #nu-mode",
        "NC background, #nu-mode",
        /*For appearance nuebar*/
        "#nu_{#mu}#rightarrow #nu_{e} CCQE, #bar{#nu}-mode",
        "#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e} CCQE,#bar{#nu}-mode",
        "#nu_{#mu}#rightarrow #nu_{#mu} CC, #bar{#nu}-mode",
        "#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu} CC,#bar{#nu}-mode",
        "beam #nu_{e}#rightarrow #nu_{e} CC, #bar{#nu}-mode",
        "beam #bar{#nu}_{e}#rightarrow #bar{#nu}_{e} CC, #bar{#nu}-mode",
        "NC #nu background, #bar{#nu}-mode",
        "NC $bar{#nu} background, #bar{#nu}-mode",
        /*For disappearance numu*/
        "#nu_{#mu}#rightarrow #nu_{#mu} CCQE,#nu-mode",
        "#nu_{#mu}#rightarrow #nu_{#mu} CCnonQE, #nu-mode",
        "#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu} CCQE, #nu-mode",
        "#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu} CCnonQE, #nu-mode",
        /*"beam #nu_{e}#rightarrow #nu_{e} CC, #nu-mode",
         "beam #bar{#nu}_{e}#rightarrow #bar{#nu}_{e} CC, #nu-mode",
         "NC background, #nu-mode",*/
        "#nu_{#mu}#rightarrow #nu_{e} CC,#nu-mode",
        /*For disappearnce numubar*/
        "#nu_{#mu}#rightarrow #nu_{#mu} CCQE, #bar{#nu}-mode",
        "#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu} CCQE,#bar{#nu}-mode",
        "#nu_{#mu}#rightarrow #nu_{#mu} CCnonQE, #bar{#nu}-mode",
        "#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu} CCnonQE, #bar{#nu}-mode",
        /*"beam #nu_{e}#rightarrow #nu_{e} CC, #bar{#nu}-mode",
         "beam #bar{#nu}_{e}#rightarrow #bar{#nu}_{e} CC, #bar{#nu}-mode",
         "NC background, #bar{#nu}-mode",*/
        "#nu_{#mu}#rightarrow #nu_{e} CC, #bar{#nu}-mode",
        "NC background, #nu-mode", /*separate for disappearance*/
        "NC $bar{#nu} background, #bar{#nu}-mode" /*separate for disappearance*/
    };
    const char *colorcode[] = {
        "#000000",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#E69F00",
        "#009E73",
        "#56B4E9",
        "#F0E442",
        "#696969"
    };
    
    Int_t ci;
    //get the event rate from GloBES
    for (Int_t irule=0; irule<NRULES; ++irule) {
        hrate_all_rule[irule] = (TH1D*)pfile->Get(Form("hrate_all_rule%d",irule));
        hrate_all_rule[irule]->SetLineWidth(2);
        hrate_sig_rule[irule] = (TH1D*)pfile->Get(Form("hrate_sig_rule%d",irule));
        hrate_sig_rule[irule]->SetLineWidth(2);
        hrate_bkg_rule[irule] = (TH1D*)pfile->Get(Form("hrate_bkg_rule%d",irule));
        hrate_bkg_rule[irule]->SetLineWidth(2);
    }
    //get contribution by channel
    /*Appearance of electron neutrino*/
    //rule0 sig
    const int Rule0_NSIG = 2;
    Int_t IndexRule0_SIG[Rule0_NSIG]={0,1};
    TH1D *hrule0_sig[Rule0_NSIG];
    for (Int_t ichan=0; ichan<Rule0_NSIG; ++ichan) {
        hrule0_sig[ichan]=(TH1D*)pfile->Get(Form("hrate_sig_rule%d_chan%d_post",0,IndexRule0_SIG[ichan]));
        hrule0_sig[ichan]->SetLineWidth(2);
    }
    //rule0 bkg
    const int Rule0_NBKG = 5;
    Int_t IndexRule0_BKG[Rule0_NBKG]={2,3,4,5,6};
    TH1D *hrule0_bkg[Rule0_NBKG];
    for (Int_t ichan=0; ichan<Rule0_NBKG; ++ichan) {
        hrule0_bkg[ichan]=(TH1D*)pfile->Get(Form("hrate_bkg_rule%d_chan%d_post",0,IndexRule0_BKG[ichan]));
        hrule0_bkg[ichan]->SetLineWidth(2);
    }
    /*Appearance of electron neutrino*/
    //rule1 sig
    const int Rule1_NSIG = 2;
    Int_t IndexRule1_SIG[Rule1_NSIG]={8,7};
    TH1D *hrule1_sig[Rule1_NSIG];
    for (Int_t ichan=0; ichan<Rule1_NSIG; ++ichan) {
        hrule1_sig[ichan]=(TH1D*)pfile->Get(Form("hrate_sig_rule%d_chan%d_post",1,IndexRule1_SIG[ichan]));
        hrule1_sig[ichan]->SetLineWidth(2);
    }
    //rule1 bkg
    const int Rule1_NBKG = 6;
    Int_t IndexRule1_BKG[Rule1_NBKG]={9,10,11,12,13,14};
    TH1D *hrule1_bkg[Rule1_NBKG];
    for (Int_t ichan=0; ichan<Rule1_NBKG; ++ichan) {
        hrule1_bkg[ichan]=(TH1D*)pfile->Get(Form("hrate_bkg_rule%d_chan%d_post",1,IndexRule1_BKG[ichan]));
        hrule1_bkg[ichan]->SetLineWidth(2);
    }
    
    /*Disappewarance of muon neutrino*/
    //rule2 sig
    const int Rule2_NSIG = 4;
    Int_t IndexRule2_SIG[Rule2_NSIG]={15,16,17,18};
    TH1D *hrule2_sig[Rule2_NSIG];
    for (Int_t ichan=0; ichan<Rule2_NSIG; ++ichan) {
        hrule2_sig[ichan]=(TH1D*)pfile->Get(Form("hrate_sig_rule%d_chan%d_post",2,IndexRule2_SIG[ichan]));
        hrule2_sig[ichan]->SetLineWidth(2);
    }
    //rule2 bkg
    const int Rule2_NBKG = 4;
    Int_t IndexRule2_BKG[Rule2_NBKG]={4,5,25,19};
    TH1D *hrule2_bkg[Rule2_NBKG];
    for (Int_t ichan=0; ichan<Rule2_NBKG; ++ichan) {
        hrule2_bkg[ichan]=(TH1D*)pfile->Get(Form("hrate_bkg_rule%d_chan%d_post",2,IndexRule2_BKG[ichan]));
        hrule2_bkg[ichan]->SetLineWidth(2);
    }
    /*Disappewarance of anti-muon neutrino*/
    //rule3 sig
    const int Rule3_NSIG = 4;
    Int_t IndexRule3_SIG[Rule3_NSIG]={21,20,22,23};
    TH1D *hrule3_sig[Rule3_NSIG];
    for (Int_t ichan=0; ichan<Rule3_NSIG; ++ichan) {
        hrule3_sig[ichan]=(TH1D*)pfile->Get(Form("hrate_sig_rule%d_chan%d_post",3,IndexRule3_SIG[ichan]));
        hrule3_sig[ichan]->SetLineWidth(2);
    }
    //rule3 bkg
    const int Rule3_NBKG = 4;
    Int_t IndexRule3_BKG[Rule3_NBKG]={11,12,26,24};
    TH1D *hrule3_bkg[Rule3_NBKG];
    for (Int_t ichan=0; ichan<Rule3_NBKG; ++ichan) {
        hrule3_bkg[ichan]=(TH1D*)pfile->Get(Form("hrate_bkg_rule%d_chan%d_post",3,IndexRule3_BKG[ichan]));
        hrule3_bkg[ichan]->SetLineWidth(2);
    }
    
    
    //Get Information from Hyper-K Design report
    TFile *pfileT2K2 = new TFile("t2k2Scale_frt2k2017.root","READ");
    //get appearance histogram
    char* pfilename[]={"t2k_fhc_app_data","t2k_fhc_app_all","t2k_fhc_app_beame_nc","t2k_fhc_app_nc"
        ,"t2k_rhc_app_data", "t2k_rhc_app_all", "t2k_rhc_app_all_0nuebar", "t2k_rhc_app_all_0nuebarnue","t2k_rhc_app_nc"};
    Int_t NHISTAPP =sizeof(pfilename)/sizeof(pfilename[0]);
    TH1D **h1app = new TH1D*[NHISTAPP];
    for (Int_t ihist=0; ihist<NHISTAPP; ++ihist) {
        h1app[ihist] = (TH1D*)pfileT2K2->Get(Form("%s",pfilename[ihist]));
        h1app[ihist]->SetLineWidth(2);
        
    }
    //get disappearance histogram
    char* pfilenamedis[]={"t2k_fhc_dis_data","t2k_fhc_dis_all","t2k_fhc_dis_all_0numubar","t2k_fhc_dis_nonccqe_nc","t2k_fhc_dis_nc","t2k_rhc_dis_data","t2k_rhc_dis_all","t2k_rhc_dis_all_0numubar","t2k_rhc_dis_all_0numubarnumu","t2k_rhc_dis_nc"};
    Int_t NHISTDIS =sizeof(pfilenamedis)/sizeof(pfilenamedis[0]);
    TH1D **h1dis = new TH1D*[NHISTDIS];
    for (Int_t ihist=0; ihist<NHISTDIS; ++ihist) {
        h1dis[ihist] = (TH1D*)pfileT2K2->Get(Form("%s",pfilenamedis[ihist]));
        h1dis[ihist]->SetLineWidth(2);
        h1dis[ihist]->GetXaxis()->SetRangeUser(0,2.9);
    }
    
    
    //check disappearance sample in neutrino mode
    new TCanvas;
    gStyle->SetOptStat(0);
    TLegend* leg0 = new TLegend(.42, .58, 0.9, .9);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(20);
    leg0->SetTextFont(43);
    h1dis[1]->SetLineColor(1);
    if(isComparetoExperimentalPred)h1dis[1]->Draw("ep");//t2k2 report for FHC
    else h1dis[1]->Draw("AXIS");
    h1dis[1]->GetYaxis()->SetRangeUser(0,300);
    h1dis[1]->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
    h1dis[1]->GetYaxis()->SetTitle("Number of events");
    titleStyle(h1dis[1]);
    TH1D* hnumudisCCQE = (TH1D*)h1dis[2]->Clone("hnumudisCCQE");
    hnumudisCCQE->Add(h1dis[3],-1);
    hnumudisCCQE->SetLineColor(2);
    h1dis[4]->SetLineColor(4);
    if(isComparetoExperimentalPred){
        hnumudisCCQE->Draw("hist same");
        h1dis[4]->Draw("hist same");
        leg0->AddEntry(h1dis[0], "T2K-2 scaled, Total", "l");
        leg0->AddEntry(hnumudisCCQE, "T2K-2 scaled, #nu_{#mu}#rightarrow #nu_{#mu} CCQE", "l");
        leg0->AddEntry(h1dis[4], "T2K-2 scaled, NC", "l");
        
    }
    //total GloBES prediction
    hrate_all_rule[2]->SetLineColor(1);
    hrate_all_rule[2]->SetLineWidth(2);
    hrate_all_rule[2]->SetLineStyle(7);
    hrate_all_rule[2]->Rebin(4);
    hrate_all_rule[2]->Draw("hist same");
    leg0->AddEntry(hrate_all_rule[2], "GloBES, Total", "l");
    TH1D* hglobes_wsmear_fhc_dis_numu = (TH1D*)hrule2_sig[0]->Clone("hglobes_wsmear_fhc_dis_numu");
    leg0->AddEntry(hglobes_wsmear_fhc_dis_numu, "GloBES, #nu_{#mu}#rightarrow #nu_{#mu} CCQE", "l");
    ci = TColor::GetColor(colorcode[1]);
    hglobes_wsmear_fhc_dis_numu->SetLineColor(ci);
    hglobes_wsmear_fhc_dis_numu->SetLineStyle(7);
    hglobes_wsmear_fhc_dis_numu->Rebin(4);
    hglobes_wsmear_fhc_dis_numu->Draw("hist same");
    hrule2_bkg[2]->SetLineStyle(7);
    hrule2_bkg[2]->Rebin(4);
    ci = TColor::GetColor(colorcode[2]);
    hrule2_bkg[2]->SetLineColor(ci);
    hrule2_bkg[2]->Draw("hist same");
    leg0->AddEntry(hrule2_bkg[2], "GloBES, NC", "l");
    leg0->Draw();
    gPad->Print(isComparetoExperimentalPred?"t2k_fhc_dis_comp_wsmear_final.pdf":"t2k_fhc_dis_comp_wsmear_final_globesonly.pdf");
    delete leg0;
    
    
    //check RHC
    new TCanvas;
    TLegend* leg1 = new TLegend(.42, .58, 0.9, .9);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(20);
    leg1->SetTextFont(43);
    h1dis[6]->SetLineColor(1);
    if(isComparetoExperimentalPred)h1dis[6]->Draw("ep");//T2K scaled
    else h1dis[6]->Draw("AXIS");
    h1dis[6]->GetYaxis()->SetRangeUser(0,120);
    h1dis[6]->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
    h1dis[6]->GetYaxis()->SetTitle("Number of events");
    titleStyle(h1dis[6]);
    TH1D* hnumubarCCQE = (TH1D*)h1dis[6]->Clone("hnumubarCCQE");
    hnumubarCCQE->Add(h1dis[7],-1);
    hnumubarCCQE->SetLineColor(2);
    h1dis[9]->SetLineColor(4);
    if (isComparetoExperimentalPred) {
        hnumubarCCQE->Draw("hist same");
        h1dis[9]->Draw("hist same");
        leg1->AddEntry(h1dis[6], "T2K-2 scaled, Total", "l");
        leg1->AddEntry(hnumubarCCQE, "T2K-2 scaled, #bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu} CCQE", "l");
        leg1->AddEntry(h1dis[9], "T2K-2 scaled, NC", "l");
    }

    
    TH1D* hglobes_wsmear_rhc_dis_numub = (TH1D*)hrule3_sig[0]->Clone("hglobes_wsmear_rhc_dis_numub");
    ci = TColor::GetColor(colorcode[1]);
    hglobes_wsmear_rhc_dis_numub->SetLineColor(ci);
    hglobes_wsmear_rhc_dis_numub->SetLineStyle(7);
    
    leg1->AddEntry(hrate_all_rule[3], "GloBES, Total", "l");
    leg1->AddEntry(hglobes_wsmear_rhc_dis_numub, "GloBES, #bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu} CCQE", "l");
    hrate_all_rule[3]->SetLineColor(1);
    hrate_all_rule[3]->SetLineStyle(7);
    hrate_all_rule[3]->Rebin(4);
    hrate_all_rule[3]->Draw("hist same");
    hglobes_wsmear_rhc_dis_numub->Rebin(4);
    hglobes_wsmear_rhc_dis_numub->Draw("hist same");
    hglobes_wsmear_rhc_dis_numub->SetLineWidth(2);
    ci = TColor::GetColor(colorcode[2]);
    hrule3_bkg[2]->SetLineColor(ci);
    hrule3_bkg[2]->SetLineStyle(7);
    hrule3_bkg[2]->Rebin(4);
    hrule3_bkg[2]->Draw("hist same");
    leg1->AddEntry(hrule3_bkg[2], "GloBES, NC", "l");
    leg1->Draw();
    gPad->Print(isComparetoExperimentalPred?"t2k_rhc_dis_comp_wsmear_final.pdf":"t2k_rhc_dis_comp_wsmear_final_globesonly.pdf");
    delete leg1;
    
    
    // for appearance of nue
    new TCanvas;
    TLegend* leg2 = new TLegend(.62, .58, 0.9, .9);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(20);
    leg2->SetTextFont(43);
    leg2->SetMargin(0.1);

    h1app[1]->SetLineColor(1);
    h1app[1]->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
    h1app[1]->GetYaxis()->SetTitle("Number of events");
    titleStyle(h1app[1]);
    if(isComparetoExperimentalPred)h1app[1]->Draw("ep");
    else h1app[1]->Draw("AXIS");
    h1app[1]->GetYaxis()->SetRangeUser(0,150);
    TH1D* hnueapp = (TH1D*) h1app[1]->Clone("hnueapp");
    hnueapp->Add(h1app[2],-1);
    hnueapp->SetLineColor(2);
    h1app[3]->SetLineColor(4);
    if (isComparetoExperimentalPred) {
        hnueapp->Draw("hist same");
        h1app[3]->Draw("hist same");
        leg2->AddEntry(h1app[1], "T2K-2 scaled, Total", "l");
        leg2->AddEntry(hnueapp, "T2K-2 scaled, #nu_{#mu}#rightarrow #nu_{e} CC", "l");
        leg2->AddEntry(h1app[3], "T2K-2 scaled, NC", "l");
    }
  
    hrate_all_rule[0]->SetLineColor(1);
    hrate_all_rule[0]->SetLineStyle(7);
    hrate_all_rule[0]->Rebin(5);
    hrate_all_rule[0]->Draw("hist same");
    leg2->AddEntry(hrate_all_rule[0], "GloBES, Total", "l");
    ci = TColor::GetColor(colorcode[1]);
    hrule0_sig[0]->SetLineColor(ci);
    hrule0_sig[0]->SetLineStyle(7);
    hrule0_sig[0]->Rebin(5);
    hrule0_sig[0]->Draw("hist same");
    leg0->AddEntry(hrule0_sig[0], "GloBES,  #nu_{#mu}#rightarrow #nu_{e} CC", "l");
    ci = TColor::GetColor(colorcode[2]);
    hrule0_bkg[4]->SetLineColor(ci);
    hrule0_bkg[4]->SetLineStyle(7);
    hrule0_bkg[4]->Rebin(5);
    hrule0_bkg[4]->Draw("hist same");
    leg2->AddEntry(hrule0_bkg[4], "GloBES,  NC", "l");
    
    leg2->Draw();
    gPad->Print(isComparetoExperimentalPred?"t2k_fhc_app_comp_wsmear_final.pdf":"t2k_fhc_app_comp_wsmear_final_globesonly.pdf");
    delete leg2;
    
    // for appearance of nuebar
    new TCanvas;
    TLegend* leg3 = new TLegend(.62, .58, 0.9, .9);
    leg3->SetFillStyle(0);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(20);
    leg3->SetTextFont(43);
    leg3->SetMargin(0.1);

    h1app[5]->SetLineColor(1);
    h1app[5]->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
    h1app[5]->GetYaxis()->SetTitle("Number of events");
    titleStyle(h1app[5]);
    if(isComparetoExperimentalPred)h1app[5]->Draw("ep");
    else h1app[5]->Draw("AXIS");
    h1app[5]->GetYaxis()->SetRangeUser(0,35);
    TH1D* hnuebarapp = (TH1D*) h1app[5]->Clone("hnuebarapp");
    hnuebarapp->Add(h1app[6],-1);
    hnuebarapp->SetLineColor(2);
    h1app[8]->SetLineColor(4);
    if (isComparetoExperimentalPred) {
        hnuebarapp->Draw("hist same");
        h1app[8]->Draw("hist same");
        leg3->AddEntry(h1app[5], "T2K-2 scaled, Total", "l");
        leg3->AddEntry(hnuebarapp, "T2K-2 scaled, #bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e} CC", "l");
        leg3->AddEntry(h1app[8], "T2K-2 scaled, NC", "l");
    }
    
    hrate_all_rule[1]->SetLineColor(1);
    hrate_all_rule[1]->SetLineStyle(7);
    hrate_all_rule[1]->Rebin(5);
    hrate_all_rule[1]->Draw("hist same");
    leg3->AddEntry(hrate_all_rule[1], "GloBES, Total", "l");
    ci = TColor::GetColor(colorcode[1]);
    hrule1_sig[0]->SetLineColor(ci);
    hrule1_sig[0]->SetLineStyle(7);
    hrule1_sig[0]->Rebin(5);
    hrule1_sig[0]->Draw("hist same");
    leg3->AddEntry(hrule1_sig[0], "GloBES,  #nu_{#mu}#rightarrow #nu_{e} CC", "l");
    TH1D* hnuebarappNC = (TH1D*)hrule1_bkg[4]->Clone("hnuebarappNC");
    hnuebarappNC->Add(hrule1_bkg[5],1.0);
    ci = TColor::GetColor(colorcode[2]);
    hnuebarappNC->SetLineColor(ci);
    hnuebarappNC->SetLineStyle(7);
    hnuebarappNC->Rebin(5);
    hnuebarappNC->Draw("hist same");
    leg3->AddEntry(hnuebarappNC, "GloBES,  NC", "l");
    
    leg3->Draw();
    gPad->Print(isComparetoExperimentalPred?"t2k_rhc_app_comp_wsmear_final.pdf":"t2k_rhc_app_comp_wsmear_final_globesonly.pdf");
    delete leg3;
    
}

void titleStyle(TH1* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.2);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.2);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.2);
    
    h1->GetYaxis()->SetTitleOffset(0.85);
    h1->GetXaxis()->SetTitleOffset(0.95);
}

