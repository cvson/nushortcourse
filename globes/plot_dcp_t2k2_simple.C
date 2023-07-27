void titleStyle(TH1* h1);
void plot_dcp_t2k2_simple(){
    char *pfilename[] = {
        "t2k2_final_nosmear_sensi_cp_simple",
        "t2k2_final_wsmear_sensi_cp_simple"
    };

    char *plegname[] = {
        "T2K-2 wo/ smearing",
        "T2K-2 w/ smearing "
    };
    
    Int_t NCOMP =sizeof(pfilename)/sizeof(pfilename[0]);
    
    char *phistname[] = {
        "hchisq_min_glob_proj",/*for projection without knowing the mass hierarchy*/
        "hchisq_min_nh_proj",/*projection with the mass hierachy known*/
        "hchisq_min_ih_proj",
        "hchisq_min_glob_sys",
        "hchisq_min_nh_nh_sys"
    };
    
    Int_t NHIST =sizeof(phistname)/sizeof(phistname[0]);
    TFile **pfile = new TFile*[NCOMP];
    TH1D *hhist[4][4];
    TLegend* leg0 = new TLegend(.42, .62, 0.9, .92);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(18);
    leg0->SetTextFont(43);
    
    
    const char *colorcode[] = { "#000000", "#0072B2","#D55E00","#CC79A7","#E69F00","#009E73","#56B4E9","#F0E442","#696969"};
    Int_t ci;
    
    for (Int_t icomp=0; icomp<NCOMP; ++icomp) {
        pfile[icomp] = new TFile(Form("%s.root",pfilename[icomp]),"READ");
        for (Int_t ihist=0; ihist<NHIST; ++ihist) {
            hhist[icomp][ihist] = (TH1D*)pfile[icomp]->Get(Form("%s",phistname[ihist]));
            ci = TColor::GetColor(colorcode[icomp+1]);
            hhist[icomp][ihist]->SetLineColor(ci);
            hhist[icomp][ihist]->SetLineWidth(2);
        }
        leg0->AddEntry(hhist[icomp][0], Form("%s",plegname[icomp]), "l");
        
    }
    
    //T2K2 official
    TFile *pfilet2k2 = new TFile("t2k2Scale_frt2k2017.root","READ");
    char* histnamecp[] = {"t2k2_dcpexcl_knownmh_th23e0d5","t2k2_dcpexcl_knownmh_th23e0d43","t2k2_dcpexcl_knownmh_th23e0d6","t2k2_dcpexcl_unknownmh_th23e0d5","t2k2_dcpexcl_unknownmh_th23e0d43","t2k2_dcpexcl_unknownmh_th23e0d6"};
    
    Int_t NHISTDCP =sizeof(histnamecp)/sizeof(histnamecp[0]);
    TGraph **dcpgraph = new TGraph*[NHISTDCP];
    for (Int_t ihist=0; ihist<NHISTDCP; ++ihist) {
        dcpgraph[ihist] = (TGraph*)pfilet2k2->Get(Form("%s",histnamecp[ihist]));
    }
    
    //convert TGraph to TH1 using Eval
    TH1D **dcphist = new TH1D*[NHISTDCP];
    for (Int_t ihist=0; ihist<NHISTDCP; ++ihist) {
        dcphist[ihist] = (TH1D*)hhist[0][0]->Clone(Form("dcphist_%d",ihist));
        for (Int_t ibin=1; ibin<=hhist[0][0]->GetNbinsX(); ++ibin) {
            double xcenter = (hhist[0][0]->GetBinCenter(ibin))*180./TMath::Pi();// convert from rad to degree
            
            double xchisqval = dcpgraph[ihist]->Eval(xcenter,0,"");//to convert to sigma
            if (xchisqval<0) {
                xchisqval=0;
            }
            if(ihist==0)cout<<"bin "<<ibin<<" center "<<xcenter<<" val "<<TMath::Sqrt(xchisqval)<<endl;
            dcphist[ihist]->SetBinContent(ibin,TMath::Sqrt(xchisqval));
        }
        ci = TColor::GetColor(colorcode[0]);
        dcphist[ihist]->SetLineColor(ci);
        
    }
    //leg0->AddEntry(dcphist[0], "T2K-2 paper", "l");
    
    for (Int_t ihist=0; ihist<NHIST; ++ihist) {
        new TCanvas;
        gStyle->SetOptStat(0);
        
        hhist[0][ihist]->Draw("hist");
        hhist[0][ihist]->GetYaxis()->SetRangeUser(0,12);
        hhist[0][ihist]->GetYaxis()->SetTitle("#sqrt{#Delta #chi^{2}}");
        hhist[0][ihist]->GetXaxis()->SetTitle("True values of #delta_{CP} [rad]");
        titleStyle(hhist[0][ihist]);
        for (Int_t icomp=1; icomp<NCOMP; ++icomp) {
            hhist[icomp][ihist]->Draw("hist same");
        }
        leg0->Draw();
        gPad->Print(Form("%s_simple_comp.pdf",phistname[ihist]));
        
    }
    
    //for compare NO and IO
    new TCanvas;
    ci = TColor::GetColor(colorcode[0]);
    hhist[1][2]->SetLineColor(ci);
    
    ci = TColor::GetColor(colorcode[1]);
    hhist[1][3]->SetLineColor(ci);
    hhist[1][3]->GetYaxis()->SetRangeUser(0,12);
    hhist[1][3]->GetYaxis()->SetTitle("#sqrt{#Delta #chi^{2}}");
    hhist[1][3]->GetXaxis()->SetTitle("Test values of #delta_{CP} [rad]");
    titleStyle(hhist[1][3]);
    hhist[1][3]->Draw("hist");
    hhist[1][2]->Draw("hist same");
    TLegend* leg1 = new TLegend(.22, .62, 0.42, .92);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(18);
    leg1->SetTextFont(43);
    leg1->AddEntry(hhist[1][2], "NO assumed", "l");
    leg1->AddEntry(hhist[1][3], "NO assumed", "l");
    leg1->Draw();
    gPad->Print(Form("%s_NOvsIO_comp.pdf",pfilename[1]));
    
}

void titleStyle(TH1* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.2);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.2);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.2);
    
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetXaxis()->SetTitleOffset(0.94);
}


