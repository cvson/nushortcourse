//
//  read_xsec.C
////////////////////////////////////////////////////
//
//  To read cross-section file and plot x-sec for different channels
//
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
void read_xsec(){
    TFile *pfile = new TFile("neut_nd_xsecs.root","READ");
    const int NCHANNEL = 26;
    TH1D *hxsec[NCHANNEL];
    char *namechan[NCHANNEL]={"tot","ccqe","ccppip","ccppi0","ccnpip","cccoh","ccgam","ccmpi","cceta","cck","ccdis","ncnpi0","ncppi0","ncppim","ncnpip","nccoh","ncngam","ncpgam","ncmpi","ncneta","ncpeta","nck0","nckp","ncdis","ncqep","ncqen"};
    for (Int_t ichan=0; ichan<NCHANNEL; ++ichan) {
        //This is for interaction of muon neutrino on Carbon, a.k.a C
        hxsec[ichan] = (TH1D*)pfile->Get(Form("C_xsec_numu_%s",namechan[ichan]));
        cout<<"channel "<<ichan<<" inte "<<hxsec[ichan]->Integral()<<endl;
    }
    const int NCHAN4SHOW = 6;
    char *namechan4show[NCHAN4SHOW] = {"Total", "CCQE","CC1#pi","CC coherent","CC other","NC"};
    TH1D *hxsec4show[NCHAN4SHOW];
    hxsec4show[0] = (TH1D*)hxsec[0]->Clone("htotal");
    hxsec4show[1] = (TH1D*)hxsec[1]->Clone("hccqe");
    
    hxsec4show[2] = (TH1D*)hxsec[2]->Clone("hcc1pi");
    hxsec4show[2]->Add(hxsec[3]);
    hxsec4show[2]->Add(hxsec[4]);
    
    hxsec4show[3] = (TH1D*)hxsec[5]->Clone("hcccoh");
    hxsec4show[4] = (TH1D*)hxsec[6]->Clone("hccother");
    for (Int_t ichan=7; ichan<11; ++ichan) {
        hxsec4show[4]->Add(hxsec[ichan]);
    }
    
    hxsec4show[5] = (TH1D*)hxsec[11]->Clone("hnc");
    for (Int_t ichan=12; ichan<26; ++ichan) {
        hxsec4show[5]->Add(hxsec[ichan]);
    }
    //color for histogram
    const char *colorcode[] = {
        "#000000",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#E69F00",
        "#009E73"
    };
    Int_t ci;
    for (Int_t ichan=0; ichan<NCHAN4SHOW; ++ichan) {
        hxsec4show[ichan]->SetLineWidth(2);//increase line widht, default is 1
        ci = TColor::GetColor(colorcode[ichan]);
        hxsec4show[ichan]->SetLineColor(ci);
    }

    //to plot
    new TCanvas;
    gStyle->SetOptStat(0);
    hxsec4show[0]->GetXaxis()->SetRangeUser(0,5);//check from 0-5 GeV
    hxsec4show[0]->Draw("hist");
    for (Int_t ichan=1; ichan<NCHAN4SHOW; ++ichan) {
        hxsec4show[ichan]->Draw("hist same");
    }
    //to add legend to the plot
    TLegend* leg0 = new TLegend(.22, .58, 0.5, .82);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(18);
    leg0->SetTextFont(43);
    for (Int_t ichan=0; ichan<NCHAN4SHOW; ++ichan) {
        leg0->AddEntry(hxsec4show[ichan],namechan4show[ichan],"l");
    }
    leg0->Draw();
    gPad->Print("C_xsec_numu.eps");
    
    //to check consistency
    //check if the sum of all channel equal to the total one
    TH1D* hsum = (TH1D*)hxsec4show[1]->Clone("hsum");
    for (Int_t ichan=2; ichan<NCHAN4SHOW; ++ichan) {
        hsum->Add(hxsec4show[ichan]);
    }
    new TCanvas;
    hxsec4show[0]->Draw("hist");
    hsum->Draw("hist same");
    gPad->Print("C_xsec_numu_checksum.eps");
    
    
    
}
