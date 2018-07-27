//
//  event_rate.C
////////////////////////////////////////////////////
//
//  To get event rate from production of flux and xsec
//
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
void event_rate(){
    TFile *pfile = new TFile("xsec/neut_nd_xsecs.root","READ");
    const int NCHANNEL = 26;
    TH1D *hxsec[NCHANNEL];
    char *namechan[NCHANNEL]={"tot","ccqe","ccppip","ccppi0","ccnpip","cccoh","ccgam","ccmpi","cceta","cck","ccdis","ncnpi0","ncppi0","ncppim","ncnpip","nccoh","ncngam","ncpgam","ncmpi","ncneta","ncpeta","nck0","nckp","ncdis","ncqep","ncqen"};
    for (Int_t ichan=0; ichan<NCHANNEL; ++ichan) {
        //This is for interaction of muon neutrino on Carbon, a.k.a C
        hxsec[ichan] = (TH1D*)pfile->Get(Form("C_xsec_numu_%s",namechan[ichan]));
        //cout<<"channel "<<ichan<<" inte "<<hxsec[ichan]->Integral()<<endl;
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

   
    TFile *pfile_flux = new TFile("expflux/T2Kflux2013/t2kflux_2013_horn250kA.root","READ");
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
            
            
        }
    }
    //need to convert xsec to graph
    TGraph *grxsec[NCHAN4SHOW];
    for (Int_t ichan=0; ichan<NCHAN4SHOW; ichan++) {
        grxsec[ichan] = new TGraph(hxsec4show[ichan]);
    }
    //for one ton of H20 for near example
    Double_t NoOfNucleon_near = 1e6*6.02*1e23/18.0;
    //for 22.5 kton of H20 for far detector
    Double_t NoOfNucleon_far = 22.5e9*6.02*1e23/18.0;
    //you also needs to define efficity detection
    Double_t Effiency_near = 0.9;
    Double_t Effiency_far = 0.7;//
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
    //clone and multiple with the xsec value at the center value of energy
    TH1D *heventrate4show[NDET][NCHAN4SHOW];
    for (Int_t idet=0; idet<NDET; ++idet) {
        for (Int_t ichan=0; ichan<NCHAN4SHOW; ++ichan) {
            heventrate4show[idet][ichan] = (TH1D*)hflux[idet][0]->Clone(Form("eventrate_det%d_chan%d",idet,ichan));
            for (Int_t ibin=1; ibin<heventrate4show[idet][ichan]->GetNbinsX()+1; ++ibin) {
                double eventrate_tmp = heventrate4show[idet][ichan]->GetBinContent(ibin)*grxsec[ichan]->Eval(heventrate4show[idet][ichan]->GetBinCenter(ibin))*(idet==0?NoOfNucleon_near:NoOfNucleon_far);
                heventrate4show[idet][ichan]->SetBinContent(ibin,eventrate_tmp);
            }
            cout<<"Det. "<<idet<<" channel "<<ichan<<" No. of events: "<<heventrate4show[idet][ichan]->Integral()<<endl;
            heventrate4show[idet][ichan]->SetLineWidth(2);//increase line widht, default is 1
            ci = TColor::GetColor(colorcode[ichan]);
            heventrate4show[idet][ichan]->SetLineColor(ci);
        }
        
        
    }
    //
    cout<<"Bin width "<< heventrate4show[0][0]->GetBinWidth(1)<<endl;
    //to plot
    new TCanvas;
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    heventrate4show[0][0]->GetXaxis()->SetRangeUser(0,5);//check from 0-5 GeV
    heventrate4show[0][0]->SetTitle("Event rate at Near Detector (1ton H20)");
    heventrate4show[0][0]->GetYaxis()->SetTitle("Events [/50 MeV/1x10^{21} p.o.t]");
    heventrate4show[0][0]->GetYaxis()->SetTitleOffset(0.7);
    heventrate4show[0][0]->Draw("hist");
    for (Int_t ichan=1; ichan<NCHAN4SHOW; ++ichan) {
        heventrate4show[0][ichan]->Draw("hist same");
    }
    //to add legend to the plot
    TLegend* leg0 = new TLegend(.62, .58, 0.85, .82);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(18);
    leg0->SetTextFont(43);
    for (Int_t ichan=0; ichan<NCHAN4SHOW; ++ichan) {
        leg0->AddEntry(heventrate4show[0][ichan],namechan4show[ichan],"l");
    }
    leg0->Draw();
    gPad->Print("eventrate_1tonH20_neardet.eps");
    
    new TCanvas;
    heventrate4show[1][0]->GetXaxis()->SetRangeUser(0,5);//check from 0-5 GeV
    heventrate4show[1][0]->SetTitle("Event rate at Far Detector (22.5 kton H20)");
    heventrate4show[1][0]->GetYaxis()->SetTitle("Events [/50 MeV/1x10^{21} p.o.t]");
    heventrate4show[1][0]->GetYaxis()->SetTitleOffset(0.7);
    heventrate4show[1][0]->Draw("hist");
    for (Int_t ichan=1; ichan<NCHAN4SHOW; ++ichan) {
        heventrate4show[1][ichan]->Draw("hist same");
    }
    leg0->Draw();
    gPad->Print("eventrate_1tonH20_fardet.eps");
    
    
    
}
