//
//  read_flux.C
////////////////////////////////////////////////////
//
//  To read and plot flux
//
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
void read_flux(){
    TFile *pfile = new TFile("t2kflux_2013_horn250kA.root","READ");
    const int NDET = 2;//for near detector and far detector
    char *namedet[NDET]={"nd280","sk"};
    
    const int NNUTYPE = 4;//numu, numubar, nue, nuebar
    char *namenutype[NNUTYPE]={"numu","numub","nue","nueb"};
    char *namenutype4leg[NNUTYPE]={"#nu_{#mu} flux","#bar{#nu}_{#mu} flux","#nu_e flux","#bar{#nu}_e flux"};
    //to retrive histogram
    TH1D* hflux[NDET][NNUTYPE];
    //color for histogram
    const char *colorcode[] = {
        "#000000",
        "#0072B2",
        "#D55E00",
        "#CC79A7"
    };
    Int_t ci;
    for (Int_t idet=0; idet<NDET; ++idet) {
        for (Int_t itype=0; itype<NNUTYPE; ++itype) {
            hflux[idet][itype] = (TH1D*)pfile->Get(Form("enu_%s_%s",namedet[idet],namenutype[itype]));
            hflux[idet][itype]->SetLineWidth(2);//increase thickness of histogram line
            //change color of histogram
            ci = TColor::GetColor(colorcode[itype]);
            hflux[idet][itype]->SetLineColor(ci);
            
        }
    }
    ////to add legend to the plot
    TLegend* leg0 = new TLegend(.62, .58, 0.9, .82);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(18);
    leg0->SetTextFont(43);
    for (Int_t itype=0; itype<NNUTYPE; ++itype) {
        leg0->AddEntry(hflux[0][itype],namenutype4leg[itype],"l");
    }
    
    //plot
    for (Int_t idet=0; idet<NDET; ++idet) {
        new TCanvas;
        gPad->SetLogy();// using log scale in y
        //draw 1st histogram
        hflux[idet][0]->Draw("hist");
        //in the root file, text in axis is big --> reduce it
        hflux[idet][0]->GetXaxis()->SetTitleSize(hflux[idet][0]->GetXaxis()->GetLabelSize()*0.7);
        hflux[idet][0]->GetYaxis()->SetTitleSize(hflux[idet][0]->GetYaxis()->GetLabelSize()*0.7);
        hflux[idet][0]->GetXaxis()->SetLabelSize(hflux[idet][0]->GetXaxis()->GetLabelSize()*0.7);
        hflux[idet][0]->GetYaxis()->SetLabelSize(hflux[idet][0]->GetYaxis()->GetLabelSize()*0.7);
        hflux[idet][0]->GetXaxis()->SetTitleOffset(0.9);
        hflux[idet][0]->GetYaxis()->SetTitleOffset(0.9);
        //draw other histogram
        for (Int_t itype=1; itype<NNUTYPE; ++itype) {
            hflux[idet][itype]->Draw("hist same");
        }
        //draw legend
        leg0->Draw();
        //save histogram
        gPad->Print(Form("t2kflux_%s.eps",namedet[idet]));//eps format, best for latex
        gPad->Print(Form("t2kflux_%s.png",namedet[idet]));//png format
        
    }
    
}
