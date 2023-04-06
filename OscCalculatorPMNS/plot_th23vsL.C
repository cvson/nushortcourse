#include <vector>
void SetStyleVariables(TStyle *t2kStyle);
std::vector<TGraph*> GetContourGraphs(TObjArray *contours, TString color, int graphNum);
void titleStyle2D(TH2* h1);
void titleStyle3D(TH3* h1);
void setTLatexStyle(TLatex *ptex);

void plot_th23vsL(){
    TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
    SetStyleVariables(t2kstyle);
    gROOT->SetStyle("T2K");
    Bool_t isDrawingORChi = true;
    TString subname = "th23vsL_nufit5d2_t2kene0d6";
    TFile *pfile = new TFile("/Users/cvson/meOffline/nushortcourse/OscCalculatorPMNS/oscprob_th23vsL_nufit5d2_2022nov_ene0d6.root","READ");
    
    /*TString subname = "th23vsL_nufit5d2_essnusbene0d23";
    TFile *pfile = new TFile("/Users/cvson/meOffline/nushortcourse/OscCalculatorPMNS/oscprob_th23vsL_nufit5d2_2022nov_essnusb_ene0d23.root","READ");*/
    
    char *phistname[] = {
        "hth23vsL_chisqs1nue",//0
        "hth23vsL_chisqs1nuebar",//1
        "hth23vsL_chisqs1nuesum",//2
        "hth23vsL_chisqs1numu",//3
        "hth23vsL_chisqs1numubar",//4
        "hth23vsL_chisqs1numusum",//5
        "hth23vsL_mu2e",//6
        "hth23vsL_mu2ebar",//7
        "hth23vsL_mu2mu",//8
        "hth23vsL_mu2mubar"//9
    };
    
    char *pTitlehistname[] = {
        "#chi^{2}_{1}=#frac{1}{L^{2}}#left(#frac{[P_{#mue}(#delta)-P_{#mue}(#theta_{23}=#pi/4)]^{2}}{P_{#mue}(#delta)}#right)",//0
        "#chi^{2}_{1}=#frac{1}{L^{2}}#left( #frac{[#bar{P}_{#mue}(#delta)-#bar{P}_{#mue}(#theta_{23}=#pi/4)]^{2}}{#bar{P}_{#mue}(#delta)}#right)",//1
        "#chi^{2}_{1}=#frac{1}{L^{2}}#left(#frac{[P_{#mue}(#delta)-P_{#mue}(#theta_{23}=#pi/4)]^{2}}{P_{#mue}(#delta)}+ #frac{[#bar{P}_{#mue}(#delta)-#bar{P}_{#mue}(#theta_{23}=#pi/4)]^{2}}{#bar{P}_{#mue}(#delta)}#right)",//2
        "#chi^{2}_{1}=#frac{1}{L^{2}}#left(#frac{[P_{#mu#mu}(#delta)-P_{#mu #mu}(#theta_{23}=#pi/4)]^{2}}{P_{#mu #mu}(#delta)}#right)",//3
        "#chi^{2}_{1}=#frac{1}{L^{2}}#left( #frac{[#bar{P}_{#mu #mu}(#delta)-#bar{P}_{#mu #mu}(#theta_{23}=#pi/4)]^{2}}{#bar{P}_{#mu #mu}(#delta)}#right)",//4
        "#chi^{2}_{1}=#frac{1}{L^{2}}#left(#frac{[P_{#mu #mu}(#delta)-P_{#mu #mu}(#theta_{23}=#pi/4)]^{2}}{P_{#mu #mu}(#delta)}+ #frac{[#bar{P}_{#mu #mu}(#delta)-#bar{P}_{#mu #mu}(#theta_{23}=#pi/4)]^{2}}{#bar{P}_{#mu #mu}(#delta)}#right)",//5
        "P(#nu_{#mu}#rightarrow #nu_{e}), Normal Ordering",//6
        "P(#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e}), Normal Ordering",//7
        "P(#nu_{#mu}#rightarrow #nu_{#mu}), Normal Ordering",//8
        "P(#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{#mu}), Normal Ordering"//9
    };
    // #frac{[P_{#mue}(#delta).#bar{P}_{#mue}(#delta=0,#pi)+ #bar{P}_{#mue}(#delta).P_{#mue}(#delta=0,#pi)]^{2}{P_{#mue}(#delta).#bar{P}^{2}_{#mue}(#delta=0,#pi)+ #bar{P}_{#mue}(#delta).P^{2}_{#mue}(#delta=0,#pi)}
    double ypos = 0.92;
    double xpos =0.5;
    TLatex *   text2kprecond = new TLatex(xpos,ypos-0.05,"(at L=295km, E_{#nu}=0.6 GeV)");
    setTLatexStyle(text2kprecond);
    
    Int_t NHIST =sizeof(phistname)/sizeof(phistname[0]);
    TH2D **hhist2d = new TH2D*[NHIST];
    
    char *colorcode4isoprob[]={"black","blue","white"};
    int colorindex[]={1, 4, 10};
    TLatex *textTitleHist;
    for (Int_t ihist=0; ihist<NHIST; ++ihist) {
        hhist2d[ihist] = (TH2D*)pfile->Get(Form("%s",phistname[ihist]));
        hhist2d[ihist]->GetXaxis()->SetTitle("Baseline [km]");
        hhist2d[ihist]->GetYaxis()->SetTitle("sin^{2}#theta_{23}");
        
        new TCanvas;
        titleStyle2D(hhist2d[ihist]);
        if (ihist<6) {
            hhist2d[ihist]->GetYaxis()->SetRangeUser(0.47,0.53);
        }
        hhist2d[ihist]->Draw("colz");
        textTitleHist = new TLatex(xpos,ypos,Form("%s",pTitlehistname[ihist]));
        setTLatexStyle(textTitleHist);
        //text2kprecond->Draw();
        textTitleHist->Draw();
        gPad->Print(Form("plots/%s_%s.pdf",subname.Data(),phistname[ihist]));
    }
       
    
    //octant resolving
    Int_t BinIndexTh23Max = hhist2d[6]->GetYaxis()->FindBin(0.5);
    Int_t NBINL = hhist2d[6]->GetNbinsX();
    Int_t NBINth23 = hhist2d[6]->GetNbinsY();
    
    TH2D *hORnue = (TH2D*)hhist2d[6]->Clone("hORnue");
    TH2D *hORnuebar = (TH2D*)hhist2d[6]->Clone("hORnuebar");
    TH2D *hORnuesum = (TH2D*)hhist2d[6]->Clone("hORnuesum");
    TH2D *hORnuesub = (TH2D*)hhist2d[6]->Clone("hORnuesub");
    Double_t minORChisqnue;
    Double_t minORChisqnuebar;
    Double_t tmpexp;
    Double_t tmpobs;
    Double_t tmpORChisq;
    if (isDrawingORChi) {
        for (int ix = 1; ix <= NBINL; ix++){
            for (int iy = 1; iy <= NBINth23; iy++){
                minORChisqnue=1e6;
                minORChisqnuebar = 1e6;
             
                if (iy==BinIndexTh23Max) {
                    hORnue->SetBinContent(ix,iy,0.0);
                    hORnuebar->SetBinContent(ix,iy,0.0);
                    hORnuesub->SetBinContent(ix,iy,0.0);
                    hORnuesum->SetBinContent(ix,iy,0.0);
                }
                else if (iy<BinIndexTh23Max){ //higher octant
                    for (int iyoct = BinIndexTh23Max+1; iyoct<=NBINth23; iyoct++) {
                        //nue
                        tmpexp = hhist2d[6]->GetBinContent(ix,iy);
                        tmpobs = hhist2d[6]->GetBinContent(ix,iyoct);
                        tmpORChisq = pow(tmpexp-tmpobs,2)/tmpexp;
                        if (tmpORChisq<minORChisqnue) {minORChisqnue = tmpORChisq;}
                        
                        //nuebar
                        tmpexp = hhist2d[7]->GetBinContent(ix,iy);
                        tmpobs = hhist2d[7]->GetBinContent(ix,iyoct);
                        tmpORChisq = pow(tmpexp-tmpobs,2)/tmpexp;
                        if (tmpORChisq<minORChisqnuebar) {minORChisqnuebar = tmpORChisq;}
                    }//end for
                    hORnue->SetBinContent(ix,iy,minORChisqnue);
                    hORnuebar->SetBinContent(ix,iy,minORChisqnuebar);
                    hORnuesub->SetBinContent(ix,iy,minORChisqnue-minORChisqnuebar);
                    hORnuesum->SetBinContent(ix,iy,minORChisqnue+minORChisqnuebar);
                }//end else
                else {//lower octant
                    for (int iyoctLo = 1; iyoctLo<=BinIndexTh23Max-1; iyoctLo++) {
                        //nue
                        tmpexp = hhist2d[6]->GetBinContent(ix,iy);
                        tmpobs = hhist2d[6]->GetBinContent(ix,iyoctLo);
                        tmpORChisq = pow(tmpexp-tmpobs,2)/tmpexp;
                        if (tmpORChisq<minORChisqnue) {minORChisqnue = tmpORChisq;}
                        
                        //nuebar
                        tmpexp = hhist2d[7]->GetBinContent(ix,iy);
                        tmpobs = hhist2d[7]->GetBinContent(ix,iyoctLo);
                        tmpORChisq = pow(tmpexp-tmpobs,2)/tmpexp;
                        if (tmpORChisq<minORChisqnuebar) {minORChisqnuebar = tmpORChisq;}
                    }
                    hORnue->SetBinContent(ix,iy,minORChisqnue);
                    hORnuebar->SetBinContent(ix,iy,minORChisqnuebar);
                    hORnuesub->SetBinContent(ix,iy,minORChisqnue-minORChisqnuebar);
                    hORnuesum->SetBinContent(ix,iy,minORChisqnue+minORChisqnuebar);
                }//end else
            }//end ixy
        }//end ix
        new TCanvas;
        hORnue->Draw("colz");
        gPad->Print(Form("plots/%s_ORnue.pdf",subname.Data()));
        
        new TCanvas;
        hORnuebar->Draw("colz");
        gPad->Print(Form("plots/%s_ORnuebar.pdf",subname.Data()));
        
        new TCanvas;
        hORnuesum->Draw("colz");
        gPad->Print(Form("plots/%s_ORnuesum.pdf",subname.Data()));
        
        new TCanvas;
        hORnuesub->Draw("colz");
        gPad->Print(Form("plots/%s_ORnuesub.pdf",subname.Data()));
    }//end if drawing
    
    
  
    
}

void setTLatexStyle(TLatex *ptex){
    ptex->SetTextAlign(22);
    ptex->SetNDC(kTRUE);
    ptex->SetTextFont(43);
    ptex->SetTextSize(24);
    ptex->SetLineWidth(3);
}
void titleStyle2D(TH2* h1){
    //h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.4);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.4);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetXaxis()->SetTitleOffset(0.9);
}


void titleStyle3D(TH3* h1){
    //h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetZaxis()->CenterTitle();
    
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.1);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.1);
    h1->GetZaxis()->SetLabelSize(h1->GetZaxis()->GetTitleSize()*1.1);
    
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.2);
    h1->GetZaxis()->SetTitleSize(h1->GetZaxis()->GetLabelSize()*1.2);
    
    h1->GetYaxis()->SetTitleOffset(1.2);
    h1->GetXaxis()->SetTitleOffset(1.2);
    h1->GetZaxis()->SetTitleOffset(1.2);
}


std::vector<TGraph*> GetContourGraphs(TObjArray *contours, TString color, int graphNum){
    std::vector<TGraph*> gVect;
    double x,y;
    int ncontours = contours->GetEntries();
    std::cout << ncontours << std::endl;
    TList *list = (TList*)contours->At(graphNum);
    std::cout << list->GetEntries() << std::endl;
    
    for(int j=0; j<list->GetEntries(); j++){
        TGraph *gc = new TGraph();
        
        TGraph *gc_temp = (TGraph*)list->At(j);
        int count=0;
        std::cout << "points" << gc_temp->GetN() << std::endl;
        for(int i=0; i<gc_temp->GetN(); i++){
            gc_temp->GetPoint(i, x, y);
            gc->SetPoint(i, x, y);
        }
        gc->SetLineWidth(2);
        gc->SetFillColor(0);
        gc->SetMarkerSize(0);
        
        //if(graphNum==0){ gc->SetLineStyle(7); }
        
        
        if(color.Contains("lue")){
            gc->SetLineColor(kBlue);
            gc->SetMarkerColor(kBlue);
        }
        else if(color.Contains("ed")){
            gc->SetLineColor(kRed);
            gc->SetMarkerColor(kRed);
        }
        else if(color.Contains("ack")){
            gc->SetLineColor(kBlack);
            gc->SetMarkerColor(kBlack);
        }
        else if(color.Contains("ite")){
            gc->SetLineColor(kWhite);
            gc->SetMarkerColor(kWhite);
        }
        
        gVect.push_back(gc);
    }
    return gVect;
}


void SetStyleVariables(TStyle *t2kStyle){
    
    t2kStyle->SetFrameBorderMode(0);
    t2kStyle->SetCanvasBorderMode(0);
    t2kStyle->SetPadBorderMode(0);
    t2kStyle->SetPadColor(0);
    t2kStyle->SetCanvasColor(0);
    t2kStyle->SetStatColor(0);
    t2kStyle->SetFillColor(0);
    t2kStyle->SetLegendBorderSize(1);
    
    t2kStyle->SetPaperSize(20,26);
    t2kStyle->SetPadTopMargin(0.18);
    t2kStyle->SetPadRightMargin(0.15); //0.05
    t2kStyle->SetPadBottomMargin(0.12);
    t2kStyle->SetPadLeftMargin(0.13);
    
    t2kStyle->SetTextFont(132);
    t2kStyle->SetTextSize(0.08);
    t2kStyle->SetLabelFont(132,"x");
    t2kStyle->SetLabelFont(132,"y");
    t2kStyle->SetLabelFont(132,"z");
    t2kStyle->SetLabelSize(0.05,"x");
    t2kStyle->SetTitleSize(0.06,"x");
    t2kStyle->SetLabelSize(0.05,"y");
    t2kStyle->SetTitleSize(0.06,"y");
    t2kStyle->SetLabelSize(0.05,"z");
    t2kStyle->SetTitleSize(0.06,"z");
    t2kStyle->SetLabelFont(132,"t");
    t2kStyle->SetTitleFont(132,"x");
    t2kStyle->SetTitleFont(132,"y");
    t2kStyle->SetTitleFont(132,"z");
    t2kStyle->SetTitleFont(132,"t");
    t2kStyle->SetTitleFillColor(0);
    t2kStyle->SetTitleX(0.25);
    t2kStyle->SetTitleFontSize(0.08);
    t2kStyle->SetTitleFont(132,"pad");
    
    //t2kStyle->SetPadGridX(true);
    //t2kStyle->SetPadGridY(true);
    
    
    //t2kStyle->SetHistLineWidth(1.85);
    t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    
    
    t2kStyle->SetOptTitle(0);
    t2kStyle->SetOptStat(0);
    t2kStyle->SetOptFit(0);
    
    t2kStyle->SetPadTickX(1);
    t2kStyle->SetPadTickY(1);
    
    t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                     NCont);
    t2kStyle->SetNumberContours(NCont);
    
}

