#include <vector>
void SetStyleVariables(TStyle *t2kStyle);
std::vector<TGraph*> GetContourGraphs(TObjArray *contours, TString color, int graphNum);
void titleStyle2D(TH2* h1);
void plot_dcpth23_prob(){
    TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
    SetStyleVariables(t2kstyle);
    gROOT->SetStyle("T2K");
    
    TString subname = "t2kprob_dcpth23";
    char *phistname[] = {
        "hdcpvssinsqth23_mu2mu",
        "hdcpvssinsqth23_mu2e"
    };
    
    
    Int_t NHIST =sizeof(phistname)/sizeof(phistname[0]);
    TH2D **hhist2d = new TH2D*[NHIST];
    
    TFile *pfile = new TFile("glbProb_4vson.root","READ");
    
    
    
    for (Int_t ihist=0; ihist<NHIST; ++ihist) {
        hhist2d[ihist] = (TH2D*)pfile->Get(Form("%s",phistname[ihist]));
        hhist2d[ihist]->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
        hhist2d[ihist]->GetYaxis()->SetTitle("#delta_{CP}");
    }
    
    new TCanvas;
    
    gStyle->SetOptStat(0);
    titleStyle2D(hhist2d[0]);
    hhist2d[0]->SetTitle("Survival probability P(#nu_{#mu}#rightarrow #nu_{#mu})");
    double ypos = 0.96;
    double xpos =0.5;
    TLatex *   text2kpre = new TLatex(xpos,ypos,"Survival probability P(#nu_{#mu}#rightarrow #nu_{#mu})");
    text2kpre->SetTextAlign(22);
    text2kpre->SetNDC(kTRUE);
    text2kpre->SetTextFont(43);
    text2kpre->SetTextSize(24);
    text2kpre->SetLineWidth(3);
    
    TLatex *   text2kprecond = new TLatex(xpos,ypos-0.05,"(at L=295km, E_{#nu}=0.6 GeV)");
    text2kprecond->SetTextAlign(22);
    text2kprecond->SetNDC(kTRUE);
    text2kprecond->SetTextFont(43);
    text2kprecond->SetTextSize(24);
    text2kprecond->SetLineWidth(3);
    
    const int nlevelsnumudis = 3;
    double levelsnumudis[3] = {0.02, 0.05,0.1};
    double xlatex[nlevelsnumudis]={0.44, 0.4, 0.34};
    int colorindex[nlevelsnumudis]={10, 1, 4};
    TLatex *platex[nlevelsnumudis];
    for (Int_t i=0; i<nlevelsnumudis; ++i) {
        platex[i] = new TLatex(xlatex[i],0,Form("P=%.2f",levelsnumudis[i]));
        platex[i]->SetNDC(kFALSE); // <- use NDC coordinate
        platex[i]->SetTextSize(0.04);
        platex[i]->SetTextColor(colorindex[i]);
        platex[i]->SetTextAngle(90);
    }
    
    TH2D* hnumu2numu_clone=(TH2D*)hhist2d[0]->Clone("hnumu2numu_clone");
    hnumu2numu_clone->SetContour(nlevelsnumudis,levelsnumudis);
    std::vector<TGraph*> gVect_hist_numu0;
    std::vector<TGraph*> gVect_hist_numu1;
    std::vector<TGraph*> gVect_hist_numu2;
    
    hnumu2numu_clone->Draw("cont list");
    gPad->Update();
    TObjArray *contours_histnumu = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    gVect_hist_numu0 = GetContourGraphs(contours_histnumu, "white",  0);
    gVect_hist_numu1 = GetContourGraphs(contours_histnumu, "black",  1);
    gVect_hist_numu2 = GetContourGraphs(contours_histnumu, "blue",  2);
    
    hhist2d[0]->Draw("colz");
    for(int j=0; j<gVect_hist_numu0.size(); j++){
        gVect_hist_numu0[j]->Draw("Lsame");
    }
    
    for(int j=0; j<gVect_hist_numu1.size(); j++){
        gVect_hist_numu1[j]->Draw("Lsame");
    }
    
    for(int j=0; j<gVect_hist_numu2.size(); j++){
        gVect_hist_numu2[j]->Draw("Lsame");
    }
    
    for (Int_t i=0; i<nlevelsnumudis; ++i) {
        platex[i]->Draw();
    }
    text2kpre->Draw();
    text2kprecond->Draw();
    gPad->Print(Form("%s_numu2numu.pdf",subname.Data()));
    
    
    TCanvas *c2 = new TCanvas();
    //gPad->SetBottomMargin(gPad->GetBottomMargin()*1.2);
    gStyle->SetOptStat(0);
    titleStyle2D(hhist2d[1]);
    hhist2d[1]->SetTitle("Appearance probability P(#nu_{#mu}#rightarrow #nu_{e})");
    //hhist2d[1]->Draw("colz");
    
    const int nlevels = 3;
    double levels[3] = {0.03, 0.05,0.07};
    double xlatexnue[nlevels]={0.32, 0.4, 0.57};
    double ylatexnue[nlevels]={0, -1., -1.};
    TLatex *platexnue[nlevelsnumudis];
    for (Int_t i=0; i<nlevels; ++i) {
        platexnue[i] = new TLatex(xlatexnue[i],ylatexnue[i],Form("P=%.2f",levels[i]));
        platexnue[i]->SetNDC(kFALSE); // <- use NDC coordinate
        platexnue[i]->SetTextSize(0.04);
        platexnue[i]->SetTextColor(colorindex[i]);
        platexnue[i]->SetTextAngle(35);
    }
    
    TH2D* hnumu2nue_clone=(TH2D*)hhist2d[1]->Clone("hnumu2nue_clone");
    hnumu2nue_clone->SetContour(nlevels,levels);
    TString gr_numbers[10] = {"0","1","2","3","4","5","6","7","8","9"};
    std::vector<TGraph*> gVect_hist_0;
    std::vector<TGraph*> gVect_hist_1;
    std::vector<TGraph*> gVect_hist_2;
    c2->cd();
    hnumu2nue_clone->Draw("cont list");
    c2->Update();
    TObjArray *contours_hist = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    gVect_hist_0 = GetContourGraphs(contours_hist, "white",  0);
    gVect_hist_1 = GetContourGraphs(contours_hist, "black",  1);
    gVect_hist_2 = GetContourGraphs(contours_hist, "blue",  2);
    
    
    
    hhist2d[1]->Draw("colz");
    for(int j=0; j<gVect_hist_0.size(); j++){
        gVect_hist_0[j]->Draw("Lsame");
    }
    
    for(int j=0; j<gVect_hist_1.size(); j++){
        gVect_hist_1[j]->Draw("Lsame");
    }
    
    for(int j=0; j<gVect_hist_2.size(); j++){
        gVect_hist_2[j]->Draw("Lsame");
    }
    
    TLatex *   text2kprenue = new TLatex(xpos,ypos,"Appearance probability P(#nu_{#mu}#rightarrow #nu_{e})");
    text2kprenue->SetTextAlign(22);
    text2kprenue->SetNDC(kTRUE);
    text2kprenue->SetTextFont(43);
    text2kprenue->SetTextSize(24);
    text2kprenue->SetLineWidth(3);
    text2kprenue->Draw();
    text2kprecond->Draw();
    for (Int_t i=0; i<nlevels; ++i) {
        platexnue[i]->Draw();
    }
    gPad->Print(Form("%s_numu2nue.pdf",subname.Data()));
    
    
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
    t2kStyle->SetPadTopMargin(0.12);
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
    
    
    t2kStyle->SetHistLineWidth(1.85);
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
