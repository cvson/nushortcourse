#include <vector>
void SetStyleVariables(TStyle *t2kStyle);
std::vector<TGraph*> GetContourGraphs(TObjArray *contours, TString color, int graphNum);
void titleStyle2D(TH2* h1);
void setTLatexStyle(TLatex *ptex);
//to check the clone solution at specific point
void plot_dcpth23_prob_root6_clonesol_v3(){
    TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
    SetStyleVariables(t2kstyle);
    gROOT->SetStyle("T2K");
    
    TString subname = "hkprob_dcpth23_clonesol_v3";
    Bool_t isDrawTheoryPrediction = false;
    Bool_t isDrawLikelihood = false;
    Bool_t isDrawSin13vsSin23 = false;
    Bool_t isDrawSin13vsdcp =false;
    Bool_t isDraw3DMarginalization = true;
    //this is with 200 bins from 0.4 to 0.6
    //TFile *pfile = new TFile("/Users/cvson/meOffline/nushortcourse/OscCalculatorPMNS/oscprob_th23vsdcp_nufit5d2_2022nov.root","READ");
    //this is with 201 bins from 0.4 to 0.6
    TFile *pfile = new TFile("/Users/cvson/meOffline/nushortcourse/OscCalculatorPMNS/oscprob_th23vsdcpvsth13_nufit5d2_2022nov_ene0d6.root","READ");
    
    /*TFile *pfile = new TFile("/Users/cvson/meOffline/nushortcourse/OscCalculatorPMNS/oscprob_th23vsdcp_v2_nufit5d2_2022nov_ene0d7.root","READ");
     subname +="_ene0d7";*/
    
    /*TFile *pfile = new TFile("/Users/cvson/meOffline/nushortcourse/OscCalculatorPMNS/oscprob_th23vsdcp_v2_nufit5d2_2022nov_ene0d5.root","READ");
     subname +="_ene0d5";*/
    
    Double_t dcp_clonesol = -1.0 * TMath::Pi()/2.;
    subname +="_dcpEmpiover2";
    /*Double_t sinsqth23_clonesol = 0.54;
     subname +="_sinsqth23E0d54";*/
    Double_t sinsqth23_clonesol = 0.52;
    subname +="_sinsqth23E0d52";
    
    Double_t sinsqth13_clonesol = 0.022;
    subname +="_sinsqth13E0d022";
    
    TMarker *mTrue = new TMarker(sinsqth23_clonesol, dcp_clonesol, 29);
    mTrue->SetMarkerColor(kBlack);
    mTrue->SetMarkerSize(1.5);
    
    TMarker *mTrue_sin13vssin23 = new TMarker(sinsqth23_clonesol, sinsqth13_clonesol, 29);
    mTrue_sin13vssin23->SetMarkerColor(kBlack);
    mTrue_sin13vssin23->SetMarkerSize(1.5);
    
    TMarker *mTrue_sin13vsdcp = new TMarker(sinsqth13_clonesol, dcp_clonesol, 29);
    mTrue_sin13vsdcp->SetMarkerColor(kBlack);
    mTrue_sin13vsdcp->SetMarkerSize(1.5);
    
    char *phistname[] = {
        "hdcpvssinsqth23_mu2mu",//0
        "hdcpvssinsqth23_mu2e",//1
        "hdcpvssinsqth23_mu2mubar",//2
        "hdcpvssinsqth23_mu2ebar",//3
        "hdcpvssinsqth23_mu2e_ih",//4
        "hdcpvssinsqth23_mu2ebar_ih",//5
        "hsinsqth13vssinsqth23_mu2e",//6
        "hsinsqth13vssinsqth23_mu2ebar",//7
        "hdcpvssinsqth13_mu2e",//8
        "hdcpvssinsqth13_mu2ebar"//9
    };
    Int_t NHIST =sizeof(phistname)/sizeof(phistname[0]);
    TH2D **hhist2d = new TH2D*[NHIST];
    
    char *p3dhistname[] ={
        "hdcpvssinsqth23vssinsqth13_mu2e",
        "hdcpvssinsqth23vssinsqth13_mu2ebar"
    };
    
    Int_t NHIST3D =sizeof(p3dhistname)/sizeof(p3dhistname[0]);
    TH3D **hhist3d = new TH3D*[NHIST3D];
    
    
    
    
    
    for (Int_t ihist=0; ihist<NHIST; ++ihist) {
        hhist2d[ihist] = (TH2D*)pfile->Get(Form("%s",phistname[ihist]));
        if(ihist<6){hhist2d[ihist]->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
            hhist2d[ihist]->GetYaxis()->SetTitle("#delta_{CP}");
        }
        else if (ihist<8){
            hhist2d[ihist]->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
            hhist2d[ihist]->GetYaxis()->SetTitle("sin^{2}#theta_{13}");
        }
        else {
            hhist2d[ihist]->GetYaxis()->SetTitle("#delta_{CP}");
            hhist2d[ihist]->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
        }
    }
    
    for (Int_t ihist3d=0; ihist3d<NHIST3D; ++ihist3d) {
        hhist3d[ihist3d] = (TH3D*)pfile->Get(Form("%s",p3dhistname[ihist3d]));
    }
    
    Int_t NBINth23 = hhist2d[0]->GetNbinsX();
    Int_t NBINdcp = hhist2d[0]->GetNbinsY();
    Int_t NBINth13 = hhist2d[8]->GetNbinsZ();
    
    new TCanvas;
    gStyle->SetOptStat(0);
    titleStyle2D(hhist2d[0]);
    //hhist2d[0]->Draw("colz");
    double ypos = 0.96;
    double xpos =0.5;
    TLatex *   text2kpre = new TLatex(xpos,ypos,"Survival probability P(#nu_{#mu}#rightarrow #nu_{#mu})");
    setTLatexStyle(text2kpre);
    
    TLatex *   text2kprecond = new TLatex(xpos,ypos-0.05,"(at L=295km, E_{#nu}=0.6 GeV)");
    setTLatexStyle(text2kprecond);
    
    double prob_numu2numu_clonesol =hhist2d[0]->GetBinContent(hhist2d[0]->FindBin(sinsqth23_clonesol,dcp_clonesol));
    cout<<" Prob. numu2numu "<<prob_numu2numu_clonesol<<endl;
    
    const int nlevelsnumudis = 3;
    /*double levelsnumudis[3] = {prob_numu2numu_clonesol-0.002, prob_numu2numu_clonesol, prob_numu2numu_clonesol+0.002};
     double xlatex[nlevelsnumudis]={0.49, 0.47, 0.44};
     int colorindex[nlevelsnumudis]={10, 1, 4};
     */
    
    double levelsnumudis[3] = {prob_numu2numu_clonesol*0.9, prob_numu2numu_clonesol, prob_numu2numu_clonesol*1.1};
    double xlatex[nlevelsnumudis]={0.51, 0.495, 0.47};
    
    int colorindex[nlevelsnumudis]={10, 1, 4};
    TLatex *platex[nlevelsnumudis];
    for (Int_t i=0; i<nlevelsnumudis; ++i) {
        platex[i] = new TLatex(xlatex[i],0,Form("P=%.3f",levelsnumudis[i]));
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
    for(int j=0; j<gVect_hist_numu0.size(); j++){gVect_hist_numu0[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_numu1.size(); j++){gVect_hist_numu1[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_numu2.size(); j++){gVect_hist_numu2[j]->Draw("Lsame");}
    
    for (Int_t i=0; i<nlevelsnumudis; ++i) {platex[i]->Draw();}
    text2kpre->Draw();
    text2kprecond->Draw();
    mTrue->Draw();
    gPad->Print(Form("plots/%s_numu2numu.eps",subname.Data()));
    
    
    
    
    
    //numu2ne
    TCanvas *c2 = new TCanvas();
    gStyle->SetOptStat(0);
    titleStyle2D(hhist2d[1]);
    
    const int nlevels = 3;
    double prob_numu2nue_clonesol =hhist2d[1]->GetBinContent(hhist2d[1]->FindBin(sinsqth23_clonesol,dcp_clonesol));
    cout<<" Prob. numu2nue "<<prob_numu2nue_clonesol<<endl;
    
    double levels[3] = {prob_numu2nue_clonesol*0.9, prob_numu2nue_clonesol, prob_numu2nue_clonesol*1.1};
    double xlatexnue[nlevels]={0.5, 0.54, 0.57};
    double ylatexnue[nlevels]={-2, -2.1, -2};
    TLatex *platexnue[nlevelsnumudis];
    for (Int_t i=0; i<nlevels; ++i) {
        platexnue[i] = new TLatex(xlatexnue[i],ylatexnue[i],Form("P=%.3f",levels[i]));
        platexnue[i]->SetNDC(kFALSE); // <- use NDC coordinate
        platexnue[i]->SetTextSize(0.04);
        platexnue[i]->SetTextColor(colorindex[i]);
        platexnue[i]->SetTextAngle(-10);
    }
    
    TH2D* hnumu2nue_clone=(TH2D*)hhist2d[1]->Clone("hnumu2nue_clone");
    hnumu2nue_clone->SetContour(nlevels,levels);
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
    hhist2d[1]->GetZaxis()->SetRangeUser(0.01,0.08);
    for(int j=0; j<gVect_hist_0.size(); j++){gVect_hist_0[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_1.size(); j++){gVect_hist_1[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_2.size(); j++){gVect_hist_2[j]->Draw("Lsame");}
    
    TLatex *   text2kprenue = new TLatex(xpos,ypos,"Appearance probability P(#nu_{#mu}#rightarrow #nu_{e}), Normal Ordering");
    setTLatexStyle(text2kprenue);
    text2kprenue->Draw();
    text2kprecond->Draw();
    for (Int_t i=0; i<nlevels; ++i) {platexnue[i]->Draw();}
    mTrue->Draw();
    
    if (isDrawTheoryPrediction) {
        Double_t valDegeneracyTheorydcp = TMath::Sin(dcp_clonesol)*TMath::Sqrt(sinsqth23_clonesol)/(TMath::Sqrt(1-sinsqth23_clonesol));
        //Double_t valDegeneracyTheorydcp = TMath::Sin(dcp_clonesol)*(sinsqth23_clonesol)/((1-sinsqth23_clonesol));
        Double_t *xvalDegTheodcp= new Double_t[NBINth23];
        Double_t *yvalDegTheodcp = new Double_t[NBINth23];
        for (Int_t ibinth23=1; ibinth23<=NBINth23+1; ibinth23++) {
            xvalDegTheodcp[ibinth23-1] = hhist2d[1]->GetXaxis()->GetBinCenter(ibinth23);
            yvalDegTheodcp[ibinth23-1] = TMath::ASin((valDegeneracyTheorydcp*TMath::Sqrt(1-xvalDegTheodcp[ibinth23-1])/TMath::Sqrt(xvalDegTheodcp[ibinth23-1])));
            //yvalDegTheodcp[ibinth23-1] = TMath::ASin((valDegeneracyTheorydcp*(1-xvalDegTheodcp[ibinth23-1])/(xvalDegTheodcp[ibinth23-1])));
            cout<<"xval "<<xvalDegTheodcp[ibinth23-1]<<" yval "<<yvalDegTheodcp[ibinth23-1]<<endl;
        }
        TGraph *pgrDegeneracyTheorydcp = new TGraph(NBINth23,xvalDegTheodcp,yvalDegTheodcp);
        pgrDegeneracyTheorydcp->SetLineColor(kRed);
        pgrDegeneracyTheorydcp->Draw("L same");
        c2->Print(Form("plots/%s_numu2nue_wtheoPred.eps",subname.Data()));
    }
    else c2->Print(Form("plots/%s_numu2nue.eps",subname.Data()));
    
    
    
    //numu2nuebar
    TCanvas *c3 = new TCanvas();
    gStyle->SetOptStat(0);
    titleStyle2D(hhist2d[3]);
    double prob_numu2nuebar_clonesol =hhist2d[3]->GetBinContent(hhist2d[3]->FindBin(sinsqth23_clonesol,dcp_clonesol));
    cout<<" Prob. numu2nuebar "<<prob_numu2nuebar_clonesol<<endl;
    
    double levelsbar[3] = {prob_numu2nuebar_clonesol*0.9, prob_numu2nuebar_clonesol, prob_numu2nuebar_clonesol*1.1};
    double xlatexnuebar[nlevels]={0.42, 0.48, 0.56};
    double ylatexnuebar[nlevels]={-1.1, -1.2, -1.2};
    TLatex *platexnuebar[nlevels];
    for (Int_t i=0; i<nlevels; ++i) {
        platexnuebar[i] = new TLatex(xlatexnuebar[i],ylatexnuebar[i],Form("P=%.3f",levelsbar[i]));
        platexnuebar[i]->SetNDC(kFALSE); // <- use NDC coordinate
        platexnuebar[i]->SetTextSize(0.04);
        platexnuebar[i]->SetTextColor(colorindex[i]);
        platexnuebar[i]->SetTextAngle(-10);
    }
    
    TH2D* hnumu2nuebar_clone=(TH2D*)hhist2d[3]->Clone("hnumu2nuebar_clone");
    hnumu2nuebar_clone->SetContour(nlevels,levelsbar);
    std::vector<TGraph*> gVect_hist_bar_0;
    std::vector<TGraph*> gVect_hist_bar_1;
    std::vector<TGraph*> gVect_hist_bar_2;
    c3->cd();
    hnumu2nuebar_clone->Draw("cont list");
    c3->Update();
    TObjArray *contours_hist_bar = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    gVect_hist_bar_0 = GetContourGraphs(contours_hist_bar, "white",  0);
    gVect_hist_bar_1 = GetContourGraphs(contours_hist_bar, "black",  1);
    gVect_hist_bar_2 = GetContourGraphs(contours_hist_bar, "blue",  2);
    
    hhist2d[3]->Draw("colz");
    hhist2d[3]->GetZaxis()->SetRangeUser(0.01,0.08);
    for(int j=0; j<gVect_hist_bar_0.size(); j++){gVect_hist_bar_0[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_bar_1.size(); j++){gVect_hist_bar_1[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_bar_2.size(); j++){gVect_hist_bar_2[j]->Draw("Lsame");}
    
    TLatex *   text2kprenuebar = new TLatex(xpos,ypos,"Appearance probability P(#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e}), Normal Ordering");
    setTLatexStyle(text2kprenuebar);
    text2kprenuebar->Draw();
    text2kprecond->Draw();
    for (Int_t i=0; i<nlevels; ++i) {platexnuebar[i]->Draw();}
    mTrue->Draw();
    gPad->Print(Form("plots/%s_numu2nuebar.eps",subname.Data()));
    
    
    //numu2nue IH
    TCanvas *c2ih = new TCanvas();
    gStyle->SetOptStat(0);
    titleStyle2D(hhist2d[4]);
    
    double prob_numu2nueIH_clonesol =hhist2d[4]->GetBinContent(hhist2d[4]->FindBin(sinsqth23_clonesol,dcp_clonesol));
    cout<<" Prob. numu2nue IH "<<prob_numu2nueIH_clonesol<<endl;
    
    double levelsIH[3] = {prob_numu2nueIH_clonesol*0.9, prob_numu2nueIH_clonesol, prob_numu2nueIH_clonesol*1.1};
    double xlatexnueIH[nlevels]={0.5, 0.54, 0.57};
    double ylatexnueIH[nlevels]={-2, -2.1, -2};
    TLatex *platexnueIH[nlevelsnumudis];
    for (Int_t i=0; i<nlevels; ++i) {
        platexnueIH[i] = new TLatex(xlatexnue[i],ylatexnue[i],Form("P=%.3f",levelsIH[i]));
        platexnueIH[i]->SetNDC(kFALSE); // <- use NDC coordinate
        platexnueIH[i]->SetTextSize(0.04);
        platexnueIH[i]->SetTextColor(colorindex[i]);
        platexnueIH[i]->SetTextAngle(-10);
    }
    
    TH2D* hnumu2nueIH_clone=(TH2D*)hhist2d[4]->Clone("hnumu2nue_clone");
    hnumu2nueIH_clone->SetContour(nlevels,levelsIH);
    std::vector<TGraph*> gVect_histIH_0;
    std::vector<TGraph*> gVect_histIH_1;
    std::vector<TGraph*> gVect_histIH_2;
    c2ih->cd();
    hnumu2nueIH_clone->Draw("cont list");
    c2ih->Update();
    TObjArray *contours_histIH = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    gVect_histIH_0 = GetContourGraphs(contours_histIH, "white",  0);
    gVect_histIH_1 = GetContourGraphs(contours_histIH, "black",  1);
    gVect_histIH_2 = GetContourGraphs(contours_histIH, "blue",  2);
    
    hhist2d[4]->Draw("colz");
    hhist2d[4]->GetZaxis()->SetRangeUser(0.01,0.08);
    for(int j=0; j<gVect_histIH_0.size(); j++){gVect_histIH_0[j]->Draw("Lsame");}
    for(int j=0; j<gVect_histIH_1.size(); j++){gVect_histIH_1[j]->Draw("Lsame");}
    for(int j=0; j<gVect_histIH_2.size(); j++){gVect_histIH_2[j]->Draw("Lsame");}
    
    TLatex *   text2kprenueIH = new TLatex(xpos,ypos,"Appearance probability P(#nu_{#mu}#rightarrow #nu_{e}), inverted ordering");
    setTLatexStyle(text2kprenueIH);
    text2kprenueIH->Draw();
    text2kprecond->Draw();
    for (Int_t i=0; i<nlevels; ++i) {platexnueIH[i]->Draw();}
    mTrue->Draw();
    gPad->Print(Form("plots/%s_numu2nue_IH.eps",subname.Data()));
    
    
    //numu2nuebar IH
    TCanvas *c3ih = new TCanvas();
    gStyle->SetOptStat(0);
    titleStyle2D(hhist2d[5]);
    double prob_numu2nuebarIH_clonesol =hhist2d[5]->GetBinContent(hhist2d[5]->FindBin(sinsqth23_clonesol,dcp_clonesol));
    cout<<" Prob. numu2nuebar IH "<<prob_numu2nuebarIH_clonesol<<endl;
    
    double levelsbarIH[3] = {prob_numu2nuebarIH_clonesol*0.9, prob_numu2nuebarIH_clonesol, prob_numu2nuebarIH_clonesol*1.1};
    /*double xlatexnuebar[nlevels]={0.42, 0.48, 0.56};
     double ylatexnuebar[nlevels]={-1, -1, -1.1};*/
    double xlatexnuebarIH[nlevels]={0.42, 0.48, 0.56};
    double ylatexnuebarIH[nlevels]={-1.1, -1.2, -1.2};
    TLatex *platexnuebarIH[nlevels];
    for (Int_t i=0; i<nlevels; ++i) {
        platexnuebarIH[i] = new TLatex(xlatexnuebarIH[i],ylatexnuebarIH[i],Form("P=%.3f",levelsbarIH[i]));
        platexnuebarIH[i]->SetNDC(kFALSE); // <- use NDC coordinate
        platexnuebarIH[i]->SetTextSize(0.04);
        platexnuebarIH[i]->SetTextColor(colorindex[i]);
        platexnuebarIH[i]->SetTextAngle(-10);
    }
    
    TH2D* hnumu2nuebarIH_clone=(TH2D*)hhist2d[5]->Clone("hnumu2nuebarIH_clone");
    hnumu2nuebarIH_clone->SetContour(nlevels,levelsbarIH);
    std::vector<TGraph*> gVect_hist_barIH_0;
    std::vector<TGraph*> gVect_hist_barIH_1;
    std::vector<TGraph*> gVect_hist_barIH_2;
    c3ih->cd();
    hnumu2nuebarIH_clone->Draw("cont list");
    c3ih->Update();
    TObjArray *contours_hist_barIH = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    gVect_hist_barIH_0 = GetContourGraphs(contours_hist_barIH, "white",  0);
    gVect_hist_barIH_1 = GetContourGraphs(contours_hist_barIH, "black",  1);
    gVect_hist_barIH_2 = GetContourGraphs(contours_hist_barIH, "blue",  2);
    
    hhist2d[5]->Draw("colz");
    hhist2d[5]->GetZaxis()->SetRangeUser(0.01,0.08);
    for(int j=0; j<gVect_hist_barIH_0.size(); j++){gVect_hist_barIH_0[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_barIH_1.size(); j++){gVect_hist_barIH_1[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_barIH_2.size(); j++){gVect_hist_barIH_2[j]->Draw("Lsame");}
    
    TLatex *   text2kprenuebarIH = new TLatex(xpos,ypos,"Appearance probability P(#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e}), Inverted Ordering");
    setTLatexStyle(text2kprenuebarIH);
    text2kprenuebarIH->Draw();
    text2kprecond->Draw();
    for (Int_t i=0; i<nlevels; ++i) {platexnuebarIH[i]->Draw();}
    mTrue->Draw();
    gPad->Print(Form("plots/%s_numu2nuebar_IH.eps",subname.Data()));
    
    //relative difference between normal and inverted
    TH2D* hAbsDiffnumu2nue=(TH2D*)hhist2d[1]->Clone("hAbsDiffnumu2nue");
    TH2D* hRelDiffnumu2nue=(TH2D*)hhist2d[1]->Clone("hRelDiffnumu2nue");
    TH2D* hOctResnumu2nue=(TH2D*)hhist2d[1]->Clone("hOctResnumu2nue");
    TH2D* hOctResnumu2nueLL=(TH2D*)hhist2d[1]->Clone("hOctResnumu2nueLL");
    
    TH2D* hAbsDiffnumu2nuebar=(TH2D*)hhist2d[1]->Clone("hAbsDiffnumu2nuebar");
    TH2D* hRelDiffnumu2nuebar=(TH2D*)hhist2d[1]->Clone("hRelDiffnumu2nuebar");
    TH2D* hOctResnumu2nuebar=(TH2D*)hhist2d[1]->Clone("hOctResnumu2nuebar");
    TH2D* hOctResnumu2nuebarLL=(TH2D*)hhist2d[1]->Clone("hOctResnumu2nuebarLL");
    
    TH2D* hOctResnumu2nueDiff=(TH2D*)hhist2d[1]->Clone("hOctResnumu2nueDiff");
    TH2D* hOctResnumu2nueDiffLL=(TH2D*)hhist2d[1]->Clone("hOctResnumu2nueDiffLL");
    
    //hRelDiffnumu2nue->Add(hhist2d[4],-1.0);
    //cout<<"Bin at 0.5 "<<hAbsDiffnumu2nuebar->GetXaxis()->FindBin(0.5)<<endl;
    Int_t th23maximalBin = hAbsDiffnumu2nuebar->GetXaxis()->FindBin(0.5);
    
    Double_t tmpOscDiffOct;
    Double_t minOscDiffOct;
    Double_t minSqrOscDiffOsc;
    Double_t tmpOscDiffOctbar;
    Double_t minOscDiffOctbar;
    Double_t minSqrOscDiffOscbar;
    
    Double_t tmpOscDiffOctLL;
    Double_t minOscDiffOctLL;
    Double_t tmpOscDiffOctbarLL;
    Double_t minOscDiffOctbarLL;
    
    Double_t tmpexp;
    Double_t tmpobs;
    for (int ix = 1; ix <= NBINth23; ix++){
        
        for (int iy = 1; iy <= NBINdcp; iy++){
            
            double binval_NH =hhist2d[1]->GetBinContent(ix,iy);
            double binval_IH =hhist2d[4]->GetBinContent(ix,iy);
            double binAbsDiff = binval_NH - binval_IH;
            double binRelDiff = (binval_NH-binval_IH)/binval_NH;
            hAbsDiffnumu2nue->SetBinContent(ix,iy,binAbsDiff);
            hRelDiffnumu2nue->SetBinContent(ix,iy,binRelDiff);
            
            binval_NH =hhist2d[3]->GetBinContent(ix,iy);
            binval_IH =hhist2d[5]->GetBinContent(ix,iy);
            binAbsDiff = binval_NH - binval_IH;
            binRelDiff = (binval_NH-binval_IH)/binval_NH;
            hAbsDiffnumu2nuebar->SetBinContent(ix,iy,binAbsDiff);
            hRelDiffnumu2nuebar->SetBinContent(ix,iy,binRelDiff);
            
            if (ix ==th23maximalBin) {
                hOctResnumu2nue->SetBinContent(ix,iy,0.0);
                hOctResnumu2nuebar->SetBinContent(ix,iy,0.0);
                hOctResnumu2nueDiff->SetBinContent(ix,iy,0.0);
                
                hOctResnumu2nueLL->SetBinContent(ix,iy,0.0);
                hOctResnumu2nuebarLL->SetBinContent(ix,iy,0.0);
                hOctResnumu2nueDiffLL->SetBinContent(ix,iy,0.0);
            }
            else if (ix<th23maximalBin){
                minOscDiffOct=1.0;
                minOscDiffOctbar = 1.0;
                minOscDiffOctLL = 1.0;
                minOscDiffOctbarLL = 1.0;
                for (int ixoct = th23maximalBin+1; ixoct<=NBINth23; ixoct++) {
                    //same iy
                    tmpexp = hhist2d[1]->GetBinContent(ix,iy);
                    tmpobs = hhist2d[1]->GetBinContent(ixoct,iy);
                    tmpOscDiffOct =TMath::Abs(tmpexp-tmpobs);
                    if (tmpOscDiffOct<minOscDiffOct) {
                        minOscDiffOct=tmpOscDiffOct;
                    }
                    //likelihood
                    tmpOscDiffOctLL = 2*(tmpexp-tmpobs+tmpobs*log(tmpobs/tmpexp));
                    if (tmpOscDiffOctLL<minOscDiffOctLL) {
                        minOscDiffOctLL = tmpOscDiffOctLL;
                    }
                    
                    //nuebar
                    tmpexp = hhist2d[3]->GetBinContent(ix,iy);
                    tmpobs = hhist2d[3]->GetBinContent(ixoct,iy);
                    tmpOscDiffOctbar =TMath::Abs(tmpexp-tmpobs);
                    if (tmpOscDiffOctbar<minOscDiffOctbar) {
                        minOscDiffOctbar=tmpOscDiffOctbar;
                    }
                    
                    //nuebar likelihood
                    tmpOscDiffOctbarLL = 2*(tmpexp-tmpobs+tmpobs*log(tmpobs/tmpexp));
                    if (tmpOscDiffOctbarLL<minOscDiffOctbarLL) {
                        minOscDiffOctbarLL = tmpOscDiffOctbarLL;
                    }
                    
                    
                    
                }//end ixoct
                minSqrOscDiffOsc = pow(minOscDiffOct,2)/pow(hhist2d[1]->GetBinContent(ix,iy),2);
                hOctResnumu2nue->SetBinContent(ix,iy,minSqrOscDiffOsc);
                hOctResnumu2nueLL->SetBinContent(ix,iy,minOscDiffOctLL);
                
                minSqrOscDiffOscbar = pow(minOscDiffOctbar,2)/pow(hhist2d[3]->GetBinContent(ix,iy),2);
                hOctResnumu2nuebar->SetBinContent(ix,iy,minSqrOscDiffOscbar);
                hOctResnumu2nuebarLL->SetBinContent(ix,iy,minOscDiffOctbarLL);
                
                hOctResnumu2nueDiff->SetBinContent(ix,iy,minSqrOscDiffOsc -minSqrOscDiffOscbar);
                hOctResnumu2nueDiffLL->SetBinContent(ix,iy,minOscDiffOctLL-minOscDiffOctbarLL);
                
                
            }
            else {
                minOscDiffOct=1.0;
                minOscDiffOctbar=1.0;
                minOscDiffOctLL = 1.0;
                minOscDiffOctbarLL = 1.0;
                for (int ixoctLo = 1; ixoctLo<=th23maximalBin-1; ixoctLo++) {
                    //same iy
                    tmpexp = hhist2d[1]->GetBinContent(ix,iy);
                    tmpobs = hhist2d[1]->GetBinContent(ixoctLo,iy);
                    tmpOscDiffOct =TMath::Abs(tmpexp-tmpobs);
                    if (tmpOscDiffOct<minOscDiffOct) {
                        minOscDiffOct=tmpOscDiffOct;
                    }
                    
                    tmpOscDiffOctLL = 2*(tmpexp-tmpobs+tmpobs*log(tmpobs/tmpexp));
                    if (tmpOscDiffOctLL<minOscDiffOctLL) {
                        minOscDiffOctLL = tmpOscDiffOctLL;
                    }
                    //nuebar
                    tmpexp = hhist2d[3]->GetBinContent(ix,iy);
                    tmpobs = hhist2d[3]->GetBinContent(ixoctLo,iy);
                    tmpOscDiffOctbar =TMath::Abs(tmpexp-tmpobs);
                    if (tmpOscDiffOctbar<minOscDiffOctbar) {
                        minOscDiffOctbar=tmpOscDiffOctbar;
                    }
                    //nuebar likelihood
                    tmpOscDiffOctbarLL = 2*(tmpexp-tmpobs+tmpobs*log(tmpobs/tmpexp));
                    if (tmpOscDiffOctbarLL<minOscDiffOctbarLL) {
                        minOscDiffOctbarLL = tmpOscDiffOctbarLL;
                    }
                    
                }//end ixoct
                
                //Chisquare (e-o)*(e-o)/(e*e);
                //Likelihood 2*(e-o + o*log(o/e));
                minSqrOscDiffOsc = pow(minOscDiffOct,2)/pow(hhist2d[1]->GetBinContent(ix,iy),2);
                hOctResnumu2nue->SetBinContent(ix,iy,minSqrOscDiffOsc);
                hOctResnumu2nueLL->SetBinContent(ix,iy,minOscDiffOctLL);
                
                minSqrOscDiffOscbar = pow(minOscDiffOctbar,2)/pow(hhist2d[3]->GetBinContent(ix,iy),2);
                hOctResnumu2nuebar->SetBinContent(ix,iy,minSqrOscDiffOscbar);
                hOctResnumu2nuebarLL->SetBinContent(ix,iy,minOscDiffOctbarLL);
                
                hOctResnumu2nueDiff->SetBinContent(ix,iy,minSqrOscDiffOsc -minSqrOscDiffOscbar);
                hOctResnumu2nueDiffLL->SetBinContent(ix,iy,minOscDiffOctLL-minOscDiffOctbarLL);
            }
            
        }
        
    }
    
    TCanvas *c4 = new TCanvas();
    gStyle->SetOptStat(0);
    hAbsDiffnumu2nue->Draw("colz");
    hAbsDiffnumu2nue->GetZaxis()->SetRangeUser(-0.02,0.02);
    c4->Print(Form("plots/%s_numu2nue_NOSubIO.eps",subname.Data()));
    
    TCanvas *c4b = new TCanvas();
    gStyle->SetOptStat(0);
    hRelDiffnumu2nue->Draw("colz");
    hRelDiffnumu2nue->GetZaxis()->SetRangeUser(-0.4,0.4);
    c4b->Print(Form("plots/%s_numu2nue_NOSubIOdivNO.eps",subname.Data()));
    
    
    TCanvas *c5 = new TCanvas();
    gStyle->SetOptStat(0);
    hAbsDiffnumu2nuebar->Draw("colz");
    hAbsDiffnumu2nuebar->GetZaxis()->SetRangeUser(-0.02,0.02);
    c5->Print(Form("plots/%s_numu2nuebar_NOSubIO.eps",subname.Data()));
    
    TCanvas *c5b = new TCanvas();
    gStyle->SetOptStat(0);
    hRelDiffnumu2nuebar->Draw("colz");
    hRelDiffnumu2nuebar->GetZaxis()->SetRangeUser(-0.4,0.4);
    c5b->Print(Form("plots/%s_numu2nuebar_NOSubIOdivNO.eps",subname.Data()));
    
    TCanvas *c6 = new TCanvas();
    gStyle->SetOptStat(0);
    hOctResnumu2nue->Draw("colz");
    hOctResnumu2nue->GetZaxis()->SetRangeUser(0.0,0.1);
    c6->Print(Form("plots/%s_numu2nue_OctRes.eps",subname.Data()));
    
    TCanvas *c6b = new TCanvas();
    gStyle->SetOptStat(0);
    hOctResnumu2nuebar->Draw("colz");
    hOctResnumu2nuebar->GetZaxis()->SetRangeUser(0.0,0.1);
    c6b->Print(Form("plots/%s_numu2nuebar_OctRes.eps",subname.Data()));
    
    TCanvas *c6c = new TCanvas();
    gStyle->SetOptStat(0);
    //hOctResnumu2nueDiff->Draw("colz");
    double level4OctRes[3] = {-0.001, 0.0, 0.001};
    /*double xlatexnuebar[nlevels]={0.42, 0.48, 0.56};
     double ylatexnuebar[nlevels]={-1, -1, -1.1};*/
    double xlatex4OctRes[3]={0.44, 0.50, 0.56};
    double ylatex4OctRes[3]={-0.5, 0.2, 2.8};
    TLatex *platex4OctRes[3];
    for (Int_t i=0; i<nlevels; ++i) {
        platex4OctRes[i] = new TLatex(xlatex4OctRes[i],ylatex4OctRes[i],Form("P=%.3f",level4OctRes[i]));
        platex4OctRes[i]->SetNDC(kFALSE); // <- use NDC coordinate
        platex4OctRes[i]->SetTextSize(0.04);
        platex4OctRes[i]->SetTextColor(colorindex[i]);
        platex4OctRes[i]->SetTextAngle(0);
    }
    
    TH2D* hOctResnumu2nueDiff_clone=(TH2D*)hOctResnumu2nueDiff->Clone("hOctResnumu2nueDiff_clone");
    hOctResnumu2nueDiff_clone->SetContour(nlevels,level4OctRes);
    std::vector<TGraph*> gVect_hist_octres_0;
    std::vector<TGraph*> gVect_hist_octres_1;
    std::vector<TGraph*> gVect_hist_octres_2;
    c6c->cd();
    hOctResnumu2nueDiff_clone->Draw("cont list");
    c6c->Update();
    TObjArray *contours_hist_octres = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    gVect_hist_octres_0 = GetContourGraphs(contours_hist_octres, "white",  0);
    gVect_hist_octres_1 = GetContourGraphs(contours_hist_octres, "black",  1);
    gVect_hist_octres_2 = GetContourGraphs(contours_hist_octres, "blue",  2);
    
    hOctResnumu2nueDiff->Draw("colz");
    //hOctResnumu2nueDiff->GetZaxis()->SetRangeUser(-0.1,0.1);
    //if 0.5GeV, the variation is larger?
    hOctResnumu2nueDiff->GetZaxis()->SetRangeUser(-0.15,0.15);
    for(int j=0; j<gVect_hist_octres_0.size(); j++){gVect_hist_octres_0[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_octres_1.size(); j++){gVect_hist_octres_1[j]->Draw("Lsame");}
    for(int j=0; j<gVect_hist_octres_2.size(); j++){gVect_hist_octres_2[j]->Draw("Lsame");}
    
    TLatex *   text2kpreOctRes = new TLatex(xpos,ypos,"#chi^{2}(P(#nu_{#mu}#rightarrow #nu_{e}))-#chi^{2}(P(#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e})), Normal Ordering");
    setTLatexStyle(text2kpreOctRes);
    text2kpreOctRes->Draw();
    //text2kprecond->Draw();
    for (Int_t i=0; i<nlevels; ++i) {platex4OctRes[i]->Draw();}
    mTrue->Draw();
    c6c->Print(Form("plots/%s_numu2nueDiff_OctRes.eps",subname.Data()));
    
    if (isDrawLikelihood) {
        TCanvas *c7 = new TCanvas();
        gStyle->SetOptStat(0);
        hOctResnumu2nueLL->Draw("colz");
        hOctResnumu2nueLL->GetZaxis()->SetRangeUser(0.0,0.004);
        c7->Print(Form("plots/%s_numu2nue_OctRes_LL.eps",subname.Data()));
        
        TCanvas *c7b = new TCanvas();
        gStyle->SetOptStat(0);
        hOctResnumu2nuebarLL->Draw("colz");
        hOctResnumu2nuebarLL->GetZaxis()->SetRangeUser(0.0,0.004);
        c7b->Print(Form("plots/%s_numu2nuebar_OctRes_LL.eps",subname.Data()));
        
        TCanvas *c7c = new TCanvas();
        gStyle->SetOptStat(0);
        hOctResnumu2nueDiffLL->Draw("colz");
        hOctResnumu2nueDiffLL->GetZaxis()->SetRangeUser(-0.004,0.004);
        c7c->Print(Form("plots/%s_numu2nueDiff_OctRes_LL.eps",subname.Data()));
    }
    
    
    
    //numu2ne sin13 vs sin23
    if (isDrawSin13vsSin23) {
        TCanvas *c8 = new TCanvas();
        gStyle->SetOptStat(0);
        titleStyle2D(hhist2d[6]);
        
        double prob_numu2nue_sin13vssin23_clonesol =hhist2d[6]->GetBinContent(hhist2d[6]->FindBin(sinsqth23_clonesol,sinsqth13_clonesol));
        cout<<" Prob. numu2nue "<<prob_numu2nue_sin13vssin23_clonesol<<endl;
        
        double levels_numu2nue_sin13vssin23[3] = {prob_numu2nue_sin13vssin23_clonesol*0.9, prob_numu2nue_sin13vssin23_clonesol, prob_numu2nue_sin13vssin23_clonesol*1.1};
        double xlatexnue_sin13vssin23[nlevels]={0.5, 0.53, 0.58};
        double ylatexnue_sin13vssin23[nlevels]={0.021, 0.021, 0.022};
        TLatex *platexnue_sin13vssin23[nlevelsnumudis];
        for (Int_t i=0; i<nlevels; ++i) {
            platexnue_sin13vssin23[i] = new TLatex(xlatexnue_sin13vssin23[i],ylatexnue_sin13vssin23[i],Form("P=%.3f",levels_numu2nue_sin13vssin23[i]));
            platexnue_sin13vssin23[i]->SetNDC(kFALSE); // <- use NDC coordinate
            platexnue_sin13vssin23[i]->SetTextSize(0.04);
            platexnue_sin13vssin23[i]->SetTextColor(colorindex[i]);
            platexnue_sin13vssin23[i]->SetTextAngle(-45);
        }
        
        TH2D* hnumu2nue_sin13vssin23_clone=(TH2D*)hhist2d[6]->Clone("hnumu2nue_sin13vssin23_clone");
        hnumu2nue_sin13vssin23_clone->SetContour(nlevels,levels_numu2nue_sin13vssin23);
        std::vector<TGraph*> gVect_hist_numu2nue_sin13vssin23_0;
        std::vector<TGraph*> gVect_hist_numu2nue_sin13vssin23_1;
        std::vector<TGraph*> gVect_hist_numu2nue_sin13vssin23_2;
        c8->cd();
        hnumu2nue_sin13vssin23_clone->Draw("cont list");
        c8->Update();
        TObjArray *contours_hist_numu2nue_sin13vssin23 = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
        gVect_hist_numu2nue_sin13vssin23_0 = GetContourGraphs(contours_hist_numu2nue_sin13vssin23, "white",  0);
        gVect_hist_numu2nue_sin13vssin23_1 = GetContourGraphs(contours_hist_numu2nue_sin13vssin23, "black",  1);
        gVect_hist_numu2nue_sin13vssin23_2 = GetContourGraphs(contours_hist_numu2nue_sin13vssin23, "blue",  2);
        
        hhist2d[6]->Draw("colz");
        hhist2d[6]->GetZaxis()->SetRangeUser(0.01,0.08);
        for(int j=0; j<gVect_hist_numu2nue_sin13vssin23_0.size(); j++){gVect_hist_numu2nue_sin13vssin23_0[j]->Draw("Lsame");}
        for(int j=0; j<gVect_hist_numu2nue_sin13vssin23_1.size(); j++){gVect_hist_numu2nue_sin13vssin23_1[j]->Draw("Lsame");}
        for(int j=0; j<gVect_hist_numu2nue_sin13vssin23_2.size(); j++){gVect_hist_numu2nue_sin13vssin23_2[j]->Draw("Lsame");}
        
        
        text2kprenue->Draw();
        text2kprecond->Draw();
        for (Int_t i=0; i<nlevels; ++i) {platexnue_sin13vssin23[i]->Draw();}
        mTrue_sin13vssin23->Draw();
        //check the TGraph, sin^22theta_13 tan^2theta_23
        if (isDrawTheoryPrediction) {
            Double_t valDegeneracyTheoryth13 = 4*sinsqth13_clonesol*(1-sinsqth13_clonesol)*TMath::Sqrt(sinsqth23_clonesol)/(TMath::Sqrt(1-sinsqth23_clonesol));
            //Double_t valDegeneracyTheoryth13 = 4*sinsqth13_clonesol*(1-sinsqth13_clonesol)*(sinsqth23_clonesol)/((1-sinsqth23_clonesol));
            Double_t *xvalDegTheoth13 = new Double_t[NBINth23];
            Double_t *yvalDegTheoth13 = new Double_t[NBINth23];
            
            for (Int_t ibinth23=1; ibinth23<=NBINth23+1; ibinth23++) {
                xvalDegTheoth13[ibinth23-1] = hhist2d[6]->GetXaxis()->GetBinCenter(ibinth23);
                yvalDegTheoth13[ibinth23-1] = pow(TMath::Sin(TMath::ASin(TMath::Sqrt(valDegeneracyTheoryth13*TMath::Sqrt(1-xvalDegTheoth13[ibinth23-1])/TMath::Sqrt(xvalDegTheoth13[ibinth23-1])))/2.),2.);
                //yvalDegTheoth13[ibinth23-1] = pow(TMath::Sin(TMath::ASin(TMath::Sqrt(valDegeneracyTheoryth13*(1-xvalDegTheoth13[ibinth23-1])/(xvalDegTheoth13[ibinth23-1])))/2.),2.);
                cout<<"xval "<<xvalDegTheoth13[ibinth23-1]<<" yval "<<yvalDegTheoth13[ibinth23-1]<<endl;
            }
            TGraph *pgrDegeneracyTheoryth13 = new TGraph(NBINth23,xvalDegTheoth13,yvalDegTheoth13);
            pgrDegeneracyTheoryth13->SetLineColor(kRed);
            pgrDegeneracyTheoryth13->Draw("L same");
            c8->Print(Form("plots/%s_numu2nue_sin13vssin23_wtheoPred.eps",subname.Data()));
        }
        else c8->Print(Form("plots/%s_numu2nue_sin13vssin23.eps",subname.Data()));
        
        
        //numu2nebar sin13 vs sin23
        TCanvas *c8b = new TCanvas();
        gStyle->SetOptStat(0);
        titleStyle2D(hhist2d[7]);
        
        double prob_numu2nuebar_sin13vssin23_clonesol =hhist2d[7]->GetBinContent(hhist2d[7]->FindBin(sinsqth23_clonesol,sinsqth13_clonesol));
        cout<<" Prob. numu2nuebar "<<prob_numu2nuebar_sin13vssin23_clonesol<<endl;
        
        double levels_numu2nuebar_sin13vssin23[3] = {prob_numu2nuebar_sin13vssin23_clonesol*0.9, prob_numu2nuebar_sin13vssin23_clonesol, prob_numu2nuebar_sin13vssin23_clonesol*1.1};
        double xlatexnuebar_sin13vssin23[nlevels]={0.5, 0.53, 0.58 };
        double ylatexnuebar_sin13vssin23[nlevels]={0.021, 0.021, 0.022};
        TLatex *platexnuebar_sin13vssin23[nlevelsnumudis];
        for (Int_t i=0; i<nlevels; ++i) {
            platexnuebar_sin13vssin23[i] = new TLatex(xlatexnuebar_sin13vssin23[i],ylatexnuebar_sin13vssin23[i],Form("P=%.3f",levels_numu2nuebar_sin13vssin23[i]));
            platexnuebar_sin13vssin23[i]->SetNDC(kFALSE); // <- use NDC coordinate
            platexnuebar_sin13vssin23[i]->SetTextSize(0.04);
            platexnuebar_sin13vssin23[i]->SetTextColor(colorindex[i]);
            platexnuebar_sin13vssin23[i]->SetTextAngle(-45);
        }
        
        TH2D* hnumu2nuebar_sin13vssin23_clone=(TH2D*)hhist2d[7]->Clone("hnumu2nuebar_sin13vssin23_clone");
        hnumu2nuebar_sin13vssin23_clone->SetContour(nlevels,levels_numu2nuebar_sin13vssin23);
        std::vector<TGraph*> gVect_hist_numu2nuebar_sin13vssin23_0;
        std::vector<TGraph*> gVect_hist_numu2nuebar_sin13vssin23_1;
        std::vector<TGraph*> gVect_hist_numu2nuebar_sin13vssin23_2;
        c8b->cd();
        hnumu2nuebar_sin13vssin23_clone->Draw("cont list");
        c8b->Update();
        TObjArray *contours_hist_numu2nuebar_sin13vssin23 = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
        gVect_hist_numu2nuebar_sin13vssin23_0 = GetContourGraphs(contours_hist_numu2nuebar_sin13vssin23, "white",  0);
        gVect_hist_numu2nuebar_sin13vssin23_1 = GetContourGraphs(contours_hist_numu2nuebar_sin13vssin23, "black",  1);
        gVect_hist_numu2nuebar_sin13vssin23_2 = GetContourGraphs(contours_hist_numu2nuebar_sin13vssin23, "blue",  2);
        
        hhist2d[7]->Draw("colz");
        hhist2d[7]->GetZaxis()->SetRangeUser(0.01,0.08);
        for(int j=0; j<gVect_hist_numu2nuebar_sin13vssin23_0.size(); j++){gVect_hist_numu2nuebar_sin13vssin23_0[j]->Draw("Lsame");}
        for(int j=0; j<gVect_hist_numu2nuebar_sin13vssin23_1.size(); j++){gVect_hist_numu2nuebar_sin13vssin23_1[j]->Draw("Lsame");}
        for(int j=0; j<gVect_hist_numu2nuebar_sin13vssin23_2.size(); j++){gVect_hist_numu2nuebar_sin13vssin23_2[j]->Draw("Lsame");}
        
        
        text2kprenue->Draw();
        text2kprecond->Draw();
        for (Int_t i=0; i<nlevels; ++i) {platexnuebar_sin13vssin23[i]->Draw();}
        mTrue_sin13vssin23->Draw();
        c8b->Print(Form("plots/%s_numu2nuebar_sin13vssin23.eps",subname.Data()));
    }
    
    if (isDrawSin13vsdcp) {
        //numu2nue sin13 vs sin23
        TCanvas *c9 = new TCanvas();
        gStyle->SetOptStat(0);
        titleStyle2D(hhist2d[8]);
        
        double prob_numu2nue_sin13vsdcp_clonesol =hhist2d[8]->GetBinContent(hhist2d[8]->FindBin(sinsqth13_clonesol,dcp_clonesol));
        cout<<" Prob. numu2nue "<<prob_numu2nue_sin13vsdcp_clonesol<<endl;
        
        double levels_numu2nue_sin13vsdcp[3] = {prob_numu2nue_sin13vsdcp_clonesol*0.9, prob_numu2nue_sin13vsdcp_clonesol, prob_numu2nue_sin13vsdcp_clonesol*1.1};
        double xlatexnue_sin13vsdcp[nlevels]={0.0196, 0.0222, 0.024};
        double ylatexnue_sin13vsdcp[nlevels]={-2, -2.1, -2};
        TLatex *platexnue_sin13vsdcp[nlevelsnumudis];
        for (Int_t i=0; i<nlevels; ++i) {
            platexnue_sin13vsdcp[i] = new TLatex(xlatexnue_sin13vsdcp[i],ylatexnue_sin13vsdcp[i],Form("P=%.3f",levels_numu2nue_sin13vsdcp[i]));
            platexnue_sin13vsdcp[i]->SetNDC(kFALSE); // <- use NDC coordinate
            platexnue_sin13vsdcp[i]->SetTextSize(0.04);
            platexnue_sin13vsdcp[i]->SetTextColor(colorindex[i]);
            platexnue_sin13vsdcp[i]->SetTextAngle(-10);
        }
        
        TH2D* hnumu2nue_sin13vsdcp_clone=(TH2D*)hhist2d[8]->Clone("hnumu2nue_sin13vsdcp_clone");
        hnumu2nue_sin13vsdcp_clone->SetContour(nlevels,levels_numu2nue_sin13vsdcp);
        std::vector<TGraph*> gVect_hist_numu2nue_sin13vsdcp_0;
        std::vector<TGraph*> gVect_hist_numu2nue_sin13vsdcp_1;
        std::vector<TGraph*> gVect_hist_numu2nue_sin13vsdcp_2;
        c9->cd();
        hnumu2nue_sin13vsdcp_clone->Draw("cont list");
        c9->Update();
        TObjArray *contours_hist_numu2nue_sin13vsdcp = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
        gVect_hist_numu2nue_sin13vsdcp_0 = GetContourGraphs(contours_hist_numu2nue_sin13vsdcp, "white",  0);
        gVect_hist_numu2nue_sin13vsdcp_1 = GetContourGraphs(contours_hist_numu2nue_sin13vsdcp, "black",  1);
        gVect_hist_numu2nue_sin13vsdcp_2 = GetContourGraphs(contours_hist_numu2nue_sin13vsdcp, "blue",  2);
        
        hhist2d[8]->Draw("colz");
        hhist2d[8]->GetZaxis()->SetRangeUser(0.01,0.08);
        for(int j=0; j<gVect_hist_numu2nue_sin13vsdcp_0.size(); j++){gVect_hist_numu2nue_sin13vsdcp_0[j]->Draw("Lsame");}
        for(int j=0; j<gVect_hist_numu2nue_sin13vsdcp_1.size(); j++){gVect_hist_numu2nue_sin13vsdcp_1[j]->Draw("Lsame");}
        for(int j=0; j<gVect_hist_numu2nue_sin13vsdcp_2.size(); j++){gVect_hist_numu2nue_sin13vsdcp_2[j]->Draw("Lsame");}
        
        
        text2kprenue->Draw();
        text2kprecond->Draw();
        for (Int_t i=0; i<nlevels; ++i) {platexnue_sin13vsdcp[i]->Draw();}
        mTrue_sin13vsdcp->Draw();
        //check the TGraph, sin^22theta_13 tan^2theta_23
        if (isDrawTheoryPrediction) {
            // Double_t valDegeneracyTheoryth13dcp = TMath::Sin(dcp_clonesol)*2*TMath::Sqrt(sinsqth13_clonesol)*TMath::Sqrt(1-sinsqth13_clonesol);
            Double_t valDegeneracyTheoryth13dcp = TMath::Sin(dcp_clonesol)*pow(sinsqth13_clonesol,5);
            //Double_t valDegeneracyTheoryth13dcp = sinsqth13_clonesol/TMath::Sin(dcp_clonesol);
            
            Double_t *xvalDegTheoth13dcp = new Double_t[NBINth13];
            Double_t *yvalDegTheoth13dcp = new Double_t[NBINth13];
            
            for (Int_t ibinth13=1; ibinth13<=NBINth13+1; ibinth13++) {
                xvalDegTheoth13dcp[ibinth13-1] = hhist2d[8]->GetXaxis()->GetBinCenter(ibinth13);
                //yvalDegTheoth13dcp[ibinth13-1] = TMath::ASin(valDegeneracyTheoryth13dcp/(2*TMath::Sqrt(xvalDegTheoth13dcp[ibinth13-1])*TMath::Sqrt(1-xvalDegTheoth13dcp[ibinth13-1])));
                
                yvalDegTheoth13dcp[ibinth13-1] = TMath::ASin(valDegeneracyTheoryth13dcp/pow(xvalDegTheoth13dcp[ibinth13-1],5));
                //yvalDegTheoth13dcp[ibinth13-1] = TMath::ASin((xvalDegTheoth13dcp[ibinth13-1])/valDegeneracyTheoryth13dcp);
                cout<<"xval "<<xvalDegTheoth13dcp[ibinth13-1]<<" yval "<<yvalDegTheoth13dcp[ibinth13-1]<<endl;
            }
            TGraph *pgrDegeneracyTheoryth13dcp = new TGraph(NBINth13,xvalDegTheoth13dcp,yvalDegTheoth13dcp);
            pgrDegeneracyTheoryth13dcp->SetLineColor(kRed);
            pgrDegeneracyTheoryth13dcp->Draw("L same");
            c9->Print(Form("plots/%s_numu2nue_sin13vsdcp_wtheoPred.eps",subname.Data()));
        }
        else c9->Print(Form("plots/%s_numu2nue_sin13vsdcp.eps",subname.Data()));
        
        
        //numu2nebar sin13 vs sin23
        TCanvas *c9b = new TCanvas();
        gStyle->SetOptStat(0);
        titleStyle2D(hhist2d[9]);
        
        double prob_numu2nuebar_sin13vsdcp_clonesol =hhist2d[9]->GetBinContent(hhist2d[9]->FindBin(sinsqth13_clonesol,dcp_clonesol));
        cout<<" Prob. numu2nuebar "<<prob_numu2nuebar_sin13vsdcp_clonesol<<endl;
        
        double levels_numu2nuebar_sin13vsdcp[3] = {prob_numu2nuebar_sin13vsdcp_clonesol*0.9, prob_numu2nuebar_sin13vsdcp_clonesol, prob_numu2nuebar_sin13vsdcp_clonesol*1.1};
        double xlatexnuebar_sin13vsdcp[nlevels]={0.0196, 0.0222, 0.024};
        double ylatexnuebar_sin13vsdcp[nlevels]={-2, -2.1, -2};
        TLatex *platexnuebar_sin13vsdcp[nlevelsnumudis];
        for (Int_t i=0; i<nlevels; ++i) {
            platexnuebar_sin13vsdcp[i] = new TLatex(xlatexnuebar_sin13vsdcp[i],ylatexnuebar_sin13vsdcp[i],Form("P=%.3f",levels_numu2nuebar_sin13vsdcp[i]));
            platexnuebar_sin13vsdcp[i]->SetNDC(kFALSE); // <- use NDC coordinate
            platexnuebar_sin13vsdcp[i]->SetTextSize(0.04);
            platexnuebar_sin13vsdcp[i]->SetTextColor(colorindex[i]);
            platexnuebar_sin13vsdcp[i]->SetTextAngle(-10);
        }
        
        TH2D* hnumu2nuebar_sin13vsdcp_clone=(TH2D*)hhist2d[9]->Clone("hnumu2nuebar_sin13vsdcp_clone");
        hnumu2nuebar_sin13vsdcp_clone->SetContour(nlevels,levels_numu2nuebar_sin13vsdcp);
        std::vector<TGraph*> gVect_hist_numu2nuebar_sin13vsdcp_0;
        std::vector<TGraph*> gVect_hist_numu2nuebar_sin13vsdcp_1;
        std::vector<TGraph*> gVect_hist_numu2nuebar_sin13vsdcp_2;
        c9b->cd();
        hnumu2nuebar_sin13vsdcp_clone->Draw("cont list");
        c9b->Update();
        TObjArray *contours_hist_numu2nuebar_sin13vsdcp = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
        gVect_hist_numu2nuebar_sin13vsdcp_0 = GetContourGraphs(contours_hist_numu2nuebar_sin13vsdcp, "white",  0);
        gVect_hist_numu2nuebar_sin13vsdcp_1 = GetContourGraphs(contours_hist_numu2nuebar_sin13vsdcp, "black",  1);
        gVect_hist_numu2nuebar_sin13vsdcp_2 = GetContourGraphs(contours_hist_numu2nuebar_sin13vsdcp, "blue",  2);
        
        hhist2d[9]->Draw("colz");
        hhist2d[9]->GetZaxis()->SetRangeUser(0.01,0.08);
        for(int j=0; j<gVect_hist_numu2nuebar_sin13vsdcp_0.size(); j++){gVect_hist_numu2nuebar_sin13vsdcp_0[j]->Draw("Lsame");}
        for(int j=0; j<gVect_hist_numu2nuebar_sin13vsdcp_1.size(); j++){gVect_hist_numu2nuebar_sin13vsdcp_1[j]->Draw("Lsame");}
        for(int j=0; j<gVect_hist_numu2nuebar_sin13vsdcp_2.size(); j++){gVect_hist_numu2nuebar_sin13vsdcp_2[j]->Draw("Lsame");}
        
        
        text2kprenue->Draw();
        text2kprecond->Draw();
        for (Int_t i=0; i<nlevels; ++i) {platexnuebar_sin13vsdcp[i]->Draw();}
        mTrue_sin13vsdcp->Draw();
        c9b->Print(Form("plots/%s_numu2nuebar_sin13vsdcp.eps",subname.Data()));
    }
    
    if (isDraw3DMarginalization) {
        TH2D* h3dMargOctResnumu2nue=(TH2D*)hhist2d[1]->Clone("h3dMargOctResnumu2nue");
        TH2D* h3dMargOctResnumu2nuebar=(TH2D*)hhist2d[1]->Clone("h3dMargOctResnumu2nuebar");
        TH2D* h3dMargOctResnumu2nueDiff=(TH2D*)hhist2d[1]->Clone("h3dMargOctResnumu2nue");
        
        Double_t tmpexp3d;
        Double_t tmpobs3d;
        Double_t tmpOscDiffOct3d;
        Double_t tmpOscDiffOct3dbar;
        Double_t minOscDiffOct3d;
        Double_t minOscDiffOct3dbar;
        Double_t minSqrOscDiffOsc3d;
        Double_t minSqrOscDiffOsc3dbar;
        
        for (int ix = 1; ix <= NBINth23; ix++){
            for (int iy = 1; iy <= NBINdcp; iy++){
                if (ix ==th23maximalBin) {
                    h3dMargOctResnumu2nue->SetBinContent(ix,iy,0.0);
                    h3dMargOctResnumu2nuebar->SetBinContent(ix,iy,0.0);
                    h3dMargOctResnumu2nueDiff->SetBinContent(ix,iy,0.0);
                }
                else if (ix<th23maximalBin){
                    minOscDiffOct3d=1.0;
                    minOscDiffOct3dbar = 1.0;
                    for (int ixoct = th23maximalBin+1; ixoct<=NBINth23; ixoct++) {
                        for (int iz = 1; iz<=NBINth13; iz++) {
                            tmpexp3d = hhist2d[1]->GetBinContent(ix,iy);
                            tmpobs3d = hhist3d[0]->GetBinContent(ixoct,iy,iz);
                            tmpOscDiffOct3d =TMath::Abs(tmpexp3d-tmpobs3d);
                            if (tmpOscDiffOct3d<minOscDiffOct3d ) {
                                minOscDiffOct3d=tmpOscDiffOct3d;
                            }
                            //nuebar
                            tmpexp3d = hhist2d[3]->GetBinContent(ix,iy);
                            tmpobs3d = hhist3d[1]->GetBinContent(ixoct,iy,iz);
                            tmpOscDiffOct3dbar =TMath::Abs(tmpexp3d-tmpobs3d);
                            if (tmpOscDiffOct3dbar<minOscDiffOct3dbar) {
                                minOscDiffOct3dbar=tmpOscDiffOct3dbar;
                            }
                        }//end iz
                    }//end ixoct
                    minSqrOscDiffOsc3d = pow(minOscDiffOct3d,2)/pow(hhist2d[1]->GetBinContent(ix,iy),2);
                    h3dMargOctResnumu2nue->SetBinContent(ix,iy,minSqrOscDiffOsc3d);
                    minSqrOscDiffOsc3dbar = pow(minOscDiffOct3dbar,2)/pow(hhist2d[3]->GetBinContent(ix,iy),2);
                    h3dMargOctResnumu2nuebar->SetBinContent(ix,iy,minSqrOscDiffOsc3dbar);
                    h3dMargOctResnumu2nueDiff->SetBinContent(ix,iy,minSqrOscDiffOsc3d -minSqrOscDiffOsc3dbar);
                }
                else {
                    minOscDiffOct3d=1.0;
                    minOscDiffOct3dbar=1.0;
                    for (int ixoctLo = 1; ixoctLo<=th23maximalBin-1; ixoctLo++) {
                        for (int izLo = 1; izLo<=NBINth13; izLo++) {
                            //same iy
                            tmpexp3d = hhist2d[1]->GetBinContent(ix,iy);
                            tmpobs3d = hhist3d[0]->GetBinContent(ixoctLo,iy,izLo);
                            tmpOscDiffOct3d =TMath::Abs(tmpexp3d-tmpobs3d);
                            if (tmpOscDiffOct3d<minOscDiffOct3d) {
                                minOscDiffOct3d=tmpOscDiffOct3d;
                            }
                            //nuebar
                            tmpexp3d = hhist2d[3]->GetBinContent(ix,iy);
                            tmpobs3d = hhist3d[1]->GetBinContent(ixoctLo,iy,izLo);
                            tmpOscDiffOct3dbar =TMath::Abs(tmpexp3d-tmpobs3d);
                            if (tmpOscDiffOct3dbar<minOscDiffOct3dbar) {
                                minOscDiffOct3dbar=tmpOscDiffOct3dbar;
                            }
                        }
                    }//end ixoct
                    minSqrOscDiffOsc3d = pow(minOscDiffOct3d,2)/pow(hhist2d[1]->GetBinContent(ix,iy),2);
                    h3dMargOctResnumu2nue->SetBinContent(ix,iy,minSqrOscDiffOsc3d);
                    minSqrOscDiffOsc3dbar = pow(minOscDiffOct3dbar,2)/pow(hhist2d[3]->GetBinContent(ix,iy),2);
                    h3dMargOctResnumu2nuebar->SetBinContent(ix,iy,minSqrOscDiffOsc3dbar);
                    h3dMargOctResnumu2nueDiff->SetBinContent(ix,iy,minSqrOscDiffOsc3d -minSqrOscDiffOsc3dbar);
                }//end else
            }//end for iy
        }//end for ix
        
        TCanvas *c10 = new TCanvas();
        gStyle->SetOptStat(0);
        h3dMargOctResnumu2nue->Draw("colz");
        h3dMargOctResnumu2nue->GetZaxis()->SetRangeUser(0.0,0.1);
        c10->Print(Form("plots/%s_numu2nue_OctRes_3dmarg.eps",subname.Data()));
        
        TCanvas *c10b = new TCanvas();
        gStyle->SetOptStat(0);
        h3dMargOctResnumu2nuebar->Draw("colz");
        h3dMargOctResnumu2nuebar->GetZaxis()->SetRangeUser(0.0,0.1);
        c10b->Print(Form("plots/%s_numu2nuebar_OctRes_3dmarg.eps",subname.Data()));
        
        TCanvas *c10c = new TCanvas();
        gStyle->SetOptStat(0);
        h3dMargOctResnumu2nueDiff->Draw("colz");
        h3dMargOctResnumu2nueDiff->GetZaxis()->SetRangeUser(-0.1,0.1);
        c10c->Print(Form("plots/%s_numu2nueDiff_OctRes_3dmarg.eps",subname.Data()));
        
        //check difference between marginalization on th13 and no-marginalization
        
        TH2D *h3dMargOctResnumu2nueDiffDC = (TH2D*)h3dMargOctResnumu2nueDiff->Clone("h3dMargOctResnumu2nueDiffDC");
        for (int ix = 1; ix <= NBINth23; ix++){
            for (int iy = 1; iy <= NBINdcp; iy++){
                double tmpDiff = hOctResnumu2nueDiff->GetBinContent(ix,iy) - h3dMargOctResnumu2nueDiff->GetBinContent(ix,iy);
                h3dMargOctResnumu2nueDiffDC->SetBinContent(ix,iy,tmpDiff);
                
            }
        }
        
        TCanvas *c10d = new TCanvas();
        gStyle->SetOptStat(0);
        h3dMargOctResnumu2nueDiffDC->Draw("colz");
        h3dMargOctResnumu2nueDiffDC->GetZaxis()->SetRangeUser(-0.01,0);
        c10d->Print(Form("plots/%s_numu2nueDiff_OctRes_3dmargDC.eps",subname.Data()));
        
    }//end if drawing marginalization
    
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
