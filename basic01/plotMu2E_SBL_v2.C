//
//  plotMu2E_SBL.C
////////////////////////////////////////////////////
//
//  Simple oscillation probability
//
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu

//howtorun: root plotMu2E_SBL.C
//as function of L/E
double model2flav_vac_Mu2E_LE(double L2E, double sinsq2theta, double dmsq);
void plotMu2E_SBL_v2(){
    //number case to check
    const Int_t NCASE= 3;
    
    double myDmsq[NCASE] = {0.6,1.0,2.0};//in eV
    double mySinsq2Theta[NCASE] = {0.01,0.004,0.002};
    TString mylegend[NCASE];
    Int_t colorLine[NCASE] = {1,2,4};
    
    //set number of steps you want to calculate
    const Int_t nL2EStep = 200;
    // set maximal energy to plot
    double L2E_min = 0.2;
    double L2E_max = 2.0;
    
    //Graph of oscillation prob. as function of energy
    //Note:formula is not work at zero energy
    TGraph *pGraphMu2E_l2e[NCASE];
    double pL2E[nL2EStep], pProb_L2E[nL2EStep];
    TLegend* leg0 = new TLegend(.12, .58, 0.5, .82);
       leg0->SetFillStyle(0);
       leg0->SetBorderSize(0);
       leg0->SetTextSize(18);
       leg0->SetTextFont(43);
    for (Int_t icase=0; icase<NCASE; ++icase) {
        mylegend[icase]=Form("#Delta m^{2}=%.2f[eV^{2}], sin^{2}2#theta=%.3f",myDmsq[icase],mySinsq2Theta[icase]);
        for (Int_t istep=1; istep<=nL2EStep; ++istep) {
            //energy value at each step
            pL2E[istep-1] = 0.2+(2.0-0.2)*istep/nL2EStep;
            //calculated osc. prob. at each step
            pProb_L2E[istep-1] = model2flav_vac_Mu2E_LE(pL2E[istep-1],mySinsq2Theta[icase], myDmsq[icase]);
            
            cout<<"step "<<istep<<" Energy "<<pL2E[istep-1]<<" Prob. "<<pProb_L2E[istep-1]<<endl;
            
        }
        pGraphMu2E_l2e[icase] = new TGraph(nL2EStep,pL2E,pProb_L2E);
        pGraphMu2E_l2e[icase]->SetLineWidth(2);
        pGraphMu2E_l2e[icase]->SetLineColor(colorLine[icase]);
        leg0->AddEntry(pGraphMu2E_l2e[icase],mylegend[icase].Data(),"l");
    }
    
    
    //create new Canvas
    new TCanvas;
    pGraphMu2E_l2e[0]->Draw("APL");//draw graph
    pGraphMu2E_l2e[0]->SetTitle("");//title of graph can be ignore
    pGraphMu2E_l2e[0]->GetXaxis()->SetTitle("L/E [m/MeV]");
    pGraphMu2E_l2e[0]->GetYaxis()->SetTitle("Appearance Probability #nu_{#mu}#rightarrow #nu_{e}");
    for (Int_t icase=1; icase<NCASE; ++icase) {
        pGraphMu2E_l2e[icase]->Draw("L same");
    }
    leg0->Draw();
    //save the graph
    gPad->Print("prob_2flav_mu2e_SBL_l2e.eps");//eps format
    gPad->Print("prob_2flav_mu2e_SBL_l2e.pdf");//pdf format
    gPad->Print("prob_2flav_mu2e_SBL_l2e.png");//pdf format
    
}


double model2flav_vac_Mu2E_LE(double L2E, double sinsq2theta, double dmsq){
    double sin_BL = sin(1.267*L2E*dmsq);//L/E term
    double probmu2e =  sinsq2theta*sin_BL*sin_BL;//survival prob.
    return probmu2e;
}
