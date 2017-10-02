//define simple oscillation prob. with four parameter
//Two-flavor simple calculation
//1. Neutrino energy
//2. Baseline
//3. Mixing angle theta_23
//4. Mass square splitting Dmsq23
//howtorun: root plotMu2Mu.C
double model2flav_vac_MuSurvive(double Ene, double baseL, double theta23, double dmsq_23);

void plotMu2Mu(){
    double myEne = 0.6;//T2K peak energy in GeV
    double myDmsq23 = 2.42e-3;//in eV
    double myBaseL = 295;//T2K km
    double myTheta23 = 0.785; //in radian 45*Pi/180
    
    //test oscillation probabity
    double myProbmu2mu = model2flav_vac_MuSurvive(myEne,myBaseL,myTheta23, myDmsq23);
    cout<<"probably at the peak"<<myProbmu2mu<<endl;
    //set number of steps you want to calculate
    const Int_t nEnStep = 200;
    // set maximal energy to plot
    double EnergyMax = 5.0;
    
    //Graph of oscillation prob. as function of energy
    //Note:formula is not work at zero energy
    TGraph *pGraphMu2Mu;
    Double_t pEne[nEnStep], pProb[nEnStep];

    for (Int_t istep=1; istep<=nEnStep; ++istep) {
        //energy value at each step
        pEne[istep-1] = EnergyMax*istep/nEnStep;
        //calculated osc. prob. at each step
        pProb[istep-1] = model2flav_vac_MuSurvive(pEne[istep-1],myBaseL,myTheta23, myDmsq23);
        cout<<"step "<<istep<<" Energy "<<pEne[istep-1]<<" Prob. "<<pProb[istep-1]<<endl;
    }
    //create new Graph based on the array
    pGraphMu2Mu = new TGraph(nEnStep,pEne,pProb);
    //set the width for graph line
    pGraphMu2Mu->SetLineWidth(2);
    
    //create new Canvas
    new TCanvas;
    pGraphMu2Mu->Draw("APL");//draw graph
    pGraphMu2Mu->SetTitle("");//title of graph can be ignore
    //set title for the graph
    pGraphMu2Mu->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    pGraphMu2Mu->GetYaxis()->SetTitle("Survival Probability #nu_{#mu}#rightarrow #nu_{#mu}");
    //save the graph
    gPad->Print("prob_2flav_mu2mu_t2k.eps");//eps format
    gPad->Print("prob_2flav_mu2mu_t2k.pdf");//pdf format
}
// formula for two-flavor oscillation in vacuum
double model2flav_vac_MuSurvive(double Ene, double baseL, double theta_23, double dmsq_23){
    double sin_BL = sin(1.267*baseL*dmsq_23/Ene);//L/E term
    double sinsq_2th23 = pow(sin(2*theta_23),2);//mixing term
    double probmu2mu = 1 - sinsq_2th23*sin_BL*sin_BL;//survival prob.
    
    return probmu2mu;
}
