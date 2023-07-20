//
//  sensitivity_ccqe_numudisapearance.C
////////////////////////////////////////////////////
//
//  To make sensitivity on atmospheric parameter with simple numu disappearance
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
double LogLikelihood(TH1* exp, TH1* obs)
{
    double llh_val = 0;
    // Hardcoded number of bins...
    for(int n = 1; n <= exp->GetNbinsX(); ++n){
        const double e_binval = exp->GetBinContent(n);
        const double o_binval = obs->GetBinContent(n);
       // if(e || o) ret += -2*log(TMath::Poisson(o,e));
        llh_val += 2*(e_binval-o_binval);
        if (o_binval!=0) {
            llh_val += -2*o_binval*log(e_binval/o_binval);
        }
    }
    return llh_val;
}

double model3flav_vac_MuToElec(bool isNeutrino, double Ene, double baseL, double theta_12, double theta_13, double theta_23, double delta, double dmsq_32, double dmsq_21){
    double dmsq_31 = dmsq_32 + dmsq_21;
    double fConvert = 1.267;
    double LT = 4*pow(cos(theta_13),2)*pow(sin(theta_13),2)*pow(sin(theta_23),2)*pow(sin(fConvert*baseL*dmsq_31/Ene),2);
    double CPC = 8*pow(cos(theta_13),2)*sin(theta_12)*sin(theta_13)*sin(theta_23)*(cos(theta_12)*cos(theta_23)*cos(delta)-sin(theta_12)*sin(theta_13)*sin(theta_23))*cos(fConvert*baseL*dmsq_32/Ene)*sin(fConvert*baseL*dmsq_31/Ene)*sin(fConvert*baseL*dmsq_21/Ene);
    double CPV = 8*pow(cos(theta_13),2)*cos(theta_12)*cos(theta_23)*sin(theta_12)*sin(theta_13)*sin(theta_23)*sin(delta)*sin(fConvert*baseL*dmsq_32/Ene)*sin(fConvert*baseL*dmsq_31/Ene)*sin(fConvert*baseL*dmsq_21/Ene);

    double Solar = 4*pow(sin(theta_12),2)*pow(cos(theta_13),2)*(pow(cos(theta_12),2)*pow(cos(theta_23),2)+pow(sin(theta_12),2)*pow(sin(theta_23),2)*pow(sin(theta_13),2)-2*cos(theta_12)*cos(theta_23)*sin(theta_12)*sin(theta_13)*sin(theta_23)*cos(delta))*pow(sin(fConvert*baseL*dmsq_21/Ene),2);
    
    double probmu2e;
    if(isNeutrino) probmu2e= LT+CPC-CPV+Solar;//appearance.
    else probmu2e = LT+CPC+CPV+Solar;//for antineutrino
    return probmu2e;
}

// formula for two-flavor oscillation in vacuum
//input theta_23 as sin^2 2theta 23 for convenience
double model2flav_vac_MuSurvive(double Ene, double baseL, double theta_23, double dmsq_23){
    double sin_BL = sin(1.267*baseL*dmsq_23/Ene);//L/E term
    double sinsq_2th23 = theta_23;//pow(sin(2*theta_23),2);//mixing term
    double ProbMu2e = 1 - sinsq_2th23*sin_BL*sin_BL;//survival prob.
    
    return ProbMu2e;
}

void sensitivity_ccqe_nueappearance(){
    TFile *pfile = new TFile("../eventpred/op_eventrate.root","READ");
    
    const int NCHAN4SHOW = 6;
    const int NDET = 2;
    TH1D *heventrate4show[NDET][NCHAN4SHOW];
    for (Int_t ichan=0; ichan<NCHAN4SHOW; ++ichan) {
        heventrate4show[0][ichan] = (TH1D*)pfile->Get(Form("ND_eventrate_chan%d",ichan));
        heventrate4show[1][ichan] = (TH1D*)pfile->Get(Form("FD_eventrate_chan%d",ichan));
    }
    //Number of event for normalization
    double NEVENTPred_nom = heventrate4show[1][1] ->Integral(0,heventrate4show[1][1]->GetNbinsX()+1 );
    //here we use the FD event rate CCQE only
    const int NEVENTMC = 10000;
    //normalization factor
    double normfactor = NEVENTPred_nom*100/NEVENTMC*1.0;
    cout<<"normalization factor "<<normfactor<<endl;
    
    TFile *popfile = new TFile("op_eventmc_4sensitivity_app.root","RECREATE");
    TTree *ptree = new TTree("Event","Neutrino event");
    double Enu, ProbMu2e;
    
    ptree->Branch("enu",&Enu,"enu/D");
    //ptree->Branch("ProbMu2e",&ProbMu2e,"ProbMu2e/D");
    
    const double T2KBASELINE = 295;//km
    
    double myEne = 0.6;//T2K peak energy in GeV
    double myBaseL = 295;//T2K km
    double myTheta12 = 0.583;
    double myTheta13 = 0.149;
    double myTheta23 = 0.785; //in radian degree*Pi/180
    double mydelta = 0;//-1.0*TMath::Pi()/2.;//-1.57;//in radian degree*Pi/180
    double myDmsq21 = 7.41e-5;
    //double myDmsq31 = 2.514e-3;
    double myDmsq32 = 2.511e-3-myDmsq21;//in eV
    double myDmsq31 = myDmsq32 + myDmsq21;
    
    
    const int NDCPBIN = 51;//50+1
    const double dcp_MIN = -TMath::Pi();
    const double dcp_MAX = TMath::Pi();
    const double dcp_BEST = -TMath::Pi()/2.;//45*Pi/180
    
    const int NTH13BIN = 51;//50+1
    const double th13_MIN = 0.020;
    const double th13_MAX = 0.024;
    const double th13_BEST = 0.022;
    
    TH1D *hpred_osc[NTH13BIN][NDCPBIN];
    for (int ith13=0; ith13<NTH13BIN; ++ith13) {
        for (int idcp=0; idcp<NDCPBIN; ++idcp) {
            hpred_osc[ith13][idcp] = new TH1D(Form("hpred_osc_dmbin%d_thetabin%d",ith13,idcp),"",100,0,5);
        }
    }
    TH1D *hpred_noosc = new TH1D("hpred_noosc","",100,0,5);
    TH1D *hdata = new TH1D("hdata","",100,0,5);
    
    for (int ievmc=0; ievmc<NEVENTMC; ++ievmc) {
        if(ievmc%1000==0) cout<<"processing "<<ievmc<<"/"<<NEVENTMC<<endl;
        Enu = heventrate4show[1][1] -> GetRandom();
        ptree->Fill();
        ProbMu2e  = 0.0001;
        hpred_noosc->Fill(Enu,ProbMu2e);
        
        for (int ith13=0; ith13<NTH13BIN; ++ith13) {
            double th13_val = th13_MIN + (th13_MAX-th13_MIN)*(ith13+0.5)/(NTH13BIN-1);
            for (int idcp=0; idcp<NDCPBIN; ++idcp) {
                double dcp_val = dcp_MIN + (dcp_MAX-dcp_MIN)*(idcp+0.5)/(NDCPBIN-1);
                //ProbMu2e = model2flav_vac_MuSurvive(Enu,T2KBASELINE,dcp_val,th13_val);
                ProbMu2e = model3flav_vac_MuToElec(true,Enu,myBaseL,myTheta23 ,th13_val,myTheta12, dcp_val, myDmsq32, myDmsq21);
                hpred_osc[ith13][idcp]->Fill(Enu,ProbMu2e);
            }
        }
        
        //best fit
        ProbMu2e = model3flav_vac_MuToElec(true,Enu,myBaseL,myTheta23 ,th13_BEST,myTheta12, dcp_BEST, myDmsq32, myDmsq21);//model2flav_vac_MuSurvive(Enu,T2KBASELINE,dcp_BEST,th13_BEST);
        hdata->Fill(Enu, ProbMu2e);
        
    }
    //take normalization
    hdata->Scale(normfactor);
    hpred_noosc->Scale(normfactor);
    for (int ith13=0; ith13<NTH13BIN; ++ith13) {
        for (int idcp=0; idcp<NDCPBIN; ++idcp) {
            hpred_osc[ith13][idcp]->Scale(normfactor);
            hpred_osc[ith13][idcp]->Write(Form("hpred_osc_th13bin%d_dcpbin%d",ith13,idcp));
        }
    }
    
    hdata->Write("hdata");
    hpred_noosc->Write("hpred_noosc");
    
    //surface of chisq
    TH2D *hsurface = new TH2D("hsurface","",NDCPBIN-1,dcp_MIN,dcp_MAX,NTH13BIN-1,th13_MIN,th13_MAX);
    for (int ith13=1; ith13<NTH13BIN; ++ith13) {
        for (int idcp=1; idcp<NDCPBIN; ++idcp) {
            double likelihood = LogLikelihood(hpred_osc[ith13-1][idcp-1],hdata);
            hsurface->SetBinContent(idcp,ith13,likelihood);
        }
    }
    hsurface->Write("hsurface");
    
    //plot the contour
    const int nlevels = 3;
    //http://www.reid.ai/2012/09/chi-squared-distribution-table-with.html
    double levels[nlevels] = {2.30, 4.61,6.18};//68%, 90% and 2 sigma
    //get minimum
    
    //double levels[nlevels] = {2.71, 6.63,9.0};//90%, 99%, 3 sigma

    
    //
    ptree->Write();
    popfile->Close();
    
    
    
}
