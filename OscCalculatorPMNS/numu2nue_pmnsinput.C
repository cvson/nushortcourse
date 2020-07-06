{

    //#include<TH1>
    //#include<TH2>
    //#include<TLatex>
    //#include<TString>
    #include <complex>
    
    gROOT->ProcessLine(".L OscCalculatorPMNS.C+");
  
    double x=0.0; //Energy
    
    Bool_t isNOVA=false;
    Bool_t isLBNE=false;//T2K if false
  
  
    //input parameters
    OscCalculatorPMNS myProb;
    double myDeltam2=2.5;//change from 2.42
    double myDeltam21=7.6;//change from 7.54
    double myTheta12=0.60126;// this is radia   TMath::ASin(TMath::Sqrt(0.8704))/2
    double myTheta23=0.7854;// TMath::ASin(TMath::Sqrt(1))/2
    //double myTheta23=0.75;
    double myTheta13=0.1479;//TMath::ASin(TMath::Sqrt(0.085))/2
    double myL=isNOVA? 810:295.;
    if(isLBNE) myL = 1250;
    double ent2k=isNOVA? 2.0:0.6;
    if(isLBNE) ent2k =2.0;
        
    
            
    //Set all parameters here
    myProb.SetL(myL);
  
    myProb.SetRho(2.71);
    myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
    myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));

    /*myProb.SetTh12(myTheta12);
    myProb.SetTh23(myTheta23);
    myProb.SetTh13(myTheta13);
    myProb.SetdCP(0.0001);*/

    //Flag to trigger antineutrino probabilities +1=neutrinos, -1=antineutrinos
    int isAnti=1;

    //Note, probability flavors are defined by pdg code:
    //nue  =12, antinue  =-12
    //numu =14, antinue  =-14
    //nutau=16, antinutau=-16
    //For instance, to compute P(numu->nue) you would call:
    // myProb.P(14,12,E)

    //Read back parameters to check
    std::cout << myProb.GetL() << endl
        << myProb.GetRho() << endl
	    /*<< myProb.GetDmsq21() << endl
	    << myProb.GetDmsq32() << endl  
	    << myProb.GetTh12() << endl
	    << myProb.GetTh23() << endl  
	    << myProb.GetTh13() << endl  
	    << myProb.GetdCP() << endl*/
	    << endl;

    std::complex<double> fPMNSInput[3][3];
    std::complex<double> idelta(0.0,0.0001);
    
    double s12 = sin(myTheta12);
    double s23 = sin(myTheta23);
    double s13 = sin(myTheta13);
    double c12 = cos(myTheta12);
    double c23 = cos(myTheta23);
    double c13 = cos(myTheta13);
    
    fPMNSInput[0][0] =  c12*c13;
    fPMNSInput[0][1] =  s12*c13;
    fPMNSInput[0][2] =  s13*exp(-idelta);
    
    fPMNSInput[1][0] = -s12*c23-c12*s23*s13*exp(idelta);
    fPMNSInput[1][1] =  c12*c23-s12*s23*s13*exp(idelta);
    fPMNSInput[1][2] =  s23*c13;
    
    fPMNSInput[2][0] =  s12*s23-c12*c23*s13*exp(idelta);
    fPMNSInput[2][1] = -c12*s23-s12*c23*s13*exp(idelta);
    fPMNSInput[2][2] =  c23*c13;
    
    
    myProb.SetMixPMNS(fPMNSInput);
    
    //This step is done automatically when computing probabilities,
    //but useful to do it here to check code works
    //myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
    myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());

    //Print mixing matrix and delta_m_sqrs.
    myProb.PrintMix();
    myProb.PrintDeltaMsqrs();

    //Keep this here for reference
    /*
     double par[9]={0.61,                       //Theta12
                   myTheta23,                  //Theta23                   
                   0.086,                       //Theta13
                   myDeltam2*TMath::Power(10,-3.), //DeltaM32^2
                   7.59*TMath::Power(10,-5.),  //DeltaM21^2
                   0.0001,                       //DeltaCP (3Pi/2)
		   2.71,                        //Rock Density  
                   myL,                        //L
		   1};                        //NuAntiNu
     */
    double epsilon=0.00001;
  

    TGraph *Mu2eDcp0=0;
    TGraph *Mu2eDcp270 =0;
    TGraph *AntiMu2eDcp0=0;
    TGraph *AntiMu2eDcp270=0;

    
    Int_t npoints=1000;
    Double_t loveremax =isNOVA?5.0:2.0;//energy range GeV
    
    Double_t musx[npoints];
    Double_t musynueDcp0[npoints];
    Double_t musynueDcp270[npoints];
    Double_t antimusynueDcp0[npoints];
    Double_t antimusynueDcp270[npoints];
    //for dcp = 0
    for (Int_t i=1; i < npoints+1; i++){
        
        x = i/(npoints*1.0/loveremax) + epsilon;//vary energy range
        musx[i-1]=x;
        
        myProb.Reset();
        myProb.SetL(myL);
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        musynueDcp0[i-1]=myProb.P(isAnti*14,isAnti*12,x);
        
    }
    
    //for dcp = 270
    for (Int_t i=1; i < npoints+1; i++){
        
        x = i/(npoints*1.0/loveremax) + epsilon;//vary energy range
        musx[i-1]=x;
        
        myProb.Reset();
        myProb.SetL(myL);
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(3*2.0*TMath::Pi()/4 + epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        musynueDcp270[i-1]=myProb.P(isAnti*14,isAnti*12,x);
        
    }
    
    //for anti-neutrino
    isAnti = -1;
    
    //for dcp = 0
    for (Int_t i=1; i < npoints+1; i++){
        
        x = i/(npoints*1.0/loveremax) + epsilon;//vary energy range
        musx[i-1]=x;
        
        myProb.Reset();
        myProb.SetL(myL);
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        antimusynueDcp0[i-1]=myProb.P(isAnti*14,isAnti*12,x);
        
    }
    
    //for dcp = 270
    for (Int_t i=1; i < npoints+1; i++){
        
        x = i/(npoints*1.0/loveremax) + epsilon;//vary energy range
        musx[i-1]=x;
        
        myProb.Reset();
        myProb.SetL(myL);
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(3*2.0*TMath::Pi()/4 + epsilon);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        antimusynueDcp270[i-1]=myProb.P(isAnti*14,isAnti*12,x);
        
    }
    
    
    
    
    Mu2eDcp0 = new TGraph(npoints,musx,musynueDcp0);
    Mu2eDcp270 = new TGraph(npoints,musx,musynueDcp270);
    AntiMu2eDcp0 = new TGraph(npoints,musx,antimusynueDcp0);
    AntiMu2eDcp270 = new TGraph(npoints,musx,antimusynueDcp270);
    
  
  
    //Plot everything

    TH2F *histo1=0;
    histo1 = new TH2F("NewOscProb","NewOscProb",npoints, 0., loveremax, 500, 0, .1);
  
    TCanvas *canvas1=0;
    canvas1 = new TCanvas("Can1","Can1",1000,800);
  
    canvas1->SetLeftMargin(canvas1->GetLeftMargin()*1.2);
    canvas1->SetBottomMargin(canvas1->GetBottomMargin()*1.2);
  
    
    canvas1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    //gPad->SetLeftMargin(gPad->GetLeftMargin()*1.4);
    double xcoord = 0.58;
    TString delMStr32=Form("%2.2f #times 10^{-3} eV^{2}",myDeltam2);
    TString legStr32="|#Delta m^{2}_{32}|=";
    TLatex* tlx=new TLatex(xcoord, 0.85,legStr32+delMStr32);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.03);
    tlx->SetTextAlign(12);
	
	TString delMStr21=Form("%2.2f #times 10^{-5} eV^{2}",myDeltam21);
    TString legStr21="|#Delta m^{2}_{21}|=";
    TLatex* tlx1=new TLatex(xcoord, 0.80,legStr21+delMStr21);
    tlx1->SetNDC(kTRUE); // <- use NDC coordinate
    tlx1->SetTextSize(0.03);
    tlx1->SetTextAlign(12);
	
	
    TString Th23Str=Form("sin^{2}2#theta_{23}=%2.2f",TMath::Sin(2*myTheta23)*TMath::Sin(2*myTheta23));
    TLatex* tlx2=new TLatex(xcoord, 0.75,Th23Str);
    tlx2->SetNDC(kTRUE); // <- use NDC coordinate
    tlx2->SetTextSize(0.03);
    tlx2->SetTextAlign(12);
	
	TString Th12Str=Form("sin^{2}2#theta_{12}=%2.2f",TMath::Sin(2*myTheta12)*TMath::Sin(2*myTheta12));
    TLatex* tlx3=new TLatex(xcoord, 0.70,Th12Str);
    tlx3->SetNDC(kTRUE); // <- use NDC coordinate
    tlx3->SetTextSize(0.03);
	tlx3->SetTextAlign(12);
    
	TString Th13Str=Form("sin^{2}2#theta_{13}=%2.2f",TMath::Sin(2*myTheta13)*TMath::Sin(2*myTheta13));
    TLatex* tlx4=new TLatex(xcoord, 0.65,Th13Str);
    tlx4->SetNDC(kTRUE); // <- use NDC coordinate
    tlx4->SetTextSize(0.03);
    tlx4->SetTextAlign(12);
	
	
    TString LStr=Form("L=%2.0f km",myL);
    TLatex* tlx5=new TLatex(xcoord, 0.60,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    tlx5->SetTextSize(0.03);
    tlx5->SetTextAlign(12);
  
  
    histo1->SetTitle("");
    
    histo1->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    histo1->GetYaxis()->SetTitle("Prob. (#nu_{#mu}#rightarrow#nu_{e} or #bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e}) ");
    histo1->GetYaxis()->CenterTitle();
    histo1->GetXaxis()->CenterTitle();
    histo1->GetXaxis()->SetLabelSize(histo1->GetXaxis()->GetTitleSize()*1.2);
    histo1->GetYaxis()->SetLabelSize(histo1->GetYaxis()->GetTitleSize()*1.2);
    histo1->GetXaxis()->SetTitleSize(histo1->GetXaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleSize(histo1->GetYaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    histo1->GetXaxis()->SetRangeUser(0.2,loveremax);
    histo1->Draw();
    Mu2eDcp0->SetLineWidth(2);
    Mu2eDcp270->SetLineWidth(2);
    AntiMu2eDcp0->SetLineWidth(2);
    AntiMu2eDcp270->SetLineWidth(2);
    
    Int_t ci;
    
    ci = TColor::GetColor("#0072B2");
    Mu2eDcp0->SetLineColor(ci);
    Mu2eDcp0->Draw("C same");
    
    ci = TColor::GetColor("#D55E00");
    Mu2eDcp270->SetLineColor(ci);
    Mu2eDcp270->Draw("C same");
    
    
    ci = TColor::GetColor("#0072B2");
    AntiMu2eDcp0->SetLineColor(ci);
    AntiMu2eDcp0->SetLineStyle(3);
    AntiMu2eDcp0->Draw("C same");
    
    ci = TColor::GetColor("#D55E00");
    AntiMu2eDcp270->SetLineColor(ci);
    AntiMu2eDcp270->SetLineStyle(3);
    AntiMu2eDcp270->Draw("C same");
    
  
    tlx->Draw("same");
    tlx1->Draw("same");
    tlx2->Draw("same");
    tlx3->Draw("same");
    tlx4->Draw("same");
    tlx5->Draw("same");
    Double_t xlegmin = 0.35;
    Double_t ylegmin = 0.65;
    TLegend* leg = new TLegend(xlegmin,ylegmin,xlegmin+0.3,ylegmin+0.23);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(18);
    leg->SetTextFont(43);;
    leg->AddEntry(Mu2eDcp0, "#nu, NH, #delta_{CP}=0","l");
    leg->AddEntry(Mu2eDcp270, "#nu,NH, #delta_{CP}=270^{#circ} ","l");
    leg->AddEntry(AntiMu2eDcp0, "#bar{#nu}, NH, #delta_{CP}=0","l");
    leg->AddEntry(AntiMu2eDcp270, "#bar{#nu}, NH, #delta_{CP}=270^{#circ}","l");
    leg->Draw();
    if(!isLBNE)canvas1->SaveAs(isNOVA?"plots/nova_threeflav_exact_probability_matter_numu2nue_pmnsinput.eps":"plots/t2k_threeflav_exact_probability_matter_numu2nue_pmnsinput.eps");
        else canvas1->SaveAs("plots/lbne_threeflav_exact_probability_matter_numu2nue_pmnsinput.eps");
            

    
    

    
        

}
