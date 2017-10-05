{

#include<TH1>
#include<TH2>
#include<TLatex>
#include<TString>
    
    gROOT->ProcessLine(".L OscCalculatorPMNS.C+");
  
    double x=0.0; //Energy
    Bool_t isNOVA=true;
    Bool_t isLBNE=false;//T2K if false
  
    /* namespace OscPar
     {
     typedef enum EOscPar{
     kTh12 = 0,
     kTh23 = 1,
     kTh13 = 2,
     kDeltaM23 = 3,
     kDeltaM12 = 4,
     kDelta = 5,
     kDensity = 6,
     kL = 7,
     kNuAntiNu = 8,
     kUnknown = 9
     } OscPar_t;
    */
  
  
    OscCalculatorPMNS myProb;

    double myDeltam2=2.42;
    double myDeltam21=7.54;
    double myTheta12=0.61;
    double myTheta23=0.78;
    double myTheta13=0.15;
    double myL=isNOVA? 810:295.;
    if(isLBNE) myL = 1250;
    double ent2k=isNOVA? 2.0:0.6;
    if(isLBNE) ent2k =2.0;
    
    double xshift = 0.06;
    double yshift = 0.02;
    TString delMStr32=Form("%2.2f #times 10^{-3}",myDeltam2);
    TString legStr32="|#Delta m^{2}_{32}|=";
    TLatex* tlx=new TLatex(0.6+xshift, 0.85+yshift,legStr32+delMStr32);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.04);
    tlx->SetTextAlign(12);
	
	TString delMStr21=Form("%2.2f #times 10^{-5}",myDeltam21);
    TString legStr21="|#Delta m^{2}_{21}|=";
    TLatex* tlx1=new TLatex(0.6+xshift, 0.80+yshift,legStr21+delMStr21);
    tlx1->SetNDC(kTRUE); // <- use NDC coordinate
    tlx1->SetTextSize(0.04);
    tlx1->SetTextAlign(12);
	
	
    TString Th23Str=Form("sin^{2}2#theta_{23}=%2.2f",TMath::Sin(2*myTheta23)*TMath::Sin(2*myTheta23));
    TLatex* tlx2=new TLatex(0.6+xshift, 0.75+yshift,Th23Str);
    tlx2->SetNDC(kTRUE); // <- use NDC coordinate
    tlx2->SetTextSize(0.04);
    tlx2->SetTextAlign(12);
	
	TString Th12Str=Form("sin^{2}2#theta_{12}=%2.2f",TMath::Sin(2*myTheta12)*TMath::Sin(2*myTheta12));
    TLatex* tlx3=new TLatex(0.6+xshift, 0.70+yshift,Th12Str);
    tlx3->SetNDC(kTRUE); // <- use NDC coordinate
    tlx3->SetTextSize(0.04);
    tlx3->SetTextAlign(12);
	
	TString Th13Str=Form("sin^{2}2#theta_{13}=%2.2f",TMath::Sin(2*myTheta13)*TMath::Sin(2*myTheta13));
    TLatex* tlx4=new TLatex(0.6+xshift, 0.65+yshift,Th13Str);
    tlx4->SetNDC(kTRUE); // <- use NDC coordinate
    tlx4->SetTextSize(0.04);
    tlx4->SetTextAlign(12);
	
    TString LStr=Form("L=%2.0f km",myL);
    TLatex* tlx5=new TLatex(0.6+xshift, 0.60+yshift,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    tlx5->SetTextSize(0.04);
    tlx5->SetTextAlign(12);
    
    TString EnStr=Form("E_{#nu}=%2.1f GeV",ent2k);
    TLatex* tlx6=new TLatex(0.6+xshift, 0.55+yshift,EnStr);
    tlx6->SetNDC(kTRUE); // <- use NDC coordinate
    tlx6->SetTextSize(0.04);
    tlx6->SetTextAlign(12);
    
    TString NorHierStr="#color[864]{#Delta m^{2}_{32}>0}";
    TLatex* tlx7=new TLatex(isNOVA? 0.67:0.72, 0.45,NorHierStr);
    tlx7->SetNDC(kTRUE); // <- use NDC coordinate
    tlx7->SetTextSize(0.04);
    tlx7->SetTextAlign(12);
    
    TString InvHierStr="#color[924]{#Delta m^{2}_{32}<0}";
    TLatex* tlx8=new TLatex(isNOVA? 0.25:0.3, 0.6,InvHierStr);
        if(isLBNE)tlx8=new TLatex(0.2, 0.6,InvHierStr);
    tlx8->SetNDC(kTRUE); // <- use NDC coordinate
    tlx8->SetTextSize(0.04);
    tlx8->SetTextAlign(12);
    
    //Set all parameters here
    myProb.SetL(myL);
  
    myProb.SetRho(2.71);
    myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
    myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));

    myProb.SetTh12(myTheta12);
    myProb.SetTh23(myTheta23);
    myProb.SetTh13(myTheta13);
    myProb.SetdCP(0.0001);

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
	    << myProb.GetDmsq21() << endl  
	    << myProb.GetDmsq32() << endl  
	    << myProb.GetTh12() << endl  
	    << myProb.GetTh23() << endl  
	    << myProb.GetTh13() << endl  
	    << myProb.GetdCP() << endl  
	    << endl;

    //This step is done automatically when computing probabilities,
    //but useful to do it here to check code works
    myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
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
    Int_t dcpstep=200;
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    // normal hierarchy
    TGraph *MuElec=0;
	Double_t pmuex[dcpstep+1], pmuey[dcpstep+1];
    //neutrino
    isAnti=1;
    for (Int_t i=0; i < dcpstep+1; i++){
        double vdcp = i*2.0*TMath::Pi()/dcpstep + epsilon;
        myProb.Reset();
        myProb.SetL(myL);
        
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(vdcp);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        
        pmuex[i]=myProb.P(isAnti*14,isAnti*12,ent2k);
        cout<< "cp "<<myProb.GetdCP() <<" prob "<< pmuex[i] <<endl;

    }
    //Anti neutrino
    isAnti=-1;
    for (Int_t i=0; i < dcpstep+1; i++){
        double vdcp = i*2.0*TMath::Pi()/dcpstep + epsilon;
        myProb.Reset();
        myProb.SetL(myL);
        
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
        
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(vdcp);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        pmuey[i]=myProb.P(isAnti*14,isAnti*12,ent2k);
        
    }
    
    MuElec = new TGraph(dcpstep,pmuex,pmuey);
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    // inverted hierarchy
    
	TGraph *MuElecInv=0;
    Double_t pmuexInv[dcpstep+1], pmueyInv[dcpstep+1];

    //neutrino
    isAnti=1;
    for (Int_t i=0; i < dcpstep+1; i++){
        double vdcp = i*2.0*TMath::Pi()/dcpstep + epsilon;
        myProb.Reset();
        myProb.SetL(myL);
        
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(-1.*myDeltam2*TMath::Power(10,-3.));
        
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(vdcp);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        pmuexInv[i]=myProb.P(isAnti*14,isAnti*12,ent2k);
    }
    //Anti neutrino
    isAnti=-1;
    for (Int_t i=0; i < dcpstep+1; i++){
        double vdcp = i*2.0*TMath::Pi()/dcpstep + epsilon;
        myProb.Reset();
        myProb.SetL(myL);
        
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(-1.*myDeltam2*TMath::Power(10,-3.));
        
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(vdcp);
        
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
        pmueyInv[i]=myProb.P(isAnti*14,isAnti*12,ent2k);
    }
    
    MuElecInv = new TGraph(dcpstep,pmuexInv,pmueyInv);
    
    //UT orange
    Int_t ci1 =  TColor::GetColor("#B45F04");
    //navy blue
    Int_t ci2 = TColor::GetColor("#006699");
    MuElec->SetLineColor(ci2);
    MuElec->SetLineWidth(2);
    MuElecInv->SetLineColor(ci1);
    MuElecInv->SetLineWidth(2);
    
    //UT orange
    Int_t ci4 =  TColor::GetColor("#CC3300");
    //navy blue
    Int_t ci3 = TColor::GetColor("#339900");
    //cp=0
    TMarker *MuEleccp0= new TMarker(pmuex[0],pmuey[0],21);
    TMarker *MuElecInvcp0= new TMarker(pmuexInv[0],pmueyInv[0],21);
    MuEleccp0->SetMarkerSize(1.3);
    MuElecInvcp0->SetMarkerSize(1.3);
    MuEleccp0->SetMarkerColor(ci3);
    MuElecInvcp0->SetMarkerColor(ci3);
    //cp=pi/2
    TMarker *MuEleccp12= new TMarker(pmuex[dcpstep/4],pmuey[dcpstep/4],8);
    TMarker *MuElecInvcp12= new TMarker(pmuexInv[dcpstep/4],pmueyInv[dcpstep/4],8);
    MuEleccp12->SetMarkerSize(1.3);
    MuElecInvcp12->SetMarkerSize(1.3);
    MuEleccp12->SetMarkerColor(ci3);
    MuElecInvcp12->SetMarkerColor(ci3);
    
    
    //cp=pi
    TMarker *MuEleccp1= new TMarker(pmuex[dcpstep/2],pmuey[dcpstep/2],25);
    TMarker *MuElecInvcp1= new TMarker(pmuexInv[dcpstep/2],pmueyInv[dcpstep/2],25);
    MuEleccp1->SetMarkerSize(1.3);
    MuElecInvcp1->SetMarkerSize(1.3);
    MuEleccp1->SetMarkerColor(ci4);
    MuElecInvcp1->SetMarkerColor(ci4);
    
    //cp=3pi/2
    TMarker *MuEleccp32= new TMarker(pmuex[dcpstep*3/4],pmuey[dcpstep*3/4],24);
    TMarker *MuElecInvcp32= new TMarker(pmuexInv[dcpstep*3/4],pmueyInv[dcpstep*3/4],24);
    MuEleccp32->SetMarkerSize(1.3);
    MuElecInvcp32->SetMarkerSize(1.3);
    MuEleccp32->SetMarkerColor(ci4);
    MuElecInvcp32->SetMarkerColor(ci4);
    
    TCanvas *canvas1=0;
    canvas1 = new TCanvas("Can1","Can1",800,600);
    
    TH2F *histo1=0;
    histo1 = new TH2F("NewOscProb","NewOscProb",100, 0., 0.09, 100, 0.0, 0.09);
	
	histo1->SetTitle("");
    //histo1->GetYaxis()->SetRangeUser(0.0,0.09);
    histo1->GetXaxis()->SetTitle("P(#nu_{#mu}#rightarrow #nu_{e})");
    histo1->GetYaxis()->SetTitle("P(#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e})");
    histo1->GetXaxis()->CenterTitle();
    histo1->GetYaxis()->CenterTitle();
    histo1->GetXaxis()->SetLabelSize(histo1->GetXaxis()->GetTitleSize()*1.2);
    histo1->GetYaxis()->SetLabelSize(histo1->GetYaxis()->GetTitleSize()*1.2);
    histo1->GetXaxis()->SetTitleSize(histo1->GetXaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleSize(histo1->GetYaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    //TGaxis::SetMaxDigits(3);
    canvas1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.03,"xy");
    gStyle->SetTitleSize(0.04,"xy");
    //gStyle->SetTitleOffset(1.2,"xy");
    
    TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
    leg->SetFillStyle(0);
	leg->SetBorderSize(0);
    leg->AddEntry(MuEleccp0,"#delta_{CP}=0","p");
    leg->AddEntry(MuEleccp12,"#delta_{CP}=#pi/2","p");
    leg->AddEntry(MuEleccp1,"#delta_{CP}=#pi","p");
    leg->AddEntry(MuEleccp32,"#delta_{CP}=3#pi/2","p");
    
    histo1->Draw();
    MuElec->Draw("same");
    MuElecInv->Draw("same");
    MuEleccp0->Draw("same");
    MuElecInvcp0->Draw("same");
    MuEleccp12->Draw("same");
    MuElecInvcp12->Draw("same");
    MuEleccp1->Draw("same");
    MuElecInvcp1->Draw("same");
    MuEleccp32->Draw("same");
    MuElecInvcp32->Draw("same");
    tlx->Draw("same");
    tlx1->Draw("same");
    tlx2->Draw("same");
    tlx3->Draw("same");
    tlx4->Draw("same");
    tlx5->Draw("same");
    tlx6->Draw("same");
    tlx7->Draw("same");
    tlx8->Draw("same");
    leg->Draw("same");
    if(!isLBNE)canvas1->SaveAs(isNOVA?"plots/nova_elec_appearance_biprobability.eps":"plots/t2k_elec_appearance_biprobability.eps");
        else canvas1->SaveAs("plots/lbne_elec_appearance_biprobability.eps");

    
        

}
