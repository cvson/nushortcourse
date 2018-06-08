double histGetLevel(TH1D* h, double p)
{
    double quantiles[1];
    double probsum[1];
    probsum[0] = p;
    
    h->GetQuantiles( 1, quantiles, probsum);
    
    return quantiles[0];
}

TH1D * GetSubHisto(TH1 *fHistogram ,double min, double max, const char * name){
    int imin = fHistogram->FindBin(min);
    int imax = fHistogram->FindBin(max);
    int nbins = fHistogram->GetNbinsX();
    double xmin = fHistogram->GetXaxis()->GetXmin();
    double xmax = fHistogram->GetXaxis()->GetXmax();
    
    double * xb = new double[nbins+3]; // nbins+1 original bin edges + 2 new bins
    int n=0; // counter
    
    int domin=1;
    int domax=1;
    
    for (int i=1;i<nbins+2;i++)
    {
        double x0 = fHistogram->GetBinLowEdge(i);
        
        if (min<x0 && domin)
        {
            xb[n++]=min;
            domin=0;
        }
        else if (min==x0)
            domin=0;
        
        if (max<x0 && domax)
        {
            xb[n++]=max;
            domax=0;
        }
        else if (max==x0)
            domax=0;
        
        xb[n++]=x0;
    }
    // now define the new histogram
    TH1D * h0 = new TH1D(name,"",n-1,xb);
    for(int i=1;i<n;i++)
    {
        double x0 = h0->GetBinCenter(i);
        if(x0<min || x0>max)
            continue;
        
        int bin=fHistogram->FindBin(x0);
        h0->SetBinContent(i, fHistogram->GetBinContent(bin));
    }
    
    return h0;
    
}
TH1D * IntegratedHistogram(TH2 *fHistogram){
    int nz = 100;
    double zmax = fHistogram -> GetMaximum();
    double dz   = zmax / double(nz);
    double nx = fHistogram -> GetNbinsX();
    double ny = fHistogram -> GetNbinsY();
    TH1D * htemp = new TH1D("htemp","",nz,0,zmax);
    htemp -> SetXTitle("z");
    htemp -> SetYTitle("Integrated probability");
    htemp -> SetStats(kFALSE);
    // loop over histogram
    for (int ix = 1; ix <= nx; ix++)
        for (int iy = 1; iy <= ny; iy++)
        {
            int binmin = int(fHistogram -> GetBinContent(ix, iy) / dz);
            for (int i = binmin; i <= nz; i++)
                htemp -> SetBinContent(i, htemp -> GetBinContent(i) + fHistogram -> GetBinContent(ix, iy));
        }
    return htemp;
}

/*TH1D * IntegratedHistogram(TH1 *fHistogram){
    int nz = 100;
    double zmax = fHistogram -> GetMaximum();
    double dz   = zmax / double(nz);
    double nx = fHistogram -> GetNbinsX();
    TH1D * htemp = new TH1D("htemp","",nz,0,zmax);
    htemp -> SetXTitle("z");
    htemp -> SetYTitle("Integrated probability");
    htemp -> SetStats(kFALSE);
    // loop over histogram
    for (int ix = 1; ix <= nx; ix++)
        for (int iy = 1; iy <= ny; iy++)
        {
            int binmin = int(fHistogram -> GetBinContent(ix, iy) / dz);
            for (int i = binmin; i <= nz; i++)
                htemp -> SetBinContent(i, htemp -> GetBinContent(i) + fHistogram -> GetBinContent(ix, iy));
        }
    return htemp;
}*/

void smooth_th2(TH2 *h){
	// Sample kernel
	const Int_t ksize_x=5;
	const Int_t ksize_y=5;
	Double_t kernel[ksize_x][ksize_y] = { { 0, 0, 1, 0, 0 },
		{ 0, 2, 2, 2, 0 },
		{ 1, 2, 5, 2, 1 },
		{ 0, 2, 2, 2, 0 },
		{ 0, 0, 1, 0, 0 } };
	
	// Determine the size of the bin buffer(s) needed
	Int_t lowest_bin=h->GetBin(0,0);
	Int_t highest_bin=h->GetBin(h->GetNbinsX()+1,h->GetNbinsY()+1);
	Int_t bufSize = highest_bin-lowest_bin+1;
	Double_t *buf = new Double_t[bufSize];
	Double_t *ebuf = new Double_t[bufSize];
	
	// Copy all the data to the temporry buffers
	for (Int_t i=0; i<=h->GetNbinsX(); i++){
		for (Int_t j=0; j<=h->GetNbinsY(); j++){
			Int_t bin = h->GetBin(i,j);
			buf[bin] =h->GetBinContent(bin);
			ebuf[bin]=h->GetBinError(bin);
		};
	};
	
	// Kernel tail sizes (kernel sizes must be odd for this to work!)
	Int_t x_push = (ksize_x-1)/2;
	Int_t y_push = (ksize_y-1)/2;
	
	// main work loop
	for (Int_t i=1; i<=h->GetNbinsX(); i++){
		for (Int_t j=1; j<=h->GetNbinsY(); j++){
			Double_t content = 0.0;
			Double_t error = 0.0;
			Double_t norm = 0.0;
			
			for (Int_t n=0; n<ksize_x; n++){
				for (Int_t m=0; m<ksize_y; m++){
					Int_t xb = i+(n-x_push);
					Int_t yb = j+(m-y_push);
					if ( (xb >= 1) && (xb <= h->GetNbinsX()) &&
						(yb >= 1) && (yb <= h->GetNbinsY()) ){
						Int_t bin = h->GetBin(xb,yb);
						Double_t k = kernel[n][m];
						if ( (k != 0.0 ) && (buf[bin] != 0.0) ) { // General version probably does not want the second condition
							content += k*buf[bin];
							error   += k*k*buf[bin]*buf[bin];
							norm    += k;
						};
					};
				};
			};
			
			content /= norm;
			error /= (norm*norm);
			if ( content != 0.0 ) { // General version probably does not want this condition
				h->SetBinContent(i,j,content);
				h->SetBinError(i,j,sqrt(error));
			};
		};
	};
	
	delete buf;
	delete ebuf;
};


void pmnsEle_glob(){
    bool isNormalHierarchy=false;
    double sinsq12=0.307;
    
    double sinsq23N=0.386;
    double sinsq23I=0.392;
    double sinsq23 = isNormalHierarchy?sinsq23N:sinsq23I;
    
    double sinsq13N=0.0241;
    double sinsq13I=0.0244;
    double sinsq13 = isNormalHierarchy?sinsq13N:sinsq13I;
    
    double deltaN=1.08*TMath::Pi();
    double deltaI=1.09*TMath::Pi();
    double delta = isNormalHierarchy?deltaN:deltaI;
    
    cout<<"The global analysis input" <<endl;
    cout<<"sinsq12: "<<sinsq12<<endl;
    cout<<"sinsq23: Normal "<<sinsq23N<<" Invert "<<sinsq23I<<endl;
    cout<<"sinsq13: Normal "<<sinsq13N<<" Invert "<<sinsq13I<<endl;
    cout<<"deltacp: Normal "<<deltaN<<" Invert "<<deltaN<<endl;
    
    //write in sin and cos
    double s12=TMath::Sqrt(sinsq12);
    double c12=TMath::Sqrt(1-sinsq12);
    
    double s23=TMath::Sqrt(sinsq23);
    double c23=TMath::Sqrt(1-sinsq23);
    
    double s13=TMath::Sqrt(sinsq13);
    double c13=TMath::Sqrt(1-sinsq13);
    //using TComplex for CP phase
    double scp=TMath::Sin(delta);
    double ccp=TMath::Cos(delta);
    TComplex cpPhasePos= TComplex(ccp,scp);
    TComplex cpPhaseNeg= TComplex(ccp,-scp);
    
    
    if(isNormalHierarchy)cout<<"Normal hierarchy" <<endl;
    else cout<<"Invert hierarchy" <<endl;
    cout<<"s12: "<<s12<<" c12 "<<c12<<endl;
    cout<<"s23: "<<s23<<" c23 "<<c23<<endl;
    cout<<"s13: "<<s13<<" c13 "<<c13<<endl;
    cout<<"scp: "<<cpPhasePos.Re()<<" ccp "<<cpPhasePos.Im()<<endl;
    
    //PMNS element
    //electron
    double Ue1=c12*c13;
    cout<<"Ue1 "<<Ue1<<endl;
    double Ue2=s12*c13;
    cout<<"Ue2 "<<Ue2<<endl;
    TComplex Ue3=cpPhaseNeg*s13;
    
    //muon
    cout<<"Ue3 Re "<<Ue3.Re()<<" Im "<<Ue3.Im()<<endl;
    TComplex Umu1=-s12*c23-c12*s23*s13*cpPhasePos;
    cout<<"Umu1 Re "<<Umu1.Re()<<" Im "<<Umu1.Im()<<endl;
    TComplex Umu2=c12*c23-s12*s23*s13*cpPhasePos;
    cout<<"Umu2 Re "<<Umu2.Re()<<" Im "<<Umu2.Im()<<endl;
    double Umu3=s23*c13;
    cout<<"Umu3  "<<Umu3<<endl;
    
    //tau
    TComplex Utau1=s12*s23-c12*c23*s13*cpPhasePos;
    cout<<"Utau1 Re "<<Utau1.Re()<<" Im "<<Utau1.Im()<<endl;
    TComplex Utau2=-c12*s23-s12*c23*s13*cpPhasePos;
    cout<<"Utau2 Re "<<Utau2.Re()<<" Im "<<Utau2.Im()<<endl;
    double Utau3=c23*c13;
    cout<<"Utau3  "<<Utau3<<endl;

    //test unitary condition
    double norm_e=Ue1**2+Ue2**2+Ue3->Rho2();
    double norm_mu=Umu1->Rho2()+Umu2->Rho2()+Umu3**2;
    double norm_tau=Utau1->Rho2()+Utau2->Rho2()+Utau3**2;
    cout<<"norm e "<<norm_e<<endl;
    cout<<"norm mu "<<norm_mu<<endl;
    cout<<"norm tau "<<norm_tau<<endl;
    
    //calculate the triangle
    TComplex Ue1Umu1=Ue1*Umu1;//Ue1 is real
    cout<<"Ue1Umu1 Re "<<Ue1Umu1.Re()<<" Im "<<Ue1Umu1.Im()<<" Amp "<<Ue1Umu1.Rho()<<endl;
    TComplex Ue2Umu2=Ue2*Umu2;//Ue2 is real
    cout<<"Ue2Umu2 Re "<<Ue2Umu2.Re()<<" Im "<<Ue2Umu2.Im()<<" Amp "<<Ue2Umu2.Rho()<<endl;
    TComplex Ue3Umu3=(TComplex::Conjugate(Ue3))*Umu3;
    cout<<"Ue3Umu3 Re "<<Ue3Umu3.Re()<<" Im "<<Ue3Umu3.Im()<<" Amp "<<Ue3Umu3.Rho()<<endl;
    
    TComplex sumCheck = Ue1Umu1+Ue2Umu2+Ue3Umu3;
    cout<<"sumCheck Re "<<sumCheck.Re()<<" Im "<<sumCheck.Im()<<endl;
    
    // sides
    TComplex sideA=Ue1Umu1/Ue2Umu2;
    TComplex sideB=Ue3Umu3/Ue2Umu2;
    
    cout<<"sideA Re "<<sideA.Re()<<" Im "<<sideA.Im()<<" Amp "<<sideA.Rho()<<" angle "<<sideA.Theta()<<endl;
    cout<<"sideB Re "<<sideB.Re()<<" Im "<<sideB.Im()<<" Amp "<<sideB.Rho()<<" angle "<<sideB.Theta()<<endl;
    
    //Draw a triangle;
    c1 = new TCanvas("c1");
    gStyle->SetOptStat(0);
    c1->Range(-1.0,-1.0,1.0,1.0);
    TH2D *h2 = new TH2D("h2","",200,-0.5,1.2,200, -0.1,0.1);
    h2->GetXaxis()->SetTitle("#rho (Real)");
    h2->GetYaxis()->SetTitle("#eta (Imagine)");
    h2->GetYaxis()->SetTitleOffset(1.2);
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->CenterTitle();
    gPad->SetGrid();
    TLine *linec= new TLine(0,0,1,0);
    TLine *linea= new TLine(1,0,1+sideA.Re(),sideA.Im());
    TLine *lineb= new TLine(0,0,-sideB.Re(),-sideB.Im());
    linec->SetLineWidth(3);
    linea->SetLineWidth(3);
    lineb->SetLineWidth(3);
    h2->Draw("AXIS");
    linec->Draw("same");
    linea->Draw("same");
    lineb->Draw("same");
    TLatex vertexA(-0.1,-0.03,"A(0,0)");
    TLatex vertexB(0.9,-0.03,"B(0,1)");
    TLatex texsideA(-0.2,0.03,"#left|#frac{U_{e1}^{*}U_{#mu1}}{U_{e2}^{*}U_{#mu2}}#right|");
    TLatex texsideB(0.8,0.035,"#left|#frac{U_{e3}^{*}U_{#mu3}}{U_{e2}^{*}U_{#mu2}}#right|");
    
    TLatex vertexC(0.22,0.07,"C(0.241,0.057)");
    TLatex angleA(0.1,0.01,"#gamma =13^{o}");
    TLatex angleB(0.6,0.005,"#beta =4^{o}");
    angleA.SetTextColor(kAzure);
    angleB.SetTextColor(kAzure);
    vertexA.SetTextColor(95);
    vertexB.SetTextColor(95);
    vertexC.SetTextColor(95);
    vertexA.Draw("same");
    vertexB.Draw("same");
    texsideA.Draw("same");
    texsideB.Draw("same");
    if(isNormalHierarchy){
        vertexC.Draw("same");
        angleA.Draw("same");
        angleB.Draw("same");
    }
    else {
        vertexC.DrawLatex(0.22,0.07,"C(0.243,0.065)");
        angleA.DrawLatex(0.1,0.01,"#gamma =15^{o}");
        angleB.DrawLatex(0.6,0.005,"#beta =5^{o}");
    }
    if(isNormalHierarchy)c1->Print("triangle_normal.eps");
        else c1->Print("triangle_invert.eps");
            
    //------Get the uncertainty from simulation
    TString filename = isNormalHierarchy?"HistToy_normal.root":"HistToy_invert.root";
    TFile *ffile = new TFile(filename);
    TH2D *hCorellation = (TH2D*)ffile->Get("hCorrelationA");
    c2 = new TCanvas;
    hCorellation->Draw("colz");
    
    //-----Integrated histogram
    TH1D * hprob = IntegratedHistogram(hCorellation);
    c3=new TCanvas;
    hprob->Draw();
    double levels[2];
    double level5 = histGetLevel(hprob,0.05);
    levels[0] = 0.;
    levels[1] = level5;
    hCorellation->SetContour(2, levels);
    
    //-----Correlation plot
    c4=new TCanvas;
    c4->SetGrid();
    hCorellation->SetTitle("");
    hCorellation->GetXaxis()->SetTitle("#rho (Real)");
    hCorellation->GetYaxis()->SetTitle("#eta (Imagine)");
    hCorellation->GetXaxis()->CenterTitle();
    hCorellation->GetYaxis()->CenterTitle();
    hCorellation->GetXaxis()->SetRangeUser(-0.4,1.0);
    hCorellation->GetYaxis()->SetRangeUser(-0.4,0.4);
    hCorellation->Draw("cont3");
    linec->Draw("same");
    linea->Draw("same");
    lineb->Draw("same");
    TLegend *leg = new TLegend(0.5,0.65,0.85,0.85);
    leg->SetFillStyle(0);
	leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    leg->AddEntry(hCorellation,"95% credible interval","l");
    leg->AddEntry(linea,"Triangle at best fit","l");
    leg->Draw("same");
    if(isNormalHierarchy)c4->Print("triangle_wUncertainty_normal.eps");
    else c4->Print("triangle_wUncertainty_invert.eps");
    
    
    //-----rho distribution
    c5 = new TCanvas;
    TH1D *hrho = (TH1D*)ffile->Get("hrhoA");
    Double_t quantileUp;
    Double_t quantileDw;
    Double_t media;
    Double_t prodMedia=0.50;
    Double_t probUp=.84;
    Double_t probDw=.16;
    hrho->GetQuantiles(1,&quantileUp,&probUp);
    hrho->GetQuantiles(1,&quantileDw,&probDw);
    hrho->GetQuantiles(1,&media,&prodMedia);
    cout<<"median "<<media<<" quantile up "<<quantileUp<<" down "<<quantileDw<<endl;
    TH1 *hrhoInt =GetSubHisto(hrho,quantileDw,quantileUp,"hrhoInt");
    hrhoInt->SetFillStyle(1001);
    hrhoInt->SetFillColor(kYellow);
    hrho->SetTitle("");
    hrho->GetXaxis()->SetTitle("#rho (real) distribution");
    hrho->GetYaxis()->SetTitleOffset(1.2);
    hrho->GetYaxis()->SetTitle("Arbitary scale");
    hrho->GetXaxis()->CenterTitle();
    hrho->GetYaxis()->CenterTitle();
    hrho->Draw();
    hrhoInt->Draw("same");
    //Tlatex
    TLatex * tmax_text = new TLatex();
    tmax_text->SetTextSize(0.05);
    tmax_text->SetTextFont(62);
    tmax_text->SetTextAlign(22); // center
    /*tmax_text->DrawLatex(-0.6,6000,
                         TString::Format( TString::Format("#eta ^{med} = %%.%dg ^{+%%.2g}_{ -%%.2g}",sd),
                                         media, quantileUp-media, media-quantileDw));*/
    tmax_text->DrawLatex(-0.5,6000,TString::Format("#rho ^{med} = %.2f^{+%.2f}_{-%.2f} ",media,quantileUp-media,media-quantileDw));
    TLegend* leg1 = new TLegend(0.16, 0.65, 0.4, 0.85);
    leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
    leg1->SetTextSize(0.05);
    TLegend* band = new TLegend(0, 0, 1, 1);
    band->SetFillColor(kYellow);
    leg1->AddEntry(band,"68% credible interval","F");
    leg1->Draw("same");
    if(isNormalHierarchy)c5->Print("rho_wUncertainty_normal.eps");
    else c5->Print("rho_wUncertainty_invert.eps");
    
    
    //-----eta distribution
    c6 = new TCanvas;
    TH1D *heta = (TH1D*)ffile->Get("hetaA");
    heta->SetTitle("");
    heta->GetXaxis()->SetTitle("#eta (imagine) distribution");
    heta->GetYaxis()->SetTitleOffset(1.2);
    heta->GetYaxis()->SetTitle("Arbitary scale");
    heta->GetXaxis()->CenterTitle();
    heta->GetYaxis()->CenterTitle();
    heta->Draw();
    Double_t etaquantileUp;
    Double_t etaquantileDw;
    Double_t etamedia;
    //68%
    heta->GetQuantiles(1,&etaquantileUp,&probUp);
    heta->GetQuantiles(1,&etaquantileDw,&probDw);
    heta->GetQuantiles(1,&etamedia,&prodMedia);
    cout<<"median "<<etamedia<<" etaquantile up "<<etaquantileUp<<" down "<<etaquantileDw<<endl;
    TH1 *hetaInt =GetSubHisto(heta,etaquantileDw,etaquantileUp,"hetaInt");
    hetaInt->SetFillStyle(1001);
    hetaInt->SetFillColor(kYellow);
    hetaInt->Draw("same");
    tmax_text->DrawLatex(-0.55,1800,TString::Format("#eta ^{med} = %.2f^{+%.2f}_{-%.2f} ",etamedia,etaquantileUp-etamedia,etamedia-etaquantileDw));
    leg1->Draw("same");
    if(isNormalHierarchy)c6->Print("eta_wUncertainty_normal.eps");
    else c6->Print("eta_wUncertainty_invert.eps");

    
}