const Int_t NEntry = 10000;
const double meanvalue = 25;
const double meanvalue2 = 20;
const double rmsvalue = 5;
const double rmsvalue2 = 4;

TMatrixD* produceSqrtMat( const TMatrixD& covMat )
{
    Double_t sum = 0;
    Int_t size = covMat.GetNrows();;
    TMatrixD* sqrtMat = new TMatrixD( size, size );
    
    for (Int_t i=0; i< size; i++) {
        
        sum = 0;
        for (Int_t j=0;j< i; j++) sum += (*sqrtMat)(i,j) * (*sqrtMat)(i,j);
        
        (*sqrtMat)(i,i) = TMath::Sqrt(TMath::Abs(covMat(i,i) - sum));
        
        for (Int_t k=i+1 ;k<size; k++) {
            
            sum = 0;
            for (Int_t l=0; l<i; l++) sum += (*sqrtMat)(k,l) * (*sqrtMat)(i,l);
            
            (*sqrtMat)(k,i) = (covMat(k,i) - sum) / (*sqrtMat)(i,i);
            
        }
    }
    return sqrtMat;
}

void getGaussRnd( TArrayD& v, const TMatrixD& sqrtMat, TRandom& R )
{
    // generate "size" correlated Gaussian random numbers
    
    // sanity check
    const Int_t size = sqrtMat.GetNrows();
    if (size != v.GetSize())
        cout << "<getGaussRnd> too short input vector: " << size << " " << v.GetSize() << endl;
    
    Double_t* tmpVec = new Double_t[size];
    
    for (Int_t i=0; i<size; i++) {
        Double_t x, y, z;
        y = R.Rndm();
        z = R.Rndm();
        x = 2*TMath::Pi()*z;
        tmpVec[i] = TMath::Sin(x) * TMath::Sqrt(-2.0*TMath::Log(y));
    }
    
    for (Int_t i=0; i<size; i++) {
        v[i] = 0;
        for (Int_t j=0; j<=i; j++) v[i] += sqrtMat(i,j) * tmpVec[j];
    }
    
    delete[] tmpVec;
}

void generate_data_2d(){
    gROOT->ProcessLine(".x ../rootlogon.C");
    ofstream opfile("data4histogram_2d.txt");
    
    opfile<<"#entry #value1 #value2"<<endl;
    //covariance matrix
    TMatrixD* covmatrix = new TMatrixD(2,2);
    /*Float_t offaxisVal = 0.45*rmsvalue*rmsvalue2;*/
    Float_t offaxisVal = 0.0*rmsvalue*rmsvalue2;
    (*covmatrix)(0,0) = rmsvalue*rmsvalue;
    (*covmatrix)(1,1) = rmsvalue2*rmsvalue2;
    (*covmatrix)(0,1) = offaxisVal*offaxisVal;
    (*covmatrix)(1,0) = (*covmatrix)(0,1);
    cout << "covariance matrix: " << endl;
    covmatrix->Print();
    
    //square-root matrix
    TMatrixD* sqrtMatrix = produceSqrtMat(  *covmatrix );
    cout << "Square root of covariance matrix: " << endl;
    sqrtMatrix->Print();

    TArrayD* v = new TArrayD( 2 );
    Float_t xvar[2];
    Float_t xS[2] = {  meanvalue,  meanvalue2 };
    
    
    
    //covariance matrix
    TRandom R( 100 );
    for (Int_t ientry=0; ientry<NEntry; ++ientry) {
        getGaussRnd( *v, *sqrtMatrix, R );
        for (Int_t ivar=0; ivar<2; ivar++) xvar[ivar] = (*v)[ivar] + xS[ivar];
        opfile<<ientry<<" "<<xvar[0]<<" "<<xvar[1]<<endl;
    }
    opfile.close();
    plot_histogram_2d();
}
void plot_histogram_2d(){
    ifstream cfile("data4histogram_2d.txt");
    string dummyline;
    getline(cfile,dummyline);
    Int_t ithentry;
    double ithvalue;
    double ithvaluey;
    double maxvalue = meanvalue+4*rmsvalue;
    double minvalue = meanvalue-4*rmsvalue;
    double maxvaluey = meanvalue2+4*rmsvalue2;
    double minvaluey = meanvalue2-4*rmsvalue2;
    TH2F* hhisto = new TH2F( "Histogram", " " ,8*rmsvalue,minvalue,maxvalue,8*rmsvalue2,minvaluey,maxvaluey);
    for (Int_t ientry=0; ientry<NEntry; ++ientry) {
        cfile >> ithentry >> ithvalue>>ithvaluey;
        //cout<<"entry "<<ithentry<<" values "<<ithvalue<<endl;
        hhisto->Fill(ithvalue,ithvaluey);
    }
    hhisto->GetXaxis()->SetTitle("Values of variable X");
    hhisto->GetYaxis()->SetTitle("Values of variable Y");
    TH1D * projh2X = hhisto->ProjectionX();
    TH1D * projh2Y = hhisto->ProjectionY();

    
    TCanvas *c1 = new TCanvas("c1", "c1",900,900);
    gStyle->SetOptStat(0);
    
    // Create the three pads
    TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
    center_pad->Draw();
    
    right_pad = new TPad("right_pad", "right_pad",0.6,0.0,1.0,0.6);
    right_pad->Draw();
    
    top_pad = new TPad("top_pad", "top_pad",0.0,0.6,0.55,1.0);
    top_pad->Draw();
    
    topright_pad = new TPad("topright_pad", "topright_pad",0.6,1.0,0.57,1.0);
    topright_pad->Draw();
    
    
    hhisto->SetMarkerStyle(8);
    hhisto->SetMarkerSize(0.3);
    
    // Drawing
    center_pad->cd();
    gStyle->SetOptStat(1);
    //gStyle->SetPalette(1);
    hhisto->Draw("SCAT");
    center_pad->Update();
    TPaveStats *st1 = (TPaveStats*)hhisto->GetListOfFunctions()->FindObject("stats");
    
    //topright_pad->cd();
    st1->SetX1NDC(0.2);
    st1->SetX2NDC(0.5);
    //st1->SetY1NDC(0.15);
    //st1->SetY2NDC(0.45);
    //st1->Draw();
    
    top_pad->cd();
    gStyle->SetOptStat(1);
    projh2X->SetFillColor(kBlue+1);
    projh2X->Draw("bar");
    top_pad->Update();
    TPaveStats *st2 = (TPaveStats*)projh2X->GetListOfFunctions()->FindObject("stats");
    st2->SetX1NDC(0.65);
    st2->SetX2NDC(0.95);
    
    right_pad->cd();
    gStyle->SetOptStat(1);
    projh2Y->SetFillColor(kRed+1);
    projh2Y->Draw("hbar");
    right_pad->Update();
    TPaveStats *st3 = (TPaveStats*)projh2Y->GetListOfFunctions()->FindObject("stats");
    st3->SetX1NDC(0.5);
    st3->SetX2NDC(0.85);
    
    //hhisto->Draw("colz");
    c1->Print("basic_histogram_2d_10k.eps");
    
    center_pad->cd();
    hhisto->Draw("COLZ");
    c1->Modified();
     c1->Print("basic_histogram_2d_10k_colz.eps");
    
    /*new TCanvas;
    hhisto->Draw();
    gStyle->SetOptFit(1111);
    hhisto->Fit("gaus");
    TPaveStats *st1 = (TPaveStats*)hhisto->GetListOfFunctions()->FindObject("stats");
    st1->SetX1NDC(0.65);
    st1->SetX2NDC(0.95);
    st1->SetY1NDC(0.65);
    st1->SetY2NDC(0.95);
    gPad->Update();
    gPad->Print("basic_histogram_100k_wfit.eps");*/
    

}
