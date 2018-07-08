const Int_t NEntry = 100000;
const double meanvalue = 25;
const double rmsvalue = 5;
void generate_data(){
    
    ofstream opfile("data4histogram.txt");
    
    opfile<<"#entry #value"<<endl;
    for (Int_t ientry=0; ientry<NEntry; ++ientry) {
        opfile<<ientry<<" "<<gRandom->Gaus(meanvalue,rmsvalue)<<endl;
    }
    opfile.close();
    plot_histogram();
}
void plot_histogram(){
    ifstream cfile("data4histogram.txt");
    string dummyline;
    getline(cfile,dummyline);
    Int_t ithentry;
    double ithvalue;
    double maxvalue = meanvalue+3*rmsvalue;
    double minvalue = meanvalue-3*rmsvalue;
    TH1F* hhisto = new TH1F( "Histogram", " " ,30,minvalue,maxvalue);
    for (Int_t ientry=0; ientry<NEntry; ++ientry) {
        cfile >> ithentry >> ithvalue;
        cout<<"entry "<<ithentry<<" values "<<ithvalue<<endl;
        hhisto->Fill(ithvalue);
    }
    new TCanvas;
    hhisto->GetXaxis()->SetTitle("Values of variable X");
    hhisto->GetYaxis()->SetTitle("Number of entries");
    hhisto->Draw();
    gPad->Print("basic_histogram_100k.eps");
    
    new TCanvas;
    hhisto->Draw();
    gStyle->SetOptFit(1111);
    hhisto->Fit("gaus");
    TPaveStats *st1 = (TPaveStats*)hhisto->GetListOfFunctions()->FindObject("stats");
    st1->SetX1NDC(0.65);
    st1->SetX2NDC(0.95);
    st1->SetY1NDC(0.65);
    st1->SetY2NDC(0.95);
    gPad->Update();
    gPad->Print("basic_histogram_100k_wfit.eps");
    

}
