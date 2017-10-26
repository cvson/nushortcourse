//try to make a toy monte carlo
//pi calculation
void pi(){
    TCanvas *c1 = new TCanvas("c1","pi",800,800);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    //c1->SetGrid();
    gBenchmark->Start("hpi");
    gSystem->Unlink("hpianim.gif"); // delete old file
    

    
    //emty histogram
    TH2* hempty = new TH2F("","",100,0,1,100,0,1);
    //condition to print
    const Int_t kUPDATE = 50;
    //number of iteration
    const Int_t kIteration = 2000;
    //random number
    double xcord, ycord;
    //estimate pi and error
    double estimatePi;
    double errorPercentage;
    //counting and print
    Int_t nhit=0;
    Int_t nhit_in=0;
    Int_t nhit_out=0;
    TString resultString;
    TLatex *   tex;
    gRandom->SetSeed();
    hempty->Draw();
    //a quarter of circle
    TArc *parc = new TArc(0,0,1,0,90);
    parc->SetLineWidth(3);
    Int_t ci;

    //brow
    ci = TColor::GetColor("#663300");
    parc->SetLineColor(ci);
    parc->Draw("same");
    gPad->RedrawAxis();
    for ( Int_t i=0; i<kIteration; i++) {
        xcord =  gRandom->Uniform(0.0,1.0);
         ycord =  gRandom->Uniform(0.0,1.0);
        TMarker *mmarker = new TMarker(xcord, ycord, 20);
        mmarker->SetMarkerSize(0.5);
        if ((pow(xcord,2.0)+pow(ycord,2.0))<=1.0) {
            ++nhit_in;
            //UT orange 180,95,4
            ci = TColor::GetColor("#B45F04");
            mmarker->SetMarkerColor(ci);
        }
        else {
            ++nhit_out;
            //navy blue 11,56,97
            ci = TColor::GetColor("#0B3861");
           mmarker->SetMarkerColor(ci);
            
        }
        mmarker->DrawClone("same");
        //to save into histogram
        if (i && (i%kUPDATE) == 0) {

            c1->Update();
            estimatePi = nhit_in*4.0/(nhit_in+nhit_out);
            errorPercentage=100*TMath::Abs((estimatePi-TMath::Pi()))/TMath::Pi();
            cout<<"Estimate Pi "<<estimatePi<<" error percentage "<<errorPercentage<<endl;
            resultString.Form("#pi = %.4g, %.2g (%%) error after n=%d", estimatePi, errorPercentage,i);
            tex = new TLatex(0.25,0.973,resultString);
            tex->SetNDC();
            tex->SetTextAlign(13);
            tex->SetTextFont(42);
            tex->SetTextSize(0.04);
            tex->SetLineWidth(2);
            tex->Draw();
            
            if (gROOT->IsBatch()) {
                c1->Print("hpianim.gif+");
                printf("i = %d\n", i);
            }
            tex->Delete();
        }


    }
    
    cout<<"nin " <<nhit_in<<" nout "<<nhit_out<<endl;
    // make infinite animation by adding "++" to the file name
    if (gROOT->IsBatch()) c1->Print("hpianim.gif++");
    gBenchmark->Show("hpi");
}