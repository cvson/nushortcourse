//
//  covid19_trans.C
////////////////////////////////////////////////////
//
//  Simple model for covid19 transmission
//  Use with ROOT developed by CERN
//  run: root -b -q covid19_trans.C
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
#include <vector>
void covid19_trans(){
    TCanvas *c1 = new TCanvas("c1","covid19trans",1200,800);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    //divide into two part
    c1->Divide(2,1);
    const Bool_t IsDynamic = true;
    gBenchmark->Start("hcovid19trans");
    //to delete the old file
    gSystem->Unlink(IsDynamic?"hcovid19transanim_dynamic.gif":"hcovid19transanim_static.gif");
    
   
    const Int_t NoPopulation = 1000;
    const Double_t DisttobeInfect = 0.05;
    const Double_t ProbGetInfect2m = 0.5;
    const Double_t maxDistanceMove = 0.2;
    
    //condition to print
    const Int_t kUPDATE = 1;
    //number of iteration
    const Int_t kIteration = 30;//
    
    Double_t xpos_t0[NoPopulation];
    Double_t ypos_t0[NoPopulation];
    Bool_t isInfected[NoPopulation];
   
    //emty histogram
    TH2* hempty = new TH2F("hempty","",100,0,1,100,0,1);
    //population postion, not used
    TH2* hpop_t0 = new TH2F("hpop_t0","",100,0,1,100,0,1);
    
    //empty histogram
    TH1* hinfected = new TH1D("hinfected","",100,0,kIteration);
    
    //distance one-to-one
    TH1* hdist_t0 = new TH1D("hdist_t0","",500,0,2);
    
    //random number
    Double_t xcord, ycord;
    
    //F0 case
    Double_t f0case_x = 0.5;
    Double_t f0case_y = 0.5;
    TMarker *mmarker_f0 = new TMarker(f0case_x, f0case_y, 20);
    mmarker_f0->SetMarkerSize(0.8);
    Int_t ci;
    
    ci = TColor::GetColor("#B45F04");
    mmarker_f0->SetMarkerColor(ci);
    for ( Int_t i=0; i<NoPopulation-1; i++) {
        xcord =  gRandom->Uniform(0.0,1.0);
         ycord =  gRandom->Uniform(0.0,1.0);
        xpos_t0[i] = xcord;
        ypos_t0[i] = ycord;
        hpop_t0->Fill(xpos_t0[i],ypos_t0[i]);
        isInfected[i] = false;
        
    }
    
    vector<Int_t> Index_of_infected;
    
    //the 1000th person is F0 infected
    xpos_t0[NoPopulation-1] = f0case_x;
    ypos_t0[NoPopulation-1]  = f0case_y;
    hpop_t0->Fill(xpos_t0[NoPopulation-1],ypos_t0[NoPopulation-1]);
    isInfected[NoPopulation-1] =true;
    Index_of_infected.push_back(NoPopulation-1);
   
    //Set Zero infected for all iteration
    Double_t aIteration[kIteration];
    Double_t aNeff[kIteration];
    for ( Int_t iter=0; iter<kIteration; iter++) {
        aIteration[iter]=iter+0.5;
        aNeff[iter] = 0;
    }
    //for the initial iteration
    aNeff[0] = Index_of_infected.size();
    
    //graph of the accumulated No. of infected
    TGraph *pgrNeff = new TGraph(kIteration,aIteration,aNeff);
    
    //calculate the distance to test the density
    Double_t rx_ij;
    Double_t ry_ij;
    Double_t dist_ij;
    for (Int_t i=0; i<NoPopulation; i++) {
        for (Int_t j=i+1; j<NoPopulation; j++) {
             rx_ij = xpos_t0[i] - xpos_t0[j];
             ry_ij = ypos_t0[i] - ypos_t0[j];
            dist_ij = TMath::Sqrt(pow(ry_ij,2.0)+pow(rx_ij,2.0));
            //cout<<"distant "<<dist_ij<<endl;
            hdist_t0->Fill(dist_ij);
        
        }
    }
    
   //initial histogram
    c1->cd(1);
    hempty->GetXaxis()->SetTitle("Person location (x)");
    hempty->GetYaxis()->SetTitle("Person location (y)");
    titleStyle2D(hempty);
    hempty->Draw();
    
    c1->cd(2);
    hinfected->GetYaxis()->SetTitle("No. of infected people");
    hinfected->GetXaxis()->SetTitle("No. of iteration");
    titleStyle1D(hinfected);
    hinfected->Draw();
    hinfected->GetYaxis()->SetRangeUser(0,1000);
    
    for ( Int_t iter=0; iter<kIteration; iter++) {
        //for iter==0, put the infected case F0
        c1->Update();
        if (iter==0) {
            c1->cd(1);
            for ( Int_t iper=0; iper<NoPopulation-1; iper++) {
                TMarker *mmarker = new TMarker(xpos_t0[iper], ypos_t0[iper], 20);
                mmarker->SetMarkerSize(0.8);
                ci = TColor::GetColor("#0B3861");
                mmarker->SetMarkerColor(ci);
                mmarker->DrawClone("same");
                delete mmarker;
            }
            mmarker_f0->DrawClone("same");
            c1->cd(2);
            pgrNeff->Draw("PL");
            
        }
        else {
            //incase of dynamic, the position of person can be changed randomly
            c1->cd(1);
            if (IsDynamic) {
                
                //new position: can be confine?
                for ( Int_t iper=0; iper<NoPopulation; iper++) {
                    //need to delete the previous position
                    TMarker *mmarker = new TMarker(xpos_t0[iper], ypos_t0[iper], 20);
                    mmarker->SetMarkerSize(0.8);
                    mmarker->SetMarkerColor(10);//white color
                    mmarker->DrawClone("same");
                    delete mmarker;
                    
                    xpos_t0[iper] += gRandom->Uniform(-1.*maxDistanceMove,maxDistanceMove);
                    //can not move out of the boundary
                    if (xpos_t0[iper]<0) {
                        xpos_t0[iper]= 0.0;
                    }
                    if (xpos_t0[iper]>1.0) {
                        xpos_t0[iper]= 1.0;
                    }
                    ypos_t0[iper] += gRandom->Uniform(-1.*maxDistanceMove,maxDistanceMove);
                    //can not move out of the boundary
                    if (ypos_t0[iper]<0.0) {
                        ypos_t0[iper]= 0.0;
                    }
                    if (ypos_t0[iper]>1.0) {
                        ypos_t0[iper]= 1.0;
                    }
                    
                    //replot new postion
                    TMarker *mmarker = new TMarker(xpos_t0[iper], ypos_t0[iper], 20);
                    mmarker->SetMarkerSize(0.8);
                    if(isInfected[iper]==false){ci = TColor::GetColor("#0B3861");
                        mmarker->SetMarkerColor(ci);}
                    else {ci = TColor::GetColor("#B45F04");
                        mmarker->SetMarkerColor(ci);}
                    mmarker->DrawClone("same");
                    delete mmarker;
                    
                    
                }
            }
            //count infected
           
            for (Int_t iper=0; iper<NoPopulation; iper++) {
                if (isInfected[iper]==false) {
                    //loop over the all infected and find minimum distance
                    Double_t mindist_to_infected=2;
                    for (int iinf=0; iinf<Index_of_infected.size();iinf++){
                        rx_ij = xpos_t0[iper] - xpos_t0[Index_of_infected[iinf]];
                        ry_ij = ypos_t0[iper] - ypos_t0[Index_of_infected[iinf]];
                       dist_ij = TMath::Sqrt(pow(ry_ij,2.0)+pow(rx_ij,2.0));
                        if (mindist_to_infected>dist_ij) {
                            mindist_to_infected = dist_ij;
                        }
                    }
                    //the probability to get infected
                    //if distant is less than 0.2, the probability to get infected is 50%
                    if (mindist_to_infected<DisttobeInfect) {
                        double prob = gRandom->Uniform(0.0,1.0);
                        if (prob>ProbGetInfect2m) {
                            isInfected[iper] = true;
                            Index_of_infected.push_back(iper);
                            //turn red for the infected person
                            TMarker *mmarker = new TMarker(xpos_t0[iper], ypos_t0[iper], 20);
                            mmarker->SetMarkerSize(0.8);
                            ci = TColor::GetColor("#B45F04");
                            mmarker->SetMarkerColor(ci);
                            mmarker->DrawClone("same");
                            delete mmarker;
                        }
                    }
                    
                }//check ineffected
            }//loop all
            //update the total infected
            aNeff[iter] = Index_of_infected.size();
            pgrNeff->SetPoint(iter,iter+0.5, aNeff[iter]);
            c1->cd(2);
            pgrNeff->Draw("PL");
            
            
        }//end else
        
           
        c1->Update();
        if (gROOT->IsBatch()) {
            c1->Print(IsDynamic?"hcovid19transanim_dynamic.gif+":"hcovid19transanim_static.gif+");
            printf("i = %d\n", iter);
        }
        
    }
    
  
    // make infinite animation by adding "++" to the file name
    if (gROOT->IsBatch()) c1->Print(IsDynamic?"hcovid19transanim_dynamic.gif++":"hcovid19transanim_static.gif++");
    gBenchmark->Show("hcovid19trans");
}


void titleStyle1D(TH1* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.1);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.1);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.1);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.1);
    
    h1->GetYaxis()->SetTitleOffset(1.3);
    h1->GetXaxis()->SetTitleOffset(1.0);
}

void titleStyle2D(TH2* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.1);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.1);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.1);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.1);
    
    h1->GetYaxis()->SetTitleOffset(1.1);
    h1->GetXaxis()->SetTitleOffset(1.0);
}
