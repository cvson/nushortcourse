void titleStyleGraph(TGraph* h1){
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.2);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.2);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetXaxis()->SetTitleOffset(0.9);
}

void vn_covid_plot(){
    //https://www.worldometers.info/coronavirus/country/viet-nam/#graph-cases-daily
    //use chrome to get source file then get data
    TString savename="vn_covid_20210819";
    Double_t newCase[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,10,1,3,5,5,0,9,4,4,5,10,9,6,3,19,10,11,14,5,10,11,20,10,8,6,15,6,1,1,4,6,0,4,2,1,4,3,1,2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,1,0,0,0,17,0,0,0,0,0,0,24,2,4,2,4,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,1,0,3,0,0,0,1,1,0,0,0,1,7,7,0,0,0,0,3,0,1,2,0,0,0,0,0,0,0,0,14,0,0,0,0,1,2,0,1,8,0,1,0,1,1,17,7,4,1,4,3,11,15,13,50,37,44,30,32,20,45,30,42,23,29,6,19,17,28,19,21,11,21,6,5,13,2,5,2,6,7,5,2,2,2,0,4,0,2,0,3,0,0,0,5,5,0,1,0,3,0,0,0,3,2,0,0,0,0,1,0,0,0,5,3,17,0,1,1,0,0,1,1,1,1,5,2,2,1,3,9,2,0,2,8,6,1,3,4,0,12,8,1,3,1,4,0,3,0,12,10,1,4,5,1,0,2,11,26,1,3,0,25,2,5,12,4,1,1,1,5,4,5,10,8,2,2,4,4,7,3,0,4,1,1,10,4,4,6,4,2,5,3,0,2,3,1,2,1,6,1,11,7,1,1,10,3,2,9,9,8,12,3,7,1,4,3,1,1,1,5,1,10,5,1,0,2,1,4,2,2,0,0,1,2,0,91,63,62,50,33,32,29,46,19,5,20,49,14,27,49,2,53,33,41,42,18,18,15,6,15,9,9,11,8,6,6,16,0,24,10,6,6,7,11,12,2,3,4,17,3,1,3,3,7,3,1,1,0,3,0,1,3,7,0,5,3,0,9,14,3,6,5,6,11,11,9,15,9,1,12,9,19,25,14,9,3,7,9,6,10,14,3,10,9,5,8,45,18,14,20,19,15,26,68,47,93,102,129,76,86,87,106,169,190,184,153,178,119,132,145,131,187,369,313,230,254,286,251,214,251,241,250,224,254,206,236,175,407,219,196,261,297,272,402,423,515,264,486,311,272,244,220,285,868,175,368,398,372,450,713,545,922,890,1102,1029,1007,1314,1625,1853,1953,2383,2301,2934,3416,3336,3718,5926,4195,4795,5357,6194,7307,9256,7531,7882,7913,6559,7594,8649,8624,8620,7455,8429,7623,7244,8324,7334,9690,9340,8390,8766,9667,9180,9716,9580,8652,9605};
    
    Double_t newDeadthCase[]={
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,3,0,2,0,2,0,0,1,4,1,1,4,1,1,1,0,1,0,0,0,1,1,0,0,2,1,0,2,0,2,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,2,2,0,1,2,0,1,1,1,0,0,0,1,1,0,2,2,0,0,2,0,0,2,1,1,0,2,0,0,1,2,2,3,0,1,2,2,0,2,0,4,1,0,3,0,2,4,7,5,3,5,2,7,6,7,6,69,18,0,29,80,0,36,26,32,34,31,31,0,106,233,298,145,192,197,376,256,393,296,234,147,360,388,342,326,275,349,337,367,331 };
    
Int_t NDATE_D =sizeof(newDeadthCase)/sizeof(newDeadthCase[0]);
    
    Int_t NDATE =sizeof(newCase)/sizeof(newCase[0]);

cout<<"Data size fr newCase "<<NDATE<<" fr Deadcase"<<NDATE_D<<endl;
    
    Double_t* dateIndex = new Double_t[NDATE];
    Double_t* average7day = new Double_t[NDATE];
    for (Int_t idate =0; idate<NDATE; ++idate) {
        //dateIndex[idate] = idate+1;
        dateIndex[idate] = idate-NDATE+1;
        newDeadthCase[idate] *=20.;
        if (idate<6) {
            average7day[idate]= newCase[idate];
        }
        else {
            average7day[idate] =  newCase[idate];
            for (Int_t ith=1; ith<7; ++ith) {
                average7day[idate] += newCase[idate-ith];
            }
            average7day[idate] /=7.;
        }
    }
    
    TGraph *pgr_dailyDead = new TGraph(NDATE,dateIndex, newDeadthCase);
    
    TGraph *pgr_dailyConfirm = new TGraph(NDATE,dateIndex, newCase);
    Double_t* exxcp2 =pgr_dailyConfirm->GetX();
    Double_t* eyycp2 =pgr_dailyConfirm->GetY();
    double xmincp2 = TMath::MinElement(pgr_dailyConfirm->GetN(), exxcp2);
    double xmaxcp2 = TMath::MaxElement(pgr_dailyConfirm->GetN(), exxcp2);

    TGraph *gabovecp2 = new TGraph(pgr_dailyConfirm->GetN()+3);
    for (Int_t i=0;i<pgr_dailyConfirm->GetN();i++) gabovecp2->SetPoint(i,exxcp2[i],eyycp2[i]);
    gabovecp2->SetPoint(pgr_dailyConfirm->GetN(),xmaxcp2,pgr_dailyConfirm->GetMaximum());
    gabovecp2->SetPoint(pgr_dailyConfirm->GetN()+1,xmincp2,pgr_dailyConfirm->GetMaximum());
    gabovecp2->SetPoint(pgr_dailyConfirm->GetN()+2,xmincp2,pgr_dailyConfirm->GetMinimum());
    gabovecp2->SetFillColor(15);
    
    
    TGraph *pgr_dailyConfirm_ave7 = new TGraph(NDATE,dateIndex,average7day);
    TCanvas *c1 = new TCanvas("c1","",1200,800);
    TGaxis::SetMaxDigits(3);
    TLegend* leg0 = new TLegend(.12, .65, 0.6, .89);
       leg0->SetFillStyle(0);
       leg0->SetBorderSize(0);
       leg0->SetTextSize(30);
       leg0->SetTextFont(43);
       leg0->AddEntry(pgr_dailyConfirm, "Daily cases", "f");
       leg0->AddEntry(pgr_dailyConfirm_ave7, "7-day average", "l");
        leg0->AddEntry(pgr_dailyDead, "Dead cases x20", "l");
    
    pgr_dailyConfirm->SetFillColor(15);
    
    pgr_dailyConfirm->Draw("AL");
    gabovecp2->Draw("f same");
    titleStyleGraph(pgr_dailyConfirm);
    pgr_dailyConfirm->GetYaxis()->SetTitle("Number of new Covid-19 case in Vietnam");
    pgr_dailyConfirm->GetXaxis()->SetTitle("Date countdown from 2021 Aug. 17th (=0)");
    pgr_dailyConfirm->SetTitle("");
    pgr_dailyConfirm->SetLineWidth(2);
    pgr_dailyConfirm_ave7->SetLineColor(2);
    pgr_dailyConfirm_ave7->Draw("L same");
    pgr_dailyConfirm_ave7->SetLineWidth(2);
    
    pgr_dailyDead->SetLineWidth(2);
    pgr_dailyDead->SetLineColor(4);
    pgr_dailyDead->Draw("L same");
    leg0->Draw();
    c1->RedrawAxis();
    c1->Print(Form("plots/%s_dailyconfirm.eps",savename.Data()));
    c1->Print(Form("plots/%s_dailyconfirm.png",savename.Data()));
    
    pgr_dailyConfirm->GetXaxis()->SetRangeUser(-60,1);
    c1->Update();
    c1->Print(Form("plots/%s_dailyconfirm_zoom.eps",savename.Data()));
    c1->Print(Form("plots/%s_dailyconfirm_zoom.png",savename.Data()));

}
