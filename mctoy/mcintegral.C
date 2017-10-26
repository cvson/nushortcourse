//
//  mcintegral.C
////////////////////////////////////////////////////
//
//  Simple Monte Carlo for integral
//
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
void mcintegral(){
    //this is to produce gif animation
    TCanvas *c1 = new TCanvas("c1","pi",800,800);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gBenchmark->Start("hmc");
    gSystem->Unlink("hmcintegral.gif"); // delete old file

    //emty histogram
    TH2* hempty = new TH2F("","",1000,-1,1,1000,-1,1);
    hempty->Draw();

    //condition to print
    const Int_t kUPDATE = 500;
    const Int_t kIteration = 100000;//1000000;
    const double PI = TMath::Pi();
    double rmax = 1.0;
    TArc *parcout = new TArc(0,0,rmax,0,360);
    parcout->SetLineWidth(1);
    
    
    double rmin = 0.98;
    TArc *parcin = new TArc(0,0,rmin,0,360);
    parcin->SetLineWidth(1);
    parcout->Draw();
    parcin->Draw();
    
    double z_level = -0.6;
    TLine *pline = new TLine(-sqrt((pow(rmin,2)-pow(z_level,2))),z_level,sqrt((pow(rmin,2)-pow(z_level,2))), z_level);
    pline->Draw();
    double rhoshell = 2.7;
    double rhowater = 1;
    //theoretical
    //shell masss
    double theomassshell=(pow(rmax,3)-pow(rmin,3))*4*PI/3.0*rhoshell;
    
    //Liquid mass
	double theomassliq=rhowater * (pow(rmin,3)*2*PI/3.0+PI*z_level*pow(rmin,2)-PI/3.0*pow(z_level,3));
    
    //central mass of liquid
	double theozliq=rhowater * (PI/4*(2*pow(rmin*z_level,2)-pow(z_level,4)-pow(rmin,4)));
	theozliq /= theomassliq;
    
    //central mass of whole object, since central mass of shell is at zero
	double theoz=theozliq*(theomassliq/(theomassliq+theomassshell));
    
    //total detection volumne in Catersian coordinates
    //cubic with edge of 2*rmax
    double voltot = pow(2*rmax,3);
    //initialize mass
    double mass = 0;
    //initialize central mass postion
    double xcm = 0;
    double ycm = 0;
    double zcm = 0;
    
    //generate number of points with mass and sum all points
    Int_t ci;
    TMarker *mmarker;
    TString resultString;
    TLatex *   tex;
    for (int i=0; i<kIteration; i++){
        //if (i%1000==0) printf("Processing %d/%d\n",i,kIteration);
        //randome coordinate
        double ix = gRandom->Uniform(-1.0,1.0)*rmax;
        double iy = gRandom->Uniform(-1.0,1.0)*rmax;
        double iz = gRandom->Uniform(-1.0,1.0)*rmax;
        //get radius
        double ir = sqrt(ix*ix + iy*iy + iz*iz);
        //ignore if outside of sphere
        if  (ir > rmax) continue;
        //check if inside shell
        if (ir >= rmin){
            mass += rhoshell;//each point have mass equal to density
            xcm  += rhoshell * ix;
            ycm  += rhoshell * iy;
            zcm  += rhoshell * iz;
            //UT orange 180,95,4
            mmarker = new TMarker(ix, iz, 20);
            mmarker->SetMarkerSize(0.5);
            ci = TColor::GetColor("#B45F04");
            mmarker->SetMarkerColor(ci);
            mmarker->DrawClone("same");

        }
        //inside the water
        else if (iz <= z_level){
            mass += rhowater;//each point have mass equal to density
            xcm  += rhowater * ix;
            ycm  += rhowater * iy;
            zcm  += rhowater * iz;
            //navy blue 11,56,97
            mmarker = new TMarker(ix, iz, 20);
            mmarker->SetMarkerSize(0.5);
            ci = TColor::GetColor("#0B3861");
            mmarker->SetMarkerColor(ci);
            mmarker->DrawClone("same");

            
        }
        if (i==1 || ((i+1)%kUPDATE) == 0) {
            c1->Update();
            double estimatezcm = zcm/mass;
            double errorPercentage=100*TMath::Abs((estimatezcm-theoz))/TMath::Abs(theoz);
            resultString.Form("z_{mc} = %.2g,z_{theo}= %.2g, %.2g (%%) error n=%d", estimatezcm, theoz, errorPercentage,i+1);
            tex = new TLatex(0.25,0.973,resultString);
            tex->SetNDC();
            tex->SetTextAlign(13);
            tex->SetTextFont(42);
            tex->SetTextSize(0.04);
            tex->SetLineWidth(2);
            tex->Draw();
            
            if (gROOT->IsBatch()) {
                c1->Print("hmcintegral.gif+");
                printf("i = %d\n", i+1);
            }
            tex->Delete();
        }
    }//end for
    cout << "Monte-carlo results drawn Cartesian : " << endl;
	cout << "mass : " << voltot*mass/kIteration << " theoretical " << theomassliq+theomassshell << endl;
	cout << " xcm : " << xcm/mass << endl;
	cout << " ycm : " << ycm/mass << endl;
	cout << " zcm : " << zcm/mass << " theoretical " << theoz << endl;

    if (gROOT->IsBatch()) c1->Print("hmcintegral.gif++");
        gBenchmark->Show("hmc");

        
}
