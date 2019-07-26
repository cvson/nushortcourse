//
//  plot_pmns_triangle.C
////////////////////////////////////////////////////
//
//  to plot PMNS triangle
//
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
void plot_pmns_triangle_2019(){
    bool isNormalHierarchy=true;
	//nufit 4.1
    double sinsq12=0.310;
    
    double sinsq23N=0.580;
    double sinsq23I=0.584;
    double sinsq23 = isNormalHierarchy?sinsq23N:sinsq23I;
    
    double sinsq13N=0.02241;
    double sinsq13I=0.02264;
    double sinsq13 = isNormalHierarchy?sinsq13N:sinsq13I;
    
    double deltaN=215.0*TMath::Pi()/180.;
    double deltaI=284.0*TMath::Pi()/180.;
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
    double LengthDP = 200;
    double Ue1=c12*c13;
    cout<<"Ue1 "<<Ue1<<" amp2 "<<pow(Ue1,2)*LengthDP<<endl;
    double Ue2=s12*c13;
    cout<<"Ue2 "<<Ue2<<" amp2 "<<pow(Ue2,2)*LengthDP<<endl;
    TComplex Ue3=cpPhaseNeg*s13;
    
    //muon
    cout<<"Ue3 Re "<<Ue3.Re()<<" Im "<<Ue3.Im()<<" amp2 "<<Ue3.Rho2()*LengthDP<<endl;
    TComplex Umu1=-s12*c23-c12*s23*s13*cpPhasePos;
    cout<<"Umu1 Re "<<Umu1.Re()<<" Im "<<Umu1.Im()<<" amp2 "<<Umu1.Rho2()*LengthDP<<endl;
    TComplex Umu2=c12*c23-s12*s23*s13*cpPhasePos;
    cout<<"Umu2 Re "<<Umu2.Re()<<" Im "<<Umu2.Im()<<" amp2 "<<Umu2.Rho2()*LengthDP<<endl;
    double Umu3=s23*c13;
    cout<<"Umu3  "<<Umu3<<" amp2 "<<pow(Umu3,2)*LengthDP<<endl;
    
    //tau
    TComplex Utau1=s12*s23-c12*c23*s13*cpPhasePos;
    cout<<"Utau1 Re "<<Utau1.Re()<<" Im "<<Utau1.Im()<<" amp2 "<<Utau1.Rho2()*LengthDP<<endl;
    TComplex Utau2=-c12*s23-s12*c23*s13*cpPhasePos;
    cout<<"Utau2 Re "<<Utau2.Re()<<" Im "<<Utau2.Im()<<" amp2 "<<Utau2.Rho2()*LengthDP<<endl;
    double Utau3=c23*c13;
    cout<<"Utau3  "<<Utau3<<" amp2 "<<pow(Utau3,2)*LengthDP<<endl;
    
    //test unitary condition
    double norm_e=Ue1**2+Ue2**2+Ue3.Rho2();
    double norm_mu=Umu1.Rho2()+Umu2.Rho2()+Umu3**2;
    double norm_tau=Utau1.Rho2()+Utau2.Rho2()+Utau3**2;
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
    
    cout<<"sideA Re "<<sideA.Re()<<" Im "<<sideA.Im()<<" Amp "<<sideA.Rho()<<" angle "<<180*(1-sideA.Theta()/TMath::Pi())<<endl;
    cout<<"sideB Re "<<sideB.Re()<<" Im "<<sideB.Im()<<" Amp "<<sideB.Rho()<<" angle "<<180*(1+sideB.Theta()/TMath::Pi())<<endl;
    
    //Draw a triangle;
    c1 = new TCanvas("c1");
    gStyle->SetOptStat(0);
    c1->Range(-1.0,-1.0,1.0,1.0);
    TH2D *h2 = new TH2D("h2","",200,-0.5,1.2,200, -0.3,0.3);
    h2->GetXaxis()->SetTitle("#rho (Real)");
    h2->GetYaxis()->SetTitle("#eta (Imagine)");
    h2->GetYaxis()->SetTitleOffset(1.2);
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->CenterTitle();
    gPad->SetGrid();
    TLine *linec= new TLine(0,0,1,0);
    TLine *linea= new TLine(1,0,1+sideA.Re(),+sideA.Im());
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
    TLatex texsideA(-0.2,0.08,"#left|#frac{U_{e1}^{*}U_{#mu1}}{U_{e2}^{*}U_{#mu2}}#right|");
    TLatex texsideB(0.8,0.08,"#left|#frac{U_{e3}^{*}U_{#mu3}}{U_{e2}^{*}U_{#mu2}}#right|");
    
    TLatex vertexC(0.22,0.22,Form("C(%.2f,%.2f)",1.0+sideA.Re(),sideA.Im()));
    TLatex angleA(0.6,0.005,Form("#gamma =%.2f^{#circ}",180*(1-sideA.Theta()/TMath::Pi())));
    TLatex angleB(0.1,0.01,Form("#beta =%.2f^{#circ}",180*(1+sideB.Theta()/TMath::Pi())));
    angleA.SetTextColor(kAzure);
    angleB.SetTextColor(kAzure);
    vertexA.SetTextColor(95);
    vertexB.SetTextColor(95);
    vertexC.SetTextColor(95);
    vertexA.Draw("same");
    vertexB.Draw("same");
    texsideA.Draw("same");
    texsideB.Draw("same");
    TLatex jcp(-0.4,-0.22,Form("J_{CP} = sign(#delta_{CP})x 2 x Area x |U_{e2}^{*}U_{#mu2}|^{2} = %.3f",-1.*sideA.Im()*pow(Ue2Umu2.Rho(),2)));
    jcp.Draw("same");
    if(isNormalHierarchy){
        vertexC.Draw("same");
        angleA.Draw("same");
        angleB.Draw("same");
    }
    else {
        vertexC.DrawLatex(0.22,0.07,Form("C(%.2f,%.2f)",1.0+sideA.Re(),sideA.Im()));
        angleA.DrawLatex(0.1,0.01,Form("#gamma =%.2f^{#circ}",sideA.Theta()));
        angleB.DrawLatex(0.6,0.005,Form("#beta =%.2f^{#circ}",sideB.Theta()));
    }
    if(isNormalHierarchy)c1->Print("triangle_normal.eps");
    else c1->Print("triangle_invert.eps");
}
