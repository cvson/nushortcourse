//
//  plot_pmns_triangle.C
////////////////////////////////////////////////////
//
//  To read cross-section file and plot x-sec for different channels
//
///////////////////////////////////////////////////
//  Created by S. Cao, cvson@utexas.edu
void plot_pmns_triangle(){
    bool isNormalHierarchy=true;
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
}
