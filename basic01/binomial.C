void binomial()
{
    
  Double_t pvalue(0) ;

    pvalue += TMath::BinomialI(6*pow(10,-4),191,5)*5*pow(10,6);

  cout << "p-value(N>=7|s=0) = " << pvalue << endl;


}


