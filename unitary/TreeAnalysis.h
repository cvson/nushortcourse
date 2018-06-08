//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 15 18:55:14 2013 by ROOT version 5.34/04
// from TTree tree/tree
// found on file: EventToy.root
//////////////////////////////////////////////////////////

#ifndef TreeAnalysis_h
#define TreeAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TreeAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        sinsq12;
   Double_t        sinsq23;
   Double_t        sinsq13;
   Double_t        cpphase;

   // List of branches
   TBranch        *b_sinsq12;   //!
   TBranch        *b_sinsq23;   //!
   TBranch        *b_sinsq13;   //!
   TBranch        *b_cpphase;   //!

   TreeAnalysis(TTree *tree=0);
   virtual ~TreeAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TreeAnalysis_cxx
TreeAnalysis::TreeAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("EventToy.root");
       //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("EventToy_normal.root");//normal
       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("EventToy_invert.root");//normal
      if (!f || !f->IsOpen()) {
         //f = new TFile("EventToy_normal.root");
          f = new TFile("EventToy_invert.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

TreeAnalysis::~TreeAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("sinsq12", &sinsq12, &b_sinsq12);
   fChain->SetBranchAddress("sinsq23", &sinsq23, &b_sinsq23);
   fChain->SetBranchAddress("sinsq13", &sinsq13, &b_sinsq13);
   fChain->SetBranchAddress("cpphase", &cpphase, &b_cpphase);
   Notify();
}

Bool_t TreeAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeAnalysis_cxx
