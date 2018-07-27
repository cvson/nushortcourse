{
TStyle* minosStyle = new  TStyle("minosStyle", "MINOS Style");

// Colors

//set the background color to white
minosStyle->SetFillColor(10);
minosStyle->SetFrameFillColor(10);
minosStyle->SetCanvasColor(10);
minosStyle->SetPadColor(10);
minosStyle->SetTitleFillColor(0);
minosStyle->SetStatColor(10);

//dont put a colored frame around the plots
minosStyle->SetFrameBorderMode(0);
minosStyle->SetCanvasBorderMode(0);
minosStyle->SetPadBorderMode(0);

//use the primary color palette
minosStyle->SetPalette(1,0);

//set the default line color for a histogram to be black
minosStyle->SetHistLineColor(kBlack);

//set the default line color for a fit function to be red
minosStyle->SetFuncColor(kRed);

//make the axis labels black
minosStyle->SetLabelColor(kBlack,"xyz");

//set the default title color to be black
minosStyle->SetTitleColor(kBlack);
 
// Sizes

//set the margins
minosStyle->SetPadBottomMargin(0.2);
minosStyle->SetPadTopMargin(0.075);
minosStyle->SetPadLeftMargin(0.15);

//set axis label and title text sizes
minosStyle->SetLabelSize(0.04,"xyz");
minosStyle->SetTitleSize(0.04,"xyz");
minosStyle->SetTitleOffset(1.2,"x");
minosStyle->SetTitleOffset(1.3,"yz");
minosStyle->SetStatFontSize(0.04);
minosStyle->SetTextSize(0.04);
minosStyle->SetTitleBorderSize(0);
minosStyle->SetStatBorderSize(0);
 
//set line widths
minosStyle->SetHistLineWidth(2);
minosStyle->SetFrameLineWidth(2);
minosStyle->SetFuncWidth(2);

// Misc

//align the titles to be centered
minosStyle->SetTitleAlign(22);

//set the number of divisions to show
minosStyle->SetNdivisions(506, "xy");

//turn off xy grids
minosStyle->SetPadGridX(0);
minosStyle->SetPadGridY(0);

//set the tick mark style
minosStyle->SetPadTickX(1);
minosStyle->SetPadTickY(1);

//show the fit parameters in a box
minosStyle->SetOptFit(1111);

//turn off all other stats
//minosStyle->SetOptStat(1111);
minosStyle->SetOptStat(0000);

//marker settings
minosStyle->SetMarkerStyle(8);
minosStyle->SetMarkerSize(0.7);

// Fonts
const int kMinosFont = 42;

minosStyle->SetStatFont(kMinosFont);
minosStyle->SetLabelFont(kMinosFont,"xyz");
minosStyle->SetTitleFont(kMinosFont,"xyz");
minosStyle->SetTextFont(kMinosFont);

//Set Legend Border Size
minosStyle->SetLegendBorderSize(0);

// Get moodier colours
const Int_t NRGBs = 5;
const Int_t NCont = 255;
Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
minosStyle->SetNumberContours(NCont);

//done
minosStyle->cd();

gROOT->ForceStyle();
gStyle->ls();

}

