//#include "AliMultGlauberNBDFitter.h"
#include "multGlauberNBDFitter.h"
#include "TList.h"
#include "TFile.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TVirtualFitter.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TTree.h"

Double_t FastIntegrate(TF1 *f1, Double_t a, Double_t b, Int_t n = 5){
  //Do fast integration with N sampling points
  const Int_t nc = n;
  Double_t x[nc], y[nc];
  Double_t lWidth = (b-a)/((double)(n-1));
  for(Int_t ii=0; ii<n; ii++){
    x[ii] = a + ((double)(ii))*lWidth;
    y[ii] = f1->Eval( x[ii] );
  }
  //Now go via trapezoids, please (this probably has a name)
  Double_t lIntegral = 0;
  for(Int_t ii=0; ii<n-1; ii++){
    lIntegral += 0.5*lWidth*(y[ii]+y[ii+1]);
  }
  return lIntegral/(b-a);
}

void glauber_fit(double min_fit_range = 70, double max_fit_range = 400, bool do_fit = kTRUE){

    // 40% - mult 50 - nPart = 70
    // 50% - mult 30 - nPart = 40

    TFile *fIn_npart_ncoll = new TFile("basehistos.root", "READ");
    TH2F *hist_npart_ncoll = (TH2F*) fIn_npart_ncoll -> Get("hNpNc");

    TFile *fIn_mult = new TFile("output_na60plus.root", "READ");
    TH1F *hist_mult = (TH1F*) fIn_mult -> Get("multiplicity selected");
    hist_mult -> SetLineColor(kBlack);
    hist_mult -> SetMarkerColor(kBlack);
    hist_mult -> SetMarkerStyle(24);
    hist_mult -> Rebin(10);
    hist_mult -> SetTitle("");
    hist_mult -> GetXaxis() -> SetRangeUser(0, max_fit_range);
    
    multGlauberNBDFitter *glauberFitter = new multGlauberNBDFitter("glauberFitter");
    glauberFitter -> SetNpartNcollCorrelation(hist_npart_ncoll);
    glauberFitter -> SetInputV0M(hist_mult);
    glauberFitter -> SetAncestorMode(2);
    glauberFitter -> SetFitNpx(100000);
    TF1 *func_glauber = glauberFitter -> GetGlauberNBD();
    if (do_fit) {
        //1  p0           3.14395e-01   1.67300e-03   1.28290e-06  -9.52161e-01
        //2  p1           6.95855e+03   5.86692e+04  -5.69103e+01  -2.87883e-08
        //3  p2          -8.64267e+00   3.13768e-04  -5.93381e-06   0.00000e+00
        //4  p3           5.17488e+04   1.05429e+03   1.05429e+03   1.04134e-06

        //func_glauber -> FixParameter(0, 3.14395e-01);
        //func_glauber -> FixParameter(1, 6.95855e+03);
        //func_glauber -> FixParameter(2, -8.64267e+00);
        //func_glauber -> FixParameter(3, 5.17488e+04);

        func_glauber -> FixParameter(0, 0.7);
        func_glauber -> FixParameter(1, 20);
        func_glauber -> FixParameter(2, 1.);
        func_glauber -> SetParameter(3, 5.17488e+04);

        //glauberFitter -> SetMu(1);
        //glauberFitter -> Setk(1.5);
        //glauberFitter -> Setf(0.8);
        //glauberFitter -> SetNorm(5.17488e+04);



        //glauberFitter -> SetMu(3.14395e-01);
        //glauberFitter -> Setk(6.95855e+03);
        //glauberFitter -> Setf(-8.64267e+00);
        //glauberFitter -> SetNorm(5.17488e+04);

        glauberFitter -> SetFitNpx(100000);
        glauberFitter -> SetFitRange(min_fit_range, max_fit_range);
        glauberFitter -> SetFitOptions("REM0");

        glauberFitter -> InitializeNpNc();
        glauberFitter -> InitAncestor();

        int fit_status = glauberFitter -> DoFit();
        int attempts = 1;
        while(attempts < 10 && fit_status == 0){
            std::cout << "Attempting fit again (" << attempts << " attempt)..." << std::endl;
            fit_status = glauberFitter -> DoFit();
        }
        std::cout << "Final fit status: " << fit_status << std::endl;
        std::cout << "OK!" << std::endl;
    } else {
        //glauberFitter -> SetMu(3.14395e-01);
        //glauberFitter -> Setk(6.95855e+03);
        //glauberFitter -> Setf(-8.64267e+00);
        //glauberFitter -> SetNorm(5.17488e+04);
        //glauberFitter -> SetFitNpx(100000);
    }


    TF1 *fNBD = (TF1*) glauberFitter -> GetNBD();
    TF1 *fGlauberNBD = (TF1*) glauberFitter -> GetGlauberNBD();

    TH1D *hist_ratio = (TH1D*) hist_mult -> Clone("hist_ratio");
    hist_ratio -> SetTitle("");

    // Evaluate the Fit / Data ratio
    for(int iBin = 1;iBin < hist_ratio -> GetNbinsX()+1;iBin++) {
        double ratio = hist_ratio -> GetBinContent(iBin);
        double ratio_integral = FastIntegrate(fGlauberNBD, hist_ratio -> GetBinLowEdge(iBin), hist_ratio -> GetBinLowEdge(iBin+1), 10);
        ratio /= ratio_integral;
        hist_ratio -> SetBinContent(iBin, ratio);
        hist_ratio -> SetBinError(iBin, hist_ratio -> GetBinError(iBin) / ratio_integral);
    }

    Double_t lXValue[80000];
    Double_t lYValue[80000];
    lXValue[0] = 0.0;
    lYValue[0] = 0.0;
    const Long_t lCutoff=5000;
    
    Double_t lIntegral = 0.0;
    Double_t correction_factor = 0;
    double mult_bins[370];
    for(Int_t ii=1; ii<lCutoff; ii++){
        Double_t lXlower = ((Double_t)(ii-1));
        Double_t lX = ((Double_t)(ii));
        
        //Calculate function delta and sum
        Double_t lFunctionDelta = FastIntegrate( fGlauberNBD, lXlower, lX, 20)/hist_mult->GetBinWidth(1);
        lIntegral+=lFunctionDelta;
        
        //calculate histogram content above
        Int_t lBinNumber = hist_mult -> FindBin(lX+1e-6);
        Double_t lRemaining = hist_mult->Integral(lBinNumber,hist_mult->GetNbinsX());


    
        //Calculate fraction elapsed of total function
        Double_t lFractionElapsed = lIntegral/(lIntegral+lRemaining);
        
        lXValue[ii] = lX;
        lYValue[ii] = lFractionElapsed;
        
        if (ii < 370) {
            mult_bins[ii] = lIntegral;
            //std::cout << correction_factor << std::endl;
            //std::cout << "centrality = " << lIntegral/(lIntegral + correction_factor) << std::endl;
            //cout<<"At ii = "<<ii<<", function delta =  "<<lFunctionDelta<<" int = "<<lIntegral<<", histo remaining: "<<lRemaining<<" fraction elapsed: "<<lFractionElapsed<<endl;
            //if(ii%500==0) cout<<"At ii = "<<ii<<", function delta =  "<<lFunctionDelta<<" int = "<<lIntegral<<", histo remaining: "<<lRemaining<<" fraction elapsed: "<<lFractionElapsed<<endl;
        }
    }
    
    double mult_selected[4];
    int  counter = 1;
    mult_selected[0] = 61;
    double mult_range = (mult_bins[360] - mult_bins[61]) / 4.; // 40%
    double tmp_mult = mult_bins[61];
    double delta_mult = 0;
    for (int i = 62;i < 370;i++) {
        if (mult_bins[i] - tmp_mult > mult_range) {
            tmp_mult = mult_bins[i];
            std::cout << tmp_mult << " -> " << i << std::endl;
            mult_selected[counter] = i;
            counter++;
        }
    }

    TLine *lineAPsVert[4];
    TLine *lineAPsVert2[4];

    for(int iMult = 0;iMult < 4;iMult++){ 
        lineAPsVert[iMult] = new TLine( mult_selected[iMult], 0, mult_selected[iMult], 2);
        lineAPsVert[iMult]->SetLineColor(kBlack);
        lineAPsVert[iMult]->SetLineWidth(1);
        lineAPsVert[iMult]->SetLineStyle(9);

        lineAPsVert2[iMult] = new TLine(mult_selected[iMult], 0, mult_selected[iMult], fGlauberNBD -> Eval(mult_selected[iMult]));
        lineAPsVert2[iMult]->SetLineColor(kBlack);
        lineAPsVert2[iMult]->SetLineWidth(1);
        lineAPsVert2[iMult]->SetLineStyle(9);
    }

    ////////////////
    /*
    TGraph *gr = new TGraph();
    TGraph *gra_inverse = new TGraph();
    gra_inverse -> SetMarkerStyle(20);
    gra_inverse -> SetMarkerColor(kBlack);
    for(Int_t ii=1; ii<lCutoff; ii++){
        gr->SetPoint(ii, lXValue[ii], lYValue[ii]);
        gra_inverse->SetPoint(ii, lYValue[ii], lXValue[ii]);
    }

    TCanvas *canvas_inverse = new TCanvas("canvas_inverse", "", 800, 600);
    TH2D *hist_grid_inverse = new TH2D("hist_grid_inverse", "", 1000, 0, 1, 1000, 0, 1000);
    hist_grid_inverse -> Draw();
    gra_inverse -> Draw("EPsame");


    Float_t lDrawAPs[] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    Float_t lRawVals[7];
    TLine *lineAPsVert[7];
    TLine *lineAPsHori[7];
    TLine *lineAPsVert2[7];
    TGraph *grAPs = new TGraph();
    TGraph *grAPsSolidCore = new TGraph();
    
    for(Int_t ii=0; ii<7; ii++){
        lRawVals[ii] = gra_inverse->Eval(lDrawAPs[ii]);
        cout<<"AP point: "<<lRawVals[ii]<<" - "<<lDrawAPs[ii]<<endl;
        lineAPsVert[ii] = new TLine( lRawVals[ii], 0, lRawVals[ii], lDrawAPs[ii]);
        lineAPsVert[ii]->SetLineColor(kBlack);
        lineAPsVert[ii]->SetLineWidth(1);
        lineAPsVert[ii]->SetLineStyle(9);
        lineAPsVert[ii]->Draw();
        lineAPsHori[ii] = new TLine( 0, lDrawAPs[ii], lRawVals[ii], lDrawAPs[ii]);
        lineAPsHori[ii]->SetLineColor(kBlack);
        lineAPsHori[ii]->SetLineWidth(1);
        lineAPsHori[ii]->SetLineStyle(9);
        lineAPsHori[ii]->Draw();
        grAPs->SetPoint(ii, lRawVals[ii], lDrawAPs[ii]);
        grAPsSolidCore->SetPoint(ii, lRawVals[ii], lDrawAPs[ii]);
        
        c1->cd(2);

        lineAPsVert[ii] = new TLine( lRawVals[ii], 0, lRawVals[ii], 2);
        lineAPsVert[ii]->SetLineColor(kBlack);
        lineAPsVert[ii]->SetLineWidth(1);
        lineAPsVert[ii]->SetLineStyle(9);

        lineAPsVert2[ii] = new TLine(lRawVals[ii], 0, lRawVals[ii], fGlauberNBD -> Eval(lRawVals[ii]));
        lineAPsVert2[ii]->SetLineColor(kBlack);
        lineAPsVert2[ii]->SetLineWidth(1);
        lineAPsVert2[ii]->SetLineStyle(9);

        if (lDrawAPs[ii] == 0.6 || lDrawAPs[ii] == 0.5) {
            lineAPsVert[ii]->SetLineColor(kRed);
            lineAPsVert2[ii]->SetLineColor(kRed);
        }
        
        c1->cd(1);
    }
    */
    ////////////////

    TLatex *latex_title = new TLatex();
    latex_title -> SetNDC();
    latex_title -> SetTextSize(0.042);
    float lPosText = 0.76;
    float lYShift = 0.25;

    TLine *line_unity = new TLine(min_fit_range, 1., 400, 1.);
    line_unity -> SetLineColor(kGray+1);
    line_unity -> SetLineWidth(2);
    line_unity -> SetLineStyle(kDashed);

    gStyle -> SetOptStat(0);
    TCanvas *canvas_fit = new TCanvas("canvas_fit ", "", 1300, 900);
    canvas_fit -> SetFrameFillStyle(0);
    canvas_fit -> SetFillStyle(0);
    canvas_fit -> Divide(1,2);
    canvas_fit -> cd(1) -> SetFrameFillStyle(0);
    canvas_fit -> cd(1) -> SetFillStyle(0);
    canvas_fit -> cd(2) -> SetFrameFillStyle(0);
    canvas_fit -> cd(2) -> SetFillStyle(0);
    
    canvas_fit -> cd(1);
    canvas_fit -> cd(1) -> SetLogy();
    canvas_fit -> cd(1) -> SetTicks(1,1);
    canvas_fit -> cd(1) -> SetPad(0,0.5,1,1);
    canvas_fit -> cd(2) -> SetPad(0,0.0,1,.5);
    
    canvas_fit -> cd(1) -> SetBottomMargin(0.001);
    canvas_fit -> cd(1) -> SetRightMargin(0.25);
    canvas_fit -> cd(1) -> SetTopMargin(0.02);
    canvas_fit -> cd(1) -> SetLeftMargin(0.07);
    
    canvas_fit -> cd(2) -> SetBottomMargin(0.14);
    canvas_fit -> cd(2) -> SetRightMargin(0.25);
    canvas_fit -> cd(2) -> SetTopMargin(0.001);
    canvas_fit -> cd(2) -> SetLeftMargin(0.07);
    canvas_fit -> cd(2) -> SetTicks(1,1);

    canvas_fit -> cd(1);
    gPad -> SetLogy(1);
    hist_mult -> Draw("EP");
    for (int i = 0;i < 4;i++) {
        lineAPsVert2[i]->Draw("same");
    }
    fGlauberNBD -> Draw("same");
    latex_title -> DrawLatex(lPosText, 0.67+lYShift, "Pb-Pb 150 GeV");
    latex_title -> DrawLatex(lPosText, 0.61+lYShift, "Glauber + NBD fit");
    latex_title -> DrawLatex(lPosText, 0.55+lYShift, "Test beam November 2022");
    latex_title -> SetTextFont(42);
    latex_title -> DrawLatex(lPosText, 0.49+lYShift, Form("Fit range: %.1f-%.1f", min_fit_range, max_fit_range));
    latex_title -> DrawLatex(lPosText, 0.43+lYShift, Form("#Chi^{2}/ndf: %.1f / %i = %.3f", fGlauberNBD -> GetChisquare(), fGlauberNBD -> GetNDF(), fGlauberNBD -> GetChisquare() / ((Double_t)(fGlauberNBD -> GetNDF()))));
    latex_title -> DrawLatex(lPosText, 0.31+lYShift, Form("Glauber f: %.3f", fGlauberNBD -> GetParameter(2)));
    latex_title -> DrawLatex(lPosText, 0.25+lYShift, Form("Glauber #mu: %.3f", fGlauberNBD -> GetParameter(0)));
    latex_title -> DrawLatex(lPosText, 0.19+lYShift, Form("Glauber k: %.3f", fGlauberNBD -> GetParameter(1)));

    latex_title -> DrawLatex(lPosText, 0.38-0.06, Form("40%% anchor point: %.0f", mult_selected[0]));
    latex_title -> DrawLatex(lPosText, 0.38-0.12, Form("30%% anchor point: %.0f", mult_selected[1]));
    latex_title -> DrawLatex(lPosText, 0.38-0.18, Form("20%% anchor point: %.0f", mult_selected[2]));
    latex_title -> DrawLatex(lPosText, 0.38-0.24, Form("10%% anchor point: %.0f", mult_selected[3]));

    canvas_fit -> cd(2);
    hist_ratio -> GetYaxis() -> SetTitle("Data/Fit");
    hist_ratio -> GetXaxis() -> SetTitle("Track multiplicity");
    hist_ratio -> GetYaxis() -> SetTitleSize(0.055);
    hist_ratio -> GetYaxis() -> SetTitleOffset(0.7);
    hist_ratio -> GetXaxis() -> SetTitleSize(0.055);
    hist_ratio -> GetYaxis() -> SetLabelSize(0.045);
    hist_ratio -> GetXaxis() -> SetLabelSize(0.045);
    hist_ratio -> GetYaxis() -> SetRangeUser(0, 2);
    hist_ratio -> GetXaxis() -> SetRangeUser(0, max_fit_range);
    hist_ratio -> SetMarkerStyle(20);
    hist_ratio -> SetMarkerColor(kBlack);
    hist_ratio -> SetLineColor(kBlack);
    hist_ratio -> SetMarkerSize(.7);
    hist_ratio -> Draw("EP");
    for (int i = 0;i < 4;i++) {
        lineAPsVert[i]->Draw("same");
    }
    line_unity -> Draw("same");

    canvas_fit -> SaveAs("glauber_fit_test_beam_november.pdf");


    TFile *fOut = new TFile("myFirstFit.root", "RECREATE");
    fNBD -> Write();
    fGlauberNBD -> Write();
    hist_mult -> Write();
    hist_npart_ncoll -> Write();
    canvas_fit -> Write("canvasFit");
    fOut -> Close();
}

void plot_results(double min_fit_range = 70){
    TFile *fIn = new TFile("myFirstFit.root", "READ");
    TF1 *fGlauberNBD = (TF1*) fIn -> Get("fGlauberNBD");
    TH1F *hist_mult = (TH1F*) fIn -> Get("multiplicity selected");

    return;

    TLatex *latex_title = new TLatex();
    latex_title -> SetNDC();
    latex_title -> SetTextSize(0.042);
    float lPosText = 0.76;
    float lYShift = 0.25;

    TLine *line_unity = new TLine(min_fit_range, 1., 400, 1.);
    line_unity -> SetLineColor(kGray+1);
    line_unity -> SetLineWidth(2);
    line_unity -> SetLineStyle(kDashed);

    gStyle -> SetOptStat(0);
    TCanvas *canvas_fit = new TCanvas("canvas_fit ", "", 1300, 900);
    canvas_fit -> SetFrameFillStyle(0);
    canvas_fit -> SetFillStyle(0);
    canvas_fit -> Divide(1,2);
    canvas_fit -> cd(1) -> SetFrameFillStyle(0);
    canvas_fit -> cd(1) -> SetFillStyle(0);
    canvas_fit -> cd(2) -> SetFrameFillStyle(0);
    canvas_fit -> cd(2) -> SetFillStyle(0);
    
    canvas_fit -> cd(1);
    canvas_fit -> cd(1) -> SetLogy();
    canvas_fit -> cd(1) -> SetTicks(1,1);
    canvas_fit -> cd(1) -> SetPad(0,0.5,1,1);
    canvas_fit -> cd(2) -> SetPad(0,0.0,1,.5);
    
    canvas_fit -> cd(1) -> SetBottomMargin(0.001);
    canvas_fit -> cd(1) -> SetRightMargin(0.25);
    canvas_fit -> cd(1) -> SetTopMargin(0.02);
    canvas_fit -> cd(1) -> SetLeftMargin(0.07);
    
    canvas_fit -> cd(2) -> SetBottomMargin(0.14);
    canvas_fit -> cd(2) -> SetRightMargin(0.25);
    canvas_fit -> cd(2) -> SetTopMargin(0.001);
    canvas_fit -> cd(2) -> SetLeftMargin(0.07);
    canvas_fit -> cd(2) -> SetTicks(1,1);

    canvas_fit -> cd(1);
    gPad -> SetLogy(1);
    hist_mult -> Draw("EP");
    fGlauberNBD -> Draw("same");
}

void set_Npart_Ncoll_correlation(){
    //cout<<"Loading library..."<<endl;
    //gSystem->Load("libmultGlauberNBDFitter.dylib");

    TFile *fIn = new TFile("GlauberMCSPS_Pb_Pb_0.00_16.00.root", "READ");
    fIn -> ls();
    TTree *tree = (TTree*) fIn -> Get("tglauber");
    tree -> Show();
    int npart, ncoll;
    tree -> SetBranchAddress("npart", &npart);
    tree -> SetBranchAddress("ncoll", &ncoll);

    // Test TH2F / TH2D
    TH2F *hist_npart_ncoll = new TH2F("hNpNc", "", 1000, -0.5, 999.5, 1000, -0.5, 999.5);
    hist_npart_ncoll -> GetXaxis() -> SetTitle("#it{N}_{part}");
    hist_npart_ncoll -> GetYaxis() -> SetTitle("#it{N}_{coll}");
    TH1F *hist_npart = new TH1F("hNp", "", 1000, -0.5, 999.5);
    hist_npart -> GetXaxis() -> SetTitle("#it{N}_{part}");
    hist_npart -> GetYaxis() -> SetTitle("counts");
    hist_npart -> SetMarkerColor(kBlack);
    hist_npart -> SetMarkerStyle(20);

    //for (int i = 0;i < 1000;i++) {
        //std::cout << hist_npart_ncoll -> GetXaxis() -> GetBinCenter(i+1) << std::endl; 
    //}


    for (int i = 0;i < tree -> GetEntries();i++) {
        tree -> GetEntry(i);
        //std::cout << npart << " " << ncoll << std::endl;
        hist_npart_ncoll -> Fill(npart, ncoll);
        hist_npart -> Fill(npart);
    }

    double integral_npart = hist_npart -> GetEntries();
    double integral = 0;
    for (int iBin = 0;iBin < 500;iBin++) {
        integral += hist_npart -> GetBinContent(iBin+1);
        if (integral/integral_npart > 0.499 && integral/integral_npart < 0.502) {
            std::cout << integral << " -> " << integral/integral_npart << " ==> " << hist_npart -> GetBinCenter(iBin+1) << std::endl;
            std::cout << "-----------> n Trk = " << (hist_npart -> GetBinCenter(iBin+1) * 360.) / 416. << std::endl;
        }
        if (integral/integral_npart > 0.599 && integral/integral_npart < 0.602) {
            std::cout << integral << " -> " << integral/integral_npart << " ==> " << hist_npart -> GetBinCenter(iBin+1) << std::endl;
            std::cout << "-----------> n Trk = " << (hist_npart -> GetBinCenter(iBin+1) * 360.) / 416. << std::endl;
        }
    }

    TFile *fOut = new TFile("basehistos.root", "RECREATE");
    hist_npart_ncoll -> Write();
    hist_npart -> Write();
    fOut -> Close();
}