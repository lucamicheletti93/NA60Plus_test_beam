import matplotlib.pyplot as plt
import array as arr
import numpy as np
import os
import sys
import argparse
import ROOT
from os import path
from ROOT import TCanvas, TLatex, TF1, TFile, TPaveText, TMath, TH1F, TString, TLegend, TRatioPlot, TGaxis
from ROOT import gROOT, gBenchmark, gPad, gStyle, kTRUE, kFALSE

def SetLatex(latex):
    latex.SetTextSize(0.035)
    latex.SetNDC()
    latex.SetTextFont(42)

def SetLegend(legend):
    legend.SetBorderSize(0)
    legend.SetFillColor(10)
    legend.SetFillStyle(1)
    legend.SetLineStyle(0)
    legend.SetLineColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)

def SetHistStyle(hist):
    hist.SetTitle("")
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitle("Entries")
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(0.05)

def SetHistStyle2(hist, title, xaxis, yaxis, color, linewidth, fillstyle, trasparency):
    hist.SetTitle(title)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitle(xaxis)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitle(yaxis)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetLineWidth(linewidth)
    #hist.SetFillStyle(fillstyle)
    hist.SetFillColorAlpha(color, trasparency)

def SetHistStyle3(hist, title, xaxis, yaxis, color, linewidth, fillstyle, trasparency):
    hist.SetTitle(title)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitle(xaxis)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitle(yaxis)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetLineWidth(linewidth)

def SetHistStyle4(hist, title, xaxis, yaxis, color, linewidth):
    hist.SetTitle(title)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitle(xaxis)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitle(yaxis)
    hist.SetLineColor(color)
    hist.SetLineWidth(linewidth)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(20)

def LoadStyle():
    gStyle.SetOptStat(0)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetEndErrorSize(0.0)
    gStyle.SetTitleSize(0.05,"X")
    gStyle.SetTitleSize(0.045,"Y")
    gStyle.SetLabelSize(0.045,"X")
    gStyle.SetLabelSize(0.045,"Y")
    gStyle.SetTitleOffset(1.2,"X")
    gStyle.SetTitleOffset(1.35,"Y")