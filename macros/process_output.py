from ast import parse
from itertools import count
from time import process_time_ns
from tkinter import Canvas
from turtle import color, right
import matplotlib.pyplot as plt
import array as arr
import numpy as np
import pandas as pd
from pandas.plotting import scatter_matrix
import os
import sys
import argparse
import ROOT
from os import path
from ROOT import TCanvas, TFile, TTree, TTreeReader, TH1F, TH2F, TLatex, TLegend, TF1
from ROOT import gStyle, gPad, kRed, kBlue, kGreen, kOrange, kAzure, kBlack
from ctypes import cdll
sys.path.append('../utils')
from plot_library import LoadStyle, SetLegend, SetHistStyle, SetLatex, SetHistStyle2

###
def read_data(fInName):
    '''
    Function to read the reduced tables and produce the output for ML training
    '''
    gStyle.SetOptStat(0)
    LoadStyle()
    ROOT.TGaxis.SetMaxDigits(3)

    vbbs = ["0V", "1V", "3V"]
    ithr = "Ithr_240"
    plotLimit = 3
    hists_0 = []
    hists_1 = []
    hists_2 = []

    fInCounter = 0
    fInNames = []
    for vbb in vbbs:
        for file in os.listdir("/Users/lucamicheletti/GITHUB/NA60Plus_test_beam/test_beam_november_2022/data"):
            d = os.path.join("/Users/lucamicheletti/GITHUB/NA60Plus_test_beam/test_beam_november_2022/data", file)
            if vbb in d and ithr in d:
                if not "ClusterCut" in d:
                    print(fInCounter)
                    fInNames.append(d)
                    fInCounter = fInCounter + 1
                    print(d)
                    break

            if fInCounter > plotLimit - 1:
                continue
                
    fInCounter = 0
    for fInName in fInNames:
        fIn = TFile.Open(fInName)
        for iPlane in range(0, 7):
            if fInCounter == 0:
                hists_0.append(fIn.Get("EventLoaderEUDAQ2/ALPIDE_{}/hPixelMultiplicityPerCorryEvent".format(iPlane)))
                hists_0[iPlane].SetDirectory(0)
                SetHistStyle2(hists_0[iPlane], "ALPIDE_{}".format(iPlane), "n. pixels", "entries", ROOT.kRed+1, 2, 3005, 0.5)
                hists_0[iPlane].Scale(1. / hists_0[iPlane].Integral())
            if fInCounter == 1:
                hists_1.append(fIn.Get("EventLoaderEUDAQ2/ALPIDE_{}/hPixelMultiplicityPerCorryEvent".format(iPlane)))
                hists_1[iPlane].SetDirectory(0)
                SetHistStyle2(hists_1[iPlane], "ALPIDE_{}".format(iPlane), "n. pixels", "entries", ROOT.kAzure+2, 2, 3005, 0.5)
                hists_1[iPlane].Scale(1. / hists_1[iPlane].Integral())
            if fInCounter == 2:
                hists_2.append(fIn.Get("EventLoaderEUDAQ2/ALPIDE_{}/hPixelMultiplicityPerCorryEvent".format(iPlane)))
                hists_2[iPlane].SetDirectory(0)
                SetHistStyle2(hists_2[iPlane], "ALPIDE_{}".format(iPlane), "n. pixels", "entries", ROOT.kGreen+1, 2, 3005, 0.5)
                hists_2[iPlane].Scale(1. / hists_2[iPlane].Integral())
        fInCounter = fInCounter + 1

    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.050)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)

    legend = ROOT.TLegend(0.30, 0.40, 0.60, 0.73, " ", "brNDC")
    SetLegend(legend)
    legend.AddEntry(hists_0[0], "Vbb = %s , iThr = %s" % (vbbs[0], ithr.replace("Ithr_", ""),), "EP")
    legend.AddEntry(hists_1[0], "Vbb = %s , iThr = %s" % (vbbs[1], ithr.replace("Ithr_", ""),), "EP")
    legend.AddEntry(hists_2[0], "Vbb = %s , iThr = %s" % (vbbs[2], ithr.replace("Ithr_", ""),), "EP")

    canvasPixelMult = ROOT.TCanvas("canvasPixelMult", "", 1200, 2400)
    canvasPixelMult.Divide(4, 2)
    for iPlane in range(0, 7):
        canvasPixelMult.cd(iPlane+1)
        gPad.SetLogx()
        gPad.SetLogy()
        hists_0[iPlane].Draw("EP")
        hists_1[iPlane].Draw("EPsame")
        hists_2[iPlane].Draw("EPsame")
        latexTitle.DrawLatex(0.2, 0.85, "#color[633]{Mean = %3.2f}" % (hists_0[iPlane].GetMean()))
        latexTitle.DrawLatex(0.2, 0.80, "#color[862]{Mean = %3.2f}" % (hists_1[iPlane].GetMean()))
        latexTitle.DrawLatex(0.2, 0.75, "#color[417]{Mean = %3.2f}" % (hists_2[iPlane].GetMean()))
    canvasPixelMult.cd(8)
    legend.Draw("same")
    canvasPixelMult.Update()

    canvasPixelMult.SaveAs("canvasPixelMult_scan.pdf")
    input()

    exit()
    

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--read_data", help="Read data and produce output histograms", action="store_true")
    args = parser.parse_args()

    if args.read_data:
        read_data("../test_beam_november_2022/data/resultscheck_0V_Ithr_240_473181837221123181842.root")

main()