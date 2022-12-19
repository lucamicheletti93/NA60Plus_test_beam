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
from plot_library import LoadStyle, SetLegend, SetHistStyle, SetLatex, SetHistStyle2, SetHistStyle3, SetHistStyle4

###
def read_data(path):
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
        for file in os.listdir(path):
            d = os.path.join(path, file)
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
        #gPad.SetLogy()
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
    
def beam_spot(fInName):
    alpidePitch = 28 #um
    fIn = TFile.Open(fInName)
    histHitMap = []
    histX = []
    histY = []
    funcX = []
    funcY = []
    histScanX = TH1F("histScanX", "", 7, 0, 7)
    SetHistStyle4(histScanX, "histScanX", "Plane", "Width", ROOT.kAzure+2, 2)
    histScanY = TH1F("histScanY", "", 7, 0, 7)
    SetHistStyle4(histScanY, "histScanY", "Plane", "Width", ROOT.kAzure+2, 2)

    for iPlane in range(0, 7):
        histHitMap.append(fIn.Get("EventLoaderEUDAQ2/ALPIDE_{}/hitmap".format(iPlane)))
        histHitMap[iPlane].SetDirectory(0)
        histX.append(histHitMap[iPlane].ProjectionX("ALPIDE_{}_X".format(iPlane)))
        SetHistStyle3(histX[iPlane], "ALPIDE_{}_X".format(iPlane), "n. pixels", "entries", ROOT.kBlack, 2, 3005, 0.5)
        funcX.append(TF1("func_ALPIDE_{}_X".format(iPlane), "gaus", 0, 500))
        histY.append(histHitMap[iPlane].ProjectionY("ALPIDE_{}_Y".format(iPlane)))
        SetHistStyle3(histY[iPlane], "ALPIDE_{}_Y".format(iPlane), "n. pixels", "entries", ROOT.kBlack, 2, 3005, 0.5)
        funcY.append(TF1("func_ALPIDE_{}_Y".format(iPlane), "gaus", 0, 500))

    canvasX = ROOT.TCanvas("canvasX", "", 1200, 2400)
    canvasX.Divide(4, 2)
    for iPlane in range(0, 7):
        canvasX.cd(iPlane+1)
        gStyle.SetOptFit(1111)
        histX[iPlane].Draw("EP")
        histY[iPlane].Fit(funcX[iPlane])
        histScanX.SetBinContent(iPlane+1, funcX[iPlane].GetParameter(2) * alpidePitch)
        histScanX.SetBinError(iPlane+1, funcX[iPlane].GetParError(2) * alpidePitch)
    canvasX.cd(8)
    gStyle.SetOptStat(0)
    histScanX.Draw("EPsame")
    canvasX.Update()
    canvasX.SaveAs("x_projection.pdf")

    canvasY = ROOT.TCanvas("canvasY", "", 1200, 2400)
    canvasY.Divide(4, 2)
    for iPlane in range(0, 7):
        canvasY.cd(iPlane+1)
        gStyle.SetOptFit(1111)
        histY[iPlane].Draw("EP")
        histY[iPlane].Fit(funcY[iPlane])
        histScanY.SetBinContent(iPlane+1, funcY[iPlane].GetParameter(2) * alpidePitch)
        histScanY.SetBinError(iPlane+1, funcY[iPlane].GetParError(2) * alpidePitch)
    canvasY.cd(8)
    gStyle.SetOptStat(0)
    histScanY.Draw("EPsame")
    canvasY.Update()
    canvasY.SaveAs("y_projection.pdf")

    input()

    



def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--read", help="Read data and produce output histograms", action="store_true")
    parser.add_argument("--beam_spot", help="Extract the beam spot via fit", action="store_true")
    args = parser.parse_args()

    if args.read:
        read_data("/Users/lucamicheletti/GITHUB/NA60Plus_test_beam/data/test_beam_november_2022")
    if args.beam_spot:
        beam_spot("/Users/lucamicheletti/GITHUB/NA60Plus_test_beam/data/test_beam_november_2022/resultscheck_0V_Ithr_240_473212644221123212649.root")


main()