from ROOT import TH1D,TH2D, TFile, TF2, TCanvas, kRed, TLegend, TFile, gStyle, kRainBow, TTree
import ROOT
import uproot
import numpy as np
import math
import argparse
import vertexer
import pandas as pd
from array import array
import yaml
import os
import sys
import time
import copy
from loading_bar import loading_bar
"""
TODO: 
- add loading bar to tell the progress
- add possibility to save the primary vertexes
- tell the execution time
- add event selections:
    - remove beam tails: x_PV, y_PV selections, may go outside the masked region
    - remove event with z_PV outside the target region 
    - only events where the vertex is reconstructed ncontributors > 0
    - beam cluster size > threshold; maybe cl1 > 20 or 30 
    - 

- add tracks selection:
    - cluster size? average cluster size?

- improve tree reading:
    - it crashes if the number of entries in the dataframe is too large

- plot results doesn't really work
"""

def read_tree(path, compute_sigmas = False, config = "", nsigma = 1, suffix =""):
    #put here the path of the file you want to analyze
    df_beam = uproot.open(path)['Tracking4D/tree/event/beam'].arrays(library="pd")
    df_tracks = uproot.open(path)['Tracking4D/tree/event/tracks'].arrays(library="pd", entry_stop=100)
    
    with open(os.path.expandvars(config), 'r') as stream:
        try:
            params = yaml.full_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    
    target = params["TARGET"]
    print(target)
    z_target = params["Z_TARGET"]
    thickness = params["THICKNESS"]
    interaction_prob = params["INTERACTION_PROB"]

    sigma_x = params["SIGMA_X_Pb"]
    sigma_y = params["SIGMA_Y_Pb"]
    mean_x = params["MEAN_X_Pb"]
    mean_y = params["MEAN_Y_Pb"]

    sigma_x_vtx = params["SIGMA_X_VTX"]
    sigma_y_vtx = params["SIGMA_Y_VTX"]
    mean_x_vtx = params["MEAN_X_VTX"]
    mean_y_vtx = params["MEAN_Y_VTX"]
    
    MIN_CL_BEAM = params["MIN_CL_BEAM"]

    old_index = 0

    vx = df_beam.at[df_beam.index[0],"beam.vx"]
    vy = df_beam.at[df_beam.index[0],"beam.vy"]
    px = df_beam.at[df_beam.index[0],"beam.px"]
    py = df_beam.at[df_beam.index[0],"beam.py"]


    x0 = vx*z_target+px
    y0 = vy*z_target+py

    ntracks = 0
    ntracks_sel = 0


    eta_list = []
    eta_list_sel = []
    #Vertex QA plots
    res_tracks_pb = TH2D("res_tracks_pb","track position - Pb position at z = z_target; x_{trk}-x_{Pb} (mm); y_{trk}-y_{Pb} (mm)",2000,-5,5,2000,-5,5)
    res_tracks_pv = TH2D("res_tracks_pv","track position - PV position at z = z_target; x_{trk}-x_{V} (mm); y_{trk}-y_{V} (mm)",2000,-5,5,2000,-5,5)
    hMultiplicity = TH1D("hMultiplicity",";track multiplicity;",500,0,499.5)
    hMultiplicitySel = TH1D("hMultiplicitySel",";track multiplicity;",500,0,499.5)
    hEta = TH1D("hEta",";#eta;",100,0,10)
    hEtaSel = TH1D("hEtaSel",";#eta;",100,0,10)
    hNcontr = TH1D("hNcontr",";number of contributors;",100,0,99.5)
    hNcontrVsZ = TH2D("hNcontrVsZ",";number of contributors;z_{v} (mm]);",100,0,99.5,1000,-2000,1000)
    hDeltaPbPVX = TH1D("hDeltaPbPVX",";x_{Pb}-x_{v} (mm);",1000,-0.5,0.5)
    hDeltaPbPVY = TH1D("hDeltaPbPVY",";y_{Pb}-y_{v} (mm);",1000,-0.5,0.5)
    hBeamVx = TH1D("hBeamVx","Beam direction x;beam v_{x};",1000,-0.1,0.1)
    hBeamVy = TH1D("hBeamVy","Beam direction y;beam v_{y};",1000,-0.1,0.1)
    hVertexX = TH1D("hVertexX",";x_{v} (mm);",100,-10,10)
    hVertexY = TH1D("hVertexY",";y_{v} (mm);",100,-10,10)
    hVertexZ = TH1D("hVertexZ",";z_{v} (mm);",100000,-2000,1000)
    hVertexZHighMult = TH1D("hVertexZHighMult",";z_{v} (mm);",100000,-2000,1000)
    file_ev = TFile(path, "r")
    hist_ev = file_ev.Get("EventLoaderEUDAQ2/ALPIDE_2/hPixelMultiplicityPerCorryEvent")
    nev_tot = hist_ev.GetEntries()
    #Vertexing
    df_vtx = pd.DataFrame(columns=df_tracks.columns.tolist())
    primary_vertexes = []
    # position (xv,yv,zv), nummber of contributors
    n_contributors = 0
    #loop over the index
    skip_event = False
    for index, row in df_tracks.iterrows():
        if skip_event:
            continue
        if old_index != index[0]:
            vertex = vertexer.fit_primary_vertex(df_vtx, df_beam.iloc[old_index],z_target)
            primary_vertexes.append([vertex[0],vertex[1],vertex[2],n_contributors])
            df_vtx.drop(df_vtx.index, inplace=True)
            hDeltaPbPVX.Fill(x0-vertex[0])
            hDeltaPbPVY.Fill(y0-vertex[1])
            old_index = index[0]
            #loading_bar(old_index,tot_event=nev_tot)
            ##iterations = old_index/nev_tot*50.
            ##loading_bar(##iterations)
            vx = df_beam.at[df_beam.index[old_index],"beam.vx"]
            vy = df_beam.at[df_beam.index[old_index],"beam.vy"]
            px = df_beam.at[df_beam.index[old_index],"beam.px"]
            py = df_beam.at[df_beam.index[old_index],"beam.py"]
            
            min_clsiz = min(df_beam.at[df_beam.index[old_index],"beam.size1"], df_beam.at[df_beam.index[old_index],"beam.size2"])

            hBeamVx.Fill(vx)
            hBeamVy.Fill(vy)
            hVertexX.Fill(vertex[0])
            hVertexY.Fill(vertex[1])
            if n_contributors > 0:
                hVertexZ.Fill(vertex[2])
            if n_contributors > 20:
                hVertexZHighMult.Fill(vertex[2])
            hNcontr.Fill(n_contributors)
            hNcontrVsZ.Fill(n_contributors, vertex[2])

            n_contributors = 0
            x0 = vx*z_target+px
            y0 = vy*z_target+py
            skip_event = False
            #print(min_clsiz)
            #if min_clsiz < MIN_CL_BEAM:
            #    skip_event = True
        
        x0 = vx*z_target+px
        y0 = vy*z_target+py
        x = row["tracks.px"]+z_target*row["tracks.vx"]
        y = row["tracks.py"]+z_target*row["tracks.vy"]
        dx = ((x-x0)-mean_x)/sigma_x
        dy = ((y-y0)-mean_y)/sigma_y
        res_tracks_pb.Fill(x-x0,y-y0)

        if dx**2+dy**2 < nsigma**2:#selecting the tracks 5 sigma from the primary vertex
            n_contributors += 1
            #row_to_add = df_tracks.iloc[0]  # Select the first row from df1
            df_vtx = df_vtx.append(row, ignore_index=True)

    #vertexer.define_PV_selection(res_tracks_pb, config, "Pb")
    
    output = TFile("../results/output_na60plus_nsigma"+str(nsigma)+"_"+suffix+".root","recreate")
    subdir_vtx = output.mkdir("vertex-qa")
    subdir_vtx.cd()
    res_tracks_pb.Write()
    hDeltaPbPVX.Write()
    hDeltaPbPVY.Write()
    hBeamVx.Write()
    hBeamVy.Write()
    hVertexX.Write()
    hVertexY.Write()
    hVertexZ.Write()
    hVertexZHighMult.Write()
    hNcontr.Write()
    hNcontrVsZ.Write()
    

    subdir_mult = output.mkdir("multiplicity")
    subdir_mult.cd()
    #analysis
    old_index = 0
    vertex = primary_vertexes[0]
    x0 = vertex[0]
    y0 = vertex[1]
    z0 = vertex[2]
    n_con = vertex[3]
    counter = 0
    #iterations = old_index/nev_tot*50
    #loading_bar(#iterations)
    
    for index, row in df_tracks.iterrows(): 
        if skip_event:
            continue

        if old_index != index[0]:
            old_index = index[0]   
            #loading_bar(old_index,tot_event=nev_tot,process="Analysis")
            
            #iterations = old_index/nev_tot*50
            #loading_bar(iterations)
            #position of the Pb ion
            #print(vertex)

            min_clsiz = min(df_beam.at[df_beam.index[old_index],"beam.size1"], df_beam.at[df_beam.index[old_index],"beam.size2"])
            skip_event = False
            #if min_clsiz < MIN_CL_BEAM:
            #    skip_event = True
            #else:
            counter += 1
            if counter == len(primary_vertexes):
                break
            vertex = primary_vertexes[counter]
            x0 = vertex[0]
            y0 = vertex[1]
            z0 = vertex[2]


            n_con = vertex[3]
            hMultiplicity.Fill(ntracks)
            hMultiplicitySel.Fill(ntracks_sel)
            ntracks = 0
            ntracks_sel = 0
        ntracks += 1

        #position of the track at z = z_target
        x = row["tracks.px"]+z0*row["tracks.vx"]
        y = row["tracks.py"]+z0*row["tracks.vy"]

        clsize_list = np.array([row["tracks.size1"], row["tracks.size2"], row["tracks.size3"], row["tracks.size4"], row["tracks.size5"], row["tracks.size6"]])
        clsize_list[clsize_list != 0]
        mean_clsize = np.mean(clsize_list)

        dx = ((x-x0)-mean_x_vtx)/sigma_x_vtx
        dy = ((y-y0)-mean_y_vtx)/sigma_y_vtx


        eta = math.atanh(1./math.sqrt(row["tracks.vx"]**2+row["tracks.vy"]**2+1))

        if dx**2+dy**2 < 5**2: #selecting the tracks 5 sigma from the primary vertex
            ntracks_sel += 1
            hEtaSel.Fill(eta)
        hEta.Fill(eta)
            #if n_con > 0:
        res_tracks_pv.Fill(x-x0, y-y0)
    
    #vertexer.define_PV_selection(res_tracks_pv, config, "VTX")
    res_tracks_pv.Write()
    hMultiplicity.Write()
    hMultiplicitySel.Write()
    hEta.Write()
    hEtaSel.Write()
    output.Close()

    sys.stdout.write('\n')
    sys.stdout.flush()

def plot_results(path_list, data_list, config_list):
    legend_multiplicity = TLegend(0.7,0.7,0.9,0.9)
    color_list = [
                   ROOT.kRed,
                   ROOT.kBlue,
                   ROOT.kGreen,
                   ROOT.kOrange,
                   ROOT.kBlack,
                      
                ]
    cv_multiplicity = TCanvas("cv_multiplicity", "cv_multiplicity")

    # Create a list to store the cloned histograms
    cloned_histograms = []

    for path, data, config, color in zip(path_list, data_list, config_list, color_list):
        fIn = TFile(path)
        hMulti = fIn.Get("multiplicity/hMultiplicitySel")
        file_ev = TFile(data, "r")
        hist_ev = file_ev.Get("EventLoaderEUDAQ2/ALPIDE_2/hPixelMultiplicityPerCorryEvent")
        nev_tot = hist_ev.GetEntries()

        with open(os.path.expandvars(config), 'r') as stream:
            try:
                params = yaml.full_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        TARGET = params["TARGET"]
        MASS_NUMBER = params["MASS_NUMBER"]

        # Create a copy of the histogram
        hMultiClone = hMulti.Clone("hMultiClone_" + TARGET)
        hMultiClone.Scale(1.0 / (MASS_NUMBER * nev_tot))
        hMultiClone.SetLineColor(color)
        hMultiClone.GetYaxis().SetTitle("event/(#trigger x A)")
        hMultiClone.GetXaxis().SetRangeUser(0, 150)

        cv_multiplicity.cd()
        hMultiClone.Draw("same")
        cloned_histograms.append(hMultiClone)  # Store the clone in the list
        legend_multiplicity.AddEntry(hMultiClone, TARGET, "lep")

    legend_multiplicity.Draw()
    cv_multiplicity.SetLogy()
    cv_multiplicity.SaveAs("multiplicity.png")





def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--read_tree", help="Read data and produce output histograms", action="store_true")
    parser.add_argument("--update", help="Update sigmas for PV selection", action="store_true")
    parser.add_argument("--plot_results", help="Produce plot from the analysis", action="store_true")
    parser.add_argument("config", help="Path to the YAML configuration file")
    args = parser.parse_args()


    if args.read_tree:
        suffix = "422231250231017232859_Be_planes_2_Unique_False"
        path = "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_"+suffix+".root"

        read_tree(path= path,compute_sigmas=True,config=args.config,nsigma=3, suffix=suffix)

    if args.plot_results:
        path_Ag = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422213212231017213217_Ag_planes_6.root"
        path_Air = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422191053231017193343_Pb_planes_6.root"
        path_Be = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422231249231017231254_Be_planes_6.root"
        path_Pb = "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_423133844231018133850_False_planes_6.root"
        path_S = "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_423140916231018140921_S_planes_6.root"
        path_list = [
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/results/output_na60plus_nsigma3_Sulfur.root",
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/results/output_na60plus_nsigma3_ARg.root"
                    ]
            
        config_list = [
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/macros/configs/linear_setup_S.yaml",
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/macros/configs/linear_setup_Ag.yaml"
                    ]
        
        data_list = [
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_423140916231018140921_S_planes_6.root",
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_422213213231017215647_Ag_planes_6.root"
                    ]
            
        plot_results(path_list, data_list, config_list)

main()
