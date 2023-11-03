from ROOT import TH1D,TH2D, TFile, TF2, TCanvas, kRed, TLegend, TFile, gStyle, kRainBow, TTree
import uproot
import numpy as np
import math
import argparse
import vertexer
import pandas as pd
from array import array
import yaml
import os

def read_tree(path, compute_sigmas = False, config = "", nsigma = 1, suffix =""):
    #put here the path of the file you want to analyze
    df_beam = uproot.open(path)['Tracking4D/tree/event/beam'].arrays(library="pd")
    df_tracks = uproot.open(path)['Tracking4D/tree/event/tracks'].arrays(library="pd")

    with open(os.path.expandvars(config), 'r') as stream:
        try:
            params = yaml.full_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    
    target = params["TARGET"]
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
    hNcontr = TH1D("hNcontr",";number of contributors;",100,0,99.5)
    hNcontrVsZ = TH2D("hNcontrVsZ",";number of contributors;z_{v} (mm]);",100,0,99.5,100,0,150)
    hDeltaPbPVX = TH1D("hDeltaPbPVX",";x_{Pb}-x_{v} (mm);",1000,-0.5,0.5)
    hDeltaPbPVY = TH1D("hDeltaPbPVY",";y_{Pb}-y_{v} (mm);",1000,-0.5,0.5)
    hVertexX = TH1D("hVertexX",";x_{v} (mm);",100,-10,10)
    hVertexY = TH1D("hVertexY",";y_{v} (mm);",100,-10,10)
    hVertexZ = TH1D("hVertexZ",";z_{v} (mm);",1000,25,125)
    hVertexZHighMult = TH1D("hVertexZHighMult",";z_{v} (mm);",1000,25,125)

    #Vertexing
    df_vtx = pd.DataFrame(columns=df_tracks.columns.tolist())
    primary_vertexes = []
    # position (xv,yv,zv), nummber of contributors
    n_contributors = 0
    for index, row in df_tracks.iterrows():
        if old_index != index[0]:
            vertex = vertexer.fit_primary_vertex(df_vtx, df_beam.iloc[old_index],z_target)
            primary_vertexes.append([vertex[0],vertex[1],vertex[2],n_contributors])
            df_vtx.drop(df_vtx.index, inplace=True)
            hDeltaPbPVX.Fill(x0-vertex[0])
            hDeltaPbPVY.Fill(y0-vertex[1])
            old_index = index[0]
            vx = df_beam.at[df_beam.index[old_index],"beam.vx"]
            vy = df_beam.at[df_beam.index[old_index],"beam.vy"]
            px = df_beam.at[df_beam.index[old_index],"beam.px"]
            py = df_beam.at[df_beam.index[old_index],"beam.py"]
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

    vertexer.define_PV_selection(res_tracks_pb, config, "Pb")
    output = TFile("../results/output_na60plus_nsigma"+str(nsigma)+"_"+suffix+".root","recreate")
    subdir_vtx = output.mkdir("vertex-qa")
    subdir_vtx.cd()
    res_tracks_pb.Write()
    hDeltaPbPVX.Write()
    hDeltaPbPVY.Write()
    hVertexX.Write()
    hVertexY.Write()
    hVertexZ.Write()
    hVertexZHighMult.Write()
    hNcontr.Write()
    hNcontrVsZ.Write()

    subdir_mult = output.mkdir("multiplicity")
    subdir_mult.cd()
    #analysis
    print(len(primary_vertexes))
    old_index = 0
    vertex = primary_vertexes[0]
    x0 = vertex[0]
    y0 = vertex[1]
    z0 = vertex[2]
    n_con = vertex[3]
    counter = 0
    for index, row in df_tracks.iterrows():
        if old_index != index[0]:
            old_index = index[0]
            counter += 1
            if counter == len(primary_vertexes):
                break
            vertex = primary_vertexes[counter]
            #position of the Pb ion
            x0 = vertex[0]
            y0 = vertex[1]
            z0 = vertex[2]

            n_con = vertex[3]
            #print(vertex)

            ntracks = 0
            ntracks_sel = 0
        ntracks += 1

        #position of the track at z = z_target
        x = row["tracks.px"]+z0*row["tracks.vx"]
        y = row["tracks.py"]+z0*row["tracks.vy"]

        clsize_list = np.array([row["tracks.size1"], row["tracks.size2"], row["tracks.size3"], row["tracks.size4"], row["tracks.size5"], row["tracks.size6"]])
        clsize_list[clsize_list != 0]
        mean_clsize = np.mean(clsize_list)

        dx = ((x-x0)-mean_x)/sigma_x_vtx
        dy = ((y-y0)-mean_y)/sigma_y_vtx


        eta = math.atanh(1./math.sqrt(row["tracks.vx"]**2+row["tracks.vy"]**2+1))

        #if dx**2+dy**2 < 5**2: #selecting the tracks 5 sigma from the primary vertex
        ntracks_sel += 1
        #if n_con > 0:
        res_tracks_pv.Fill(x-x0, y-y0)
    
    vertexer.define_PV_selection(res_tracks_pv, config, "VTX")
    res_tracks_pv.Write()
    output.Close()

def plot_results():
    fIn = TFile("../results/output_na60plus.root")
    th1vxy = fIn.Get("tracks_vxy")
    th1pxy = fIn.Get("tracks_pxy")
    proj_eta = fIn.Get("central_eta")
    th1res_pxy = fIn.Get("res_tracks_pxy")
    th1res_pxy_sel = fIn.Get("res_tracks_pxy_sel")
    th2_eta_vs_mult = fIn.Get("eta_vs_mult")
    th2_eta_vs_mult_sel = fIn.Get("eta_vs_mult_sel")
    th1_multi = fIn.Get("multiplicity")
    th1_multi_sel = fIn.Get("multiplicity selected")
    th1_eta = fIn.Get("eta")
    th1_eta_sel = fIn.Get("eta_sel")
    th1_mean_clsize = fIn.Get("th1_mean_clsize")
    th1_mean_clsize_sel = fIn.Get("th1_mean_clsize_sel")
    th2_eta_vs_mean_clsize = fIn.Get("eta_vs_mean_clsize")
    th2_eta_vs_mean_clsize_sel = fIn.Get("eta_vs_mean_clsize_sel")

    legend_eta = TLegend(0.7,0.7,0.9,0.9)
    legend_eta.AddEntry(th1_eta,"Before selection","lep")
    legend_eta.AddEntry(th1_eta_sel,"After selection","lep")

    legend_multi = TLegend(0.7,0.7,0.9,0.9)
    legend_multi.AddEntry(th1_multi,"Before selection","lep")
    legend_multi.AddEntry(th1_multi_sel,"After selection","lep")

    cv = TCanvas("cv","cv",1600,1200)
    th1pxy.Draw("colz")
    cv.SaveAs("../results/figure/tracksAtZtarget.png")

    proj_eta.Draw()
    cv.SaveAs("../results/figure/centralEtaSelected.png")

    th1res_pxy.Draw("colz")
    cv.SaveAs("../results/figure/resAtZtarget.png")

    th1_multi.SetLineColor(kRed)
    th1_multi.Draw()
    legend_multi.Draw()
    th1_multi_sel.Draw("same")
    cv.SetLogy()
    cv.SaveAs("../results/figure/multiplicity.png")

    th1_eta.SetLineColor(kRed)
    th1_eta.Draw()
    cv.SetLogy(0)
    legend_eta.Draw()
    th1_eta_sel.Draw("same")
    cv.SaveAs("../results/figure/eta.png")

    proj_eta.Draw()
    cv.SaveAs("../results/figure/centralEtaSelected.png")

    cv2 = TCanvas("cv2","cv2",800,600)
    gStyle.SetPalette(kRainBow)
    th2_eta_vs_mean_clsize.Draw("COLZ")
    cv2.SaveAs("../results/figure/etaVsMeanCLsize.pdf")

    th2_eta_vs_mean_clsize_sel.Draw("COLZ")
    cv2.SaveAs("../results/figure/etaVsMeanCLsizeSelected.pdf")


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--read_tree", help="Read data and produce output histograms", action="store_true")
    parser.add_argument("--update", help="Update sigmas for PV selection", action="store_true")
    parser.add_argument("--plot_results", help="Produce plot from the analysis", action="store_true")
    parser.add_argument("config", help="Path to the YAML configuration file")
    args = parser.parse_args()


    if args.read_tree:
        sigma_list = [3]
        for sigma in sigma_list:
            path_new = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422213212231017213217_Ag_planes_6.root"
            path_old = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422191053231017193343_Pb_planes_6.root"
            path_be = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422231249231017231254_Be_planes_6.root"
            read_tree(path= path_be,compute_sigmas=args.update,config=args.config,nsigma=sigma, suffix="422213212231017213217_Be_planes_6")

    if args.plot_results:
        plot_results()


main()