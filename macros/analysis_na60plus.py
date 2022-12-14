from ROOT import TH1D,TH2D, TFile, TF2, TCanvas, kRed, TLegend, TFile
import uproot
import math

#put here the path of the file you want to analyze
path = "/home/giacomo/its-corryvreckan-tools/output/SPSNovember22/alignment_475230516221125231020.root"

df_beam = uproot.open(path)['Tracking4D/tree/event/beam'].arrays(library="pd")
df_tracks = uproot.open(path)['Tracking4D/tree/event/tracks'].arrays(library="pd")

th1vxy = TH2D("tracks_vxy",";;",1000,-0.1,0.1,1000,-0.1,0.1)
th1pxy = TH2D("tracks_pxy","track position at z = z_target; x (mm); y(mm)",2000,-5,5,2000,-5,5)
th1res_pxy = TH2D("res_tracks_pxy","track position - Pb position at z = z_target; x_{trk}-x_{Pb} (mm); y_{trk}-y_{Pb} (mm)",2000,-5,5,2000,-5,5)
th2_eta_vs_mult = TH2D("eta_vs_mult",";#eta;multiplicity;counts",100,0,10,400,-0.5,399.5)
th2_eta_vs_mult_sel = TH2D("eta_vs_mult_sel",";#eta;multiplicity;counts",100,0,10,400,-0.5,399.5)
th1_multi = TH1D("multiplicity","multiplicity;number of tracks per event;counts",400,-0.5,399.5)
th1_multi_sel = TH1D("multiplicity selected","multiplicity;number of tracks per event;counts",400,-0.5,399.5)
th1_eta = TH1D("eta","eta;#eta;counts",100,0,10)
th1_eta_sel = TH1D("eta_sel","eta;#eta;counts",100,0,10)

#position of the target in the telescope
z_target = 125 #mm

old_index = 0

x0 = df_beam.at[df_beam.index[0],"beam.x"]
y0 = df_beam.at[df_beam.index[0],"beam.y"]

ntracks = 0
ntracks_sel = 0

sigma_x = 0.05292421774392007
sigma_y = 0.040370864336865826
mean_x = -0.015625873902035948
mean_y = -0.024083787738233064
eta_list = []
eta_list_sel = []
for index, row in df_tracks.iterrows():
    if old_index != index[0]:
        old_index = index[0]
        #position of the Pb ion
        x0 = df_beam.at[df_beam.index[old_index],"beam.x"]
        y0 = df_beam.at[df_beam.index[old_index],"beam.y"]
        th1_multi.Fill(ntracks)
        th1_multi_sel.Fill(ntracks_sel)
        for eta in eta_list:
            th2_eta_vs_mult.Fill(eta, ntracks)
        for eta in eta_list_sel:
            th2_eta_vs_mult_sel.Fill(eta, ntracks_sel)
        eta_list = []
        eta_list_sel = []
        ntracks = 0
        ntracks_sel = 0
    ntracks += 1
    #position of the track at z = z_target
    x = row["tracks.px"]+z_target*row["tracks.vx"]
    y = row["tracks.py"]+z_target*row["tracks.vy"]
    
    th1vxy.Fill(row["tracks.vx"],row["tracks.vy"])
    th1pxy.Fill(x,y)
    th1res_pxy.Fill(x-x0,y-y0)

    dx = ((x-x0)-mean_x)/sigma_x
    dy = ((y-y0)-mean_y)/sigma_y

    p =  math.sqrt(row["tracks.vx"]**2+row["tracks.vy"]**2+1)
    eta = math.log((p+1)/(p-1))/2.

    if dx**2+dy**2 < 5**2:#selecting the tracks 5 sigma from the primary vertex
        ntracks_sel += 1
        th1_eta_sel.Fill(eta)
        eta_list_sel.append(eta)

    th1_eta.Fill(eta)
    eta_list.append(eta)

tf2 = TF2("xygaus","xygaus",-0.25,0.25,-0.2,0.2)
th1res_pxy.Fit(tf2,"MR0")

sigma_x = tf2.GetParameter(2)
sigma_y = tf2.GetParameter(4)
mean_x = tf2.GetParameter(1)
mean_y = tf2.GetParameter(3)

#values used before to reject the tracks far from the primary vertex
print("sigma_x: ",sigma_x)
print("sigma_y: ",sigma_y)
print("mean_x: ",mean_x)
print("mean_y: ",mean_y)


#first estimation of the eta in the 0-5% centrality class
#number of events in the run
corry_output = TFile(path,"read") 
n_events = corry_output.Get("EventLoaderEUDAQ2/ALPIDE_1/hPixelMultiplicityPerCorryEvent").GetEntries()
#interaction probability of the beam with the target
pint = 0.15
#centrality class
centrality = [0.0,0.05]
#expected number of events
n_ev = (centrality[1]-centrality[0])*pint*n_events 
#select the bin after which we have the x most central event 
integral = 0
bin_min = 0
for i in range(th1_multi_sel.GetNbinsX(), 0, -1):
    integral += th1_multi_sel.GetBinContent(i)
    if integral > n_ev:
        bin_min = i
        break

nbins = th2_eta_vs_mult_sel.GetNbinsY()
proj_eta = th2_eta_vs_mult_sel.ProjectionX("central_eta",bin_min, nbins)

legend_eta = TLegend(0.7,0.7,0.9,0.9)
legend_eta.AddEntry(th1_eta,"Before selection","lep")
legend_eta.AddEntry(th1_eta_sel,"After selection","lep")

legend_multi = TLegend(0.7,0.7,0.9,0.9)
legend_multi.AddEntry(th1_multi,"Before selection","lep")
legend_multi.AddEntry(th1_multi_sel,"After selection","lep")

#save the results
output = TFile("../results/output_na60plus.root","recreate")
cv = TCanvas("cv","cv",1600,1200)
th1vxy.Write()
th1pxy.Write()
th1pxy.Draw("colz")
cv.SaveAs("../results/figure/tracksAtZtarget.png")
proj_eta.Write()

proj_eta.Draw()
cv.SaveAs("../results/figure/centralEtaSelected.png")

th1res_pxy.Write()
th1res_pxy.Draw("colz")
cv.SaveAs("../results/figure/resAtZtarget.png")

th1_multi.Write()
th1_multi_sel.Write()

th1_multi.SetLineColor(kRed)
th1_multi.Draw()
legend_multi.Draw()
th1_multi_sel.Draw("same")
cv.SetLogy()
cv.SaveAs("../results/figure/multiplicity.png")

th1_eta.Write()
th1_eta_sel.Write()

th1_eta.SetLineColor(kRed)
th1_eta.Draw()
cv.SetLogy(0)
legend_eta.Draw()
th1_eta_sel.Draw("same")
cv.SaveAs("../results/figure/eta.png")

th2_eta_vs_mult_sel.Write()
th2_eta_vs_mult.Write()
proj_eta.Draw()
cv.SaveAs("../results/figure/centralEtaSelected.png")
output.Close()