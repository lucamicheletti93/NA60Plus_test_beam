import ROOT

#path of the file with the hitmap
path = "alignment_475230516221125231020.root"
#alpide that we whant to mask
ALPIDE = "5"

file = ROOT.TFile(path)
hitmap = file.Get("ClusteringSpatial/ALPIDE_"+ALPIDE+"/clusterPositionLocal")

mean_x = hitmap.GetMean(1)
mean_y = hitmap.GetMean(2)

pixel_pitch = [29.24,26.88]
radius = 1000

theta_list = []

X = radius/pixel_pitch[0]
Y = radius/pixel_pitch[1]

x0 = int(mean_x)
y0 = int(mean_y)

f = open("mask_ALPIDE_"+ALPIDE+"_2.txt", "w")
for ix in range(0,1024):
    for iy in range(0,512):
        r2 = (ix-x0)**2/X**2+(iy-y0)**2/Y**2
        if r2 < 1:
            f.write("p	"+str(ix)+"  "+str(iy)+"\n")
f.close()