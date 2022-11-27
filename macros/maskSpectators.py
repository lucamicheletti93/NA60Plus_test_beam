import ROOT

#path of the file with the hitmap
path = "/home/giacomo/its-corryvreckan-tools/output/SPSNovember22/alignment_475230516221125231020.root"
#alpide that we whant to mask
alpide_list = [1,2,3,4,5]
z_target = 125
z_position = [150,175,200,225,250]
pixel_pitch = [29.24,26.88]
file = ROOT.TFile(path)
radius0 = 1000

for alpide,z_alpide in zip(alpide_list,z_position):
    hitmap = file.Get("ClusteringSpatial/ALPIDE_"+str(alpide)+"/clusterPositionLocal")

    mean_x = hitmap.GetMean(1)
    mean_y = hitmap.GetMean(2)
    
    theta_list = []

    radius = radius0*(z_alpide-z_target)/125.

    X = radius/pixel_pitch[0]
    Y = radius/pixel_pitch[1]

    x0 = int(mean_x)
    y0 = int(mean_y)

    f = open("masking_"+str(alpide)+".txt", "w")
    for ix in range(0,1024):
        for iy in range(0,512):
            r2 = (ix-x0)**2/X**2+(iy-y0)**2/Y**2
            if r2 < 1:
                f.write("p	"+str(ix)+"  "+str(iy)+"\n")
    f.close()