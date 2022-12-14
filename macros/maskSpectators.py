import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--merge', help='Merge with other masking', action='store_true')
parser.add_argument('-c', '--center', help='Mask the center of the ALPIDEs', action='store_true')
parser.add_argument('-e', '--edges', help='Mask the edge of the ALPIDEs', action='store_true')
args = parser.parse_args()


MERGE = args.merge
MASK_CENTER = args.center
MASK_EDGES = args.edges

#path of the file with the hitmap
path = "/home/giacomo/its-corryvreckan-tools/output/SPSNovember22/for_masking_475230516221125231020.root"
#alpide that we whant to mask
alpide_list = [1,2,3,4,5]
z_target = 125
target_thickness = 5
z_position = [150,175,200,225,250]
pixel_pitch = [29.24,26.88]
file = ROOT.TFile(path)
radius0 = 1000
target_distance_last = z_position[-1] -z_target -target_thickness/2.
npixels_x = 1024
npixels_y = 512
file_name_list = [
                    "masking_1.txt",
                    "masking_2.txt",
                    "masking_3.txt",
                    "masking_4.txt",
                    "masking_5.txt",
                ]

test_file = ROOT.TFile("test.root","recreate")
for alpide,z_alpide,file_name in zip(alpide_list,z_position,file_name_list):
    #if MERGE:
    #    f = open(file_name, "a")
    #else:
    if MASK_CENTER:
        th2 = ROOT.TH2D("th2_"+str(alpide),";x;y",1024,-0.5,1023.5,512,-0.5,511.5)
        f = open("masking_"+str(alpide)+".txt", "w")
        print("masking_"+str(alpide)+".txt")
        hitmap = file.Get("ClusteringSpatial/ALPIDE_"+str(alpide)+"/clusterPositionLocal")

        mean_x = hitmap.GetMean(1)
        mean_y = hitmap.GetMean(2)
        print("ALPIDE:",alpide)
        print("mean x: ",mean_x)
        print("mean y: ",mean_y)
        theta_list = []

        radius = radius0*(z_alpide-z_target)/125.
        X = radius/pixel_pitch[0]
        Y = radius/pixel_pitch[1]

        for ix in range(0,npixels_x):
            for iy in range(0,npixels_y):
                r2 = (ix-mean_x)**2/X**2+(iy-mean_y)**2/Y**2
                if r2 < 1:
                    f.write("p	"+str(ix)+"  "+str(iy)+"\n")
                    th2.SetBinContent(ix,iy,1)
        f.close()
    #work in progress    
    if False:
        if alpide == 5:
            d1 = npixels_y - mean_y
            if d1 > mean_y:
                d1 = mean_y
            
            X = d1*pixel_pitch[1]/pixel_pitch[0]
            Y = d1

            for ix in range(0,npixels_x):
                for iy in range(0,npixels_y):
                    r2 = (ix-mean_x)**2/X**2+(iy-mean_y)**2/Y**2
                    if r2 > 1:
                        f.write("p	"+str(ix)+"  "+str(iy)+"\n")
                        th2.SetBinContent(ix,iy,1)
            th2.Write()

    if MASK_EDGES and alpide != alpide_list[-1]:
        if False:
            target_distance = z_alpide -z_target -target_thickness/2.
            Lxi = int(target_distance/target_distance_last * npixels_x / 2.)
            Lyi = int(target_distance/target_distance_last * npixels_y / 2.)
            print("ALPIDE:",alpide)
            print(Lxi)
            print(Lyi)
            for iy in range(0,npixels_y):
                for ix in range(0, Lxi+1):
                    f.write("p	"+str(ix)+"  "+str(iy)+"\n")
                for ix in range(int(npixels_x/2.)+Lxi,npixels_x):
                    f.write("p	"+str(ix)+"  "+str(iy)+"\n")

            for ix in range(Lxi+1, int(npixels_x/2.)+Lxi+1):
                for iy in range(0,Lyi+1):
                    f.write("p	"+str(ix)+"  "+str(iy)+"\n")
                for iy in range(int(npixels_y/2.)+Lyi, npixels_y):
                    f.write("p	"+str(ix)+"  "+str(iy)+"\n")
test_file.Close()
