import random
from scipy import optimize
from ROOT import TF2

def fit_primary_vertex(tracks, beam, z_target):

    x0 = beam["beam.vx"]*z_target+beam["beam.px"]
    y0 = beam["beam.vy"]*z_target+beam["beam.py"]

    #define the function
    def distToPvSquared(parameters):
        xv, yv, zv = parameters
        dist2 = 0
        for _, row in tracks.iterrows():
            dist2 += (row["tracks.vx"]*zv+row["tracks.px"]-xv)**2+(row["tracks.vy"]*zv+row["tracks.py"]-yv)**2
        dist2 += (beam["beam.vx"]*zv+beam["beam.px"]-xv)**2+(beam["beam.vy"]*zv+beam["beam.py"]-yv)**2
        return dist2

    # Initial guess
    initial_guess = [x0, y0 , z_target]
    #print("fit:")
    #print(tracks)
    #print(initial_guess)

    # Optimize the function using the BFGS method
    result = optimize.minimize(distToPvSquared, initial_guess, method='BFGS', tol=1e-5)

    #print(result.x)

    # Retrieve the results
    return result.x

def define_PV_selection(residuals, conf_file, type = "Pb"):

    tf2 = TF2("xygaus","xygaus",-0.20,0.20,-0.2,0.2)
    residuals.Fit(tf2,"MR0")

    sigma_x = tf2.GetParameter(2)
    sigma_y = tf2.GetParameter(4)
    mean_x = tf2.GetParameter(1)
    mean_y = tf2.GetParameter(3)
    write_selection(conf_file, [sigma_x,sigma_y,mean_x,mean_y], type)


def write_selection(conf_file, parameters = [0,0,0,0], type = "Pb"):
    counter = 0
    with open(conf_file, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break

            if f"SIGMA_X_{type}" in line.strip():
                counter_SIGMA_X = counter
            if f"SIGMA_Y_{type}" in line.strip():
                counter_SIGMA_Y = counter
            if f"MEAN_X_{type}" in line.strip():
                counter_MEAN_X = counter
            if f"MEAN_Y_{type}" in line.strip():
                counter_MEAN_Y = counter
            counter += 1
    
    # Open the file in read mode and read all its lines into a list
    with open(conf_file, 'r') as file:
        lines = file.readlines()

    lines[counter_SIGMA_X] = f"SIGMA_X_{type} : {parameters[0]} \n"
    lines[counter_SIGMA_Y] = f"SIGMA_Y_{type} : {parameters[1]} \n"
    lines[counter_MEAN_X] = f"MEAN_X_{type} : {parameters[2]} \n"
    lines[counter_MEAN_Y] = f"MEAN_Y_{type} : {parameters[3]} \n"
    # Open the file in write mode and write the modified list to the file
    with open(conf_file, 'w') as file:
        file.writelines(lines)