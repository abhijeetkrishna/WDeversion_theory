from methods import *

#print("entered array_wd.py")

#Taking arguments from command line
map_index_dest=sys.argv[1]
task_id=int(sys.argv[2])
#print("map_index_dest : " + str(map_index_dest))
#print("task id : " + str(task_id))

#getting the dictionary of variables and their default values
var_dict = {
    'seed':0,
    'nb_iterations':20,
    'thickness':0.1,
    'R_initial':1,
    'theta_DV':0.1931,
    'DV_present':True,
    #'outDV_gradient':True,
    #'volume_conservation':False,
    'k_type':'k_c',
    'theta_max': 0.8662, #32*np.pi/180,
    'theta_ref': 0.8662,
    'overwrite_old_simulation':True,
    #'tol':1e-6,
    'dt':0.01,
    'angle_of_rotation':0.1,
    'mesh_refine_factor':30,
    'tol_by_dt':1e-8,
    'isotropic_contribution':"all", #no or all
    'anisotropic_contribution':"all", #all or no or elongation or rearrangement
    'height_contribution':"no", #no or all
    'genotype':"ecadGFPnbG4",#"fit_lambdas_df_ecadGFPnbG4.pkl", #fit_lambdas_df_ecadGFPnbG4myoVI.pkl
    'disc_name':'isotropic_homogeneous',
    #'inDV_lambda_anisotropic_coeffs' : [0.5, 1]
}

#reading the map_index csv
#print('reading the map_index csv')
map_index=pd.read_csv(map_index_dest)
input_var=list(map_index.columns) #variables for which values are being imported
for i in range(len(input_var)):
    if input_var[i] == 'folder_name':
        dirname = map_index.loc[task_id,str(input_var[i])]
        continue
    var_dict[input_var[i]]=map_index.loc[task_id,str(input_var[i])]

#initializing variables
#print('initializing variables')
for key,val in var_dict.items():
    exec(key + '=val')
#set the seed
np.random.seed(seed)
#set the tol 
tol = dt*tol_by_dt

###############
#Generating mesh #
###############

angle_of_rotation = np.round(np.random.uniform(low = 0, high = 2*np.pi),4)
#print("angle_of_rotation : " + str(angle_of_rotation))

#read spherical mesh
sphere_pickle = f"/projects/project-krishna/WingDiscEversion_theory/Simulations/meshes/icosaspherical_cap/pickle/IcoSph_mesh_refine_factor_{mesh_refine_factor}.pkl"
[sphere_balls_df,sphere_edges_df,sphere_cells_df] = pickle.load(open(sphere_pickle,"rb"))
#if this file is not there then we can create it using get_meshzoo_icosa_sphere function in methods.py 

# Crop and add thickness to mesh #
[balls_df, springs_df] = add_thickness_to_mesh(v_df = sphere_balls_df[['x', 'y', 'z']], edges_df = sphere_edges_df, #cells_df = cells_df,
                                                            rotate_mesh_bool=True, angle_of_rotation = angle_of_rotation,
                                                            crop_mesh_bool=True, theta_max=theta_max,thickness = thickness,
                                                            )

#add basis vectors to mesh
[balls_df, springs_df] = add_basis_vectors_to_Sph(balls_df, springs_df, theta_DV=theta_DV)
pickle.dump([balls_df, springs_df], open(dirname+'runfiles/init_mesh.pickle', 'wb'))

####################################
# Choosing the lambda coefficients #
####################################
#TO DO - height info is discontinued so no need to include lambda_h

props = [
"lambda_isotropic_coeffs", "lambda_anisotropic_coeffs",
'inDV_lambda_isotropic_coeffs', 'inDV_lambda_anisotropic_coeffs',
#'lambda_height_coeffs', 'inDV_lambda_height_coeffs'
]
def read_lambda_values_from_pkl(stages_df = None, stage = None,
    isotropic_contribution = "all", anisotropic_contribution = "all", #height_contribution = "no",
    ):
    #this function will read the lambda values from the stages_df
    #first we take filter the stage
    #first we make the array props
    #for each of the value first we check if we want to take its real value or keep it 1
    #To Do - user should be able to assign a value for each of the prop

    temp_stages_df = stages_df[stages_df["stage"] == stage].reset_index(drop = True)
    #print("reading " + temp_stages_df.iloc[0]["stage_name"])

    values = [temp_stages_df[temp_stages_df["prop"] == prop]["value"].values[0] for prop in props]

    #if isotropic contribution is None then we put values as [1]
    if isotropic_contribution == "no":
        values[0] = np.array([1]) #lambda_isotropic_coeffs
        values[2] = np.array([1]) #inDV_lambda_isotropic_coeffs

    #if anisotropic contribution is None then we put values as [1]
    if anisotropic_contribution == "no":
        values[1] = np.array([1]) #lambda_anisotropic_coeffs
        values[3] = np.array([1]) #inDV_lambda_anisotropic_coeffs
    elif anisotropic_contribution == "elongation":
        values[1] = temp_stages_df[temp_stages_df["prop"] == "lambda_Q_coeffs"]["value"].values[0]
        values[3] = temp_stages_df[temp_stages_df["prop"] == "inDV_lambda_Q_coeffs"]["value"].values[0]
    elif anisotropic_contribution == "rearrangement":
        values[1] = temp_stages_df[temp_stages_df["prop"] == "lambda_rearrangement_coeffs"]["value"].values[0]
        values[3] = temp_stages_df[temp_stages_df["prop"] == "inDV_lambda_rearrangement_coeffs"]["value"].values[0]

    #if volume conservation is false then we put lambda_h coefficients as 1
    #if volume_conservation is True and lambda_height_coeffs is None then 
    #we will put lambda_height_coeffs as 1/(lambda_isotropic^2)
    #To make sure that the lambda_height_coeffs read from the file is implemented
    #height_contribution has to be set to  "all" or "yes", just not "no"
    #if height_contribution == "no":
    #    values[4] = [1] #lambda_height_coeffs
    #    values[5] = [1] #inDV_lambda_height_coeffs

    return(values)

###############################
# Reading a lambdas dataframe #
###############################

#Read

import shutil
#shutil.copy('/projects/project-krishna/WingDiscEversion_theory/Experiments/analysis/fit_lambda_files/fit_lambdas_df_'+genotype+'.pkl', dirname+'runfiles/fit_lambdas_df.pkl')
#shutil.copy('/projects/project-krishna/WingDiscEversion_theory/Experiments/analysis/fit_lambda_files/fit_lambdas_df_'+genotype+'.csv', dirname+'runfiles/fit_lambdas_df.csv')
shutil.copy('input_lambda_df.pkl', dirname+'runfiles/fit_lambdas_df.pkl')
shutil.copy('input_lambda_df.csv', dirname+'runfiles/fit_lambdas_df.csv')

stages_df = pd.read_pickle(dirname+'runfiles/fit_lambdas_df.pkl') #pd.read_csv("stages_df.csv")
temp_stages_df = stages_df.copy(deep = True)

#stages_df = pd.read_pickle(dirname+'runfiles/fit_lambdas_df.pkl') #pd.read_csv("stages_df.csv")
#stages_df = stages_df.sort_values(by = 'stage')
stages = np.unique(stages_df["stage"])
#stage_names = np.unique(stages_df["stage_name"])

for j in range(len(stages)):
    stage = stages[j]
    #temp_stages_df = stages_df[stages_df["stage"] == stage].reset_index(drop = True)

    values = read_lambda_values_from_pkl(stages_df,stage,
        isotropic_contribution, anisotropic_contribution, #height_contribution,
        )
    for prop,value in zip(props,values):
        #print(prop + " : " + str(value))
        exec(prop+"=value")
        exec(prop.replace("_coeffs","_obj") + "=np.poly1d(value)")
        #save values in lambdas_df dataframe
        #temp_stages_df.at[np.where( (temp_stages_df["prop"] == prop) & (temp_stages_df["stage"] == stage))[0][0], "value"] = value

    #changing natural lengths of springs
    for index, row in springs_df.iterrows():
        if k_type == 'k_c':

            lambda_alpha = get_lambda(balls_df.iloc[row['ball1']], DV_present = DV_present, 
                                              lambda_anisotropic_obj = lambda_anisotropic_obj, lambda_isotropic_obj = lambda_isotropic_obj,
                                              inDV_lambda_anisotropic_obj = inDV_lambda_anisotropic_obj, inDV_lambda_isotropic_obj = inDV_lambda_isotropic_obj,
                                              #lambda_height_obj = lambda_height_obj,inDV_lambda_height_obj = inDV_lambda_height_obj,
                                              )
            lambda_beta = get_lambda(balls_df.iloc[row['ball2']], DV_present = DV_present, 
                                              lambda_anisotropic_obj = lambda_anisotropic_obj, lambda_isotropic_obj = lambda_isotropic_obj,
                                              inDV_lambda_anisotropic_obj = inDV_lambda_anisotropic_obj, inDV_lambda_isotropic_obj = inDV_lambda_isotropic_obj,
                                              #lambda_height_obj = lambda_height_obj,inDV_lambda_height_obj = inDV_lambda_height_obj,
                                              )            

        lambda_alpha_beta = 0.5*(lambda_alpha + lambda_beta)
        spring_vector = np.array([ row['x1'] - row['x2'], row['y1'] - row['y2'], row['z1'] - row['z2']])
        virtual_spring_vector = np.matmul(lambda_alpha_beta, spring_vector)
        l0 = np.linalg.norm(virtual_spring_vector)
        springs_df.loc[index, 'l0_stage_'+str(stage)] = l0

    springs_df['l0_target_final'] = springs_df['l0_stage_'+str(stage)]
    springs_df['l0_target'] = springs_df['l0_target_final']
    springs_df['l1_initial'] = springs_df['l1']

    #making plots
    #title = k_type + ', thickness:' + str(thickness)
    #plot_shell(balls_df, springs_df, x = 'y', y = 'z', filename = dirname + 'sim_output/DV_parallel_' + str(stage) + '.pdf', cbar_name = r'$\frac{l_{o}}{l_{initial}}$', title = title)
    #plot_shell(balls_df, springs_df, x = 'x', y = 'z', filename = dirname + 'sim_output/DV_across_' + str(stage) + '.pdf', cbar_name = r'$\frac{l_{o}}{l_{initial}}$', title = title)
    #plot_shell(balls_df, springs_df, x = 'x', y = 'y', filename = dirname + 'sim_output/top_view_' + str(stage) + '.pdf', cbar_name = r'$\frac{l_{o}}{l_{initial}}$', title = title)

balls_df.to_csv(dirname + 'sim_output/init_balls.csv', index = False)
springs_df.to_csv(dirname + 'sim_output/init_springs.csv', index = False)
#save lambdas
#remove rows not needed
#temp_stages_df = temp_stages_df.loc[temp_stages_df["prop"].isin(["inDV_lambda_isotropic_coeffs", "inDV_lambda_anisotropic_coeffs", "lambda_isotropic_coeffs", "lambda_anisotropic_coeffs"])].reset_index(drop = True)
#temp_stages_df.to_pickle("input_lambda_df.pkl")
#temp_stages_df.to_csv('input_lambda_df.csv', index = False)

# doing simulation

vtk_filename = dirname + 'sim_output/' + k_type
vtk_filename = vtk_filename + '_thickness_' + str(thickness)
vtk_filename = vtk_filename + '_'

#dfToVtk(balls_df,springs_df,filename=vtk_filename + '0.vtk', add_polygons = True)

springs_df['l0'] = springs_df['l1']


for i in range(len(stages)*nb_iterations):

    #print('growth step : ' + str(i))

    l0_stage_2 = "l0_stage_" + str(int(i/nb_iterations))
    if int(i/nb_iterations) == 0:
        l0_stage_1 = 'l1_initial'
    else:
        l0_stage_1 = "l0_stage_" + str(int(i/nb_iterations)-1)

    springs_df["l0_target"] = springs_df[l0_stage_2]
    springs_df['l0'] = springs_df['l0'] + (springs_df[l0_stage_2] - springs_df[l0_stage_1])/nb_iterations


    [balls_df, springs_df] = initialize_cpp_simulation(balls_df, springs_df, dt = dt,
                                                       csv_t_save = 10000, tol = tol, path = dirname)

    #save the first timepoint
    if i == 0:
        os.makedirs(dirname + 'sim_output/growth_0', exist_ok=True)
        temp_df = pd.DataFrame(balls_df[["x","y","z"]].values, columns = ['x[0]','x[1]','x[2]'])
        temp_df.to_csv(dirname + 'sim_output/growth_0/final_0_0.csv', index = False)
        dfToVtk(balls_df,springs_df,filename=dirname + 'sim_output/stage_0.vtk', add_polygons = True)

    if os.path.exists('files/'):
        shutil.rmtree('files/')

    #print('$$$$$$$ Running openfpm $$$$$$$')
    os.system("cd " + dirname + " && source ~/openfpm_vars && make && grid")
    #print('$$$$ Exit OpenFPM $$$$')
    
    #delete everything that is not needed
    #import glob, os
    for f in glob.glob(dirname + "files/*.vtk"):
        os.remove(f)
    for f in glob.glob(dirname + "files/Spring_*.csv"):
        os.remove(f)

    #rename the folder
    sim_folder = dirname + 'sim_output/growth_'+str(i+1)
    if os.path.exists(sim_folder):
        shutil.rmtree(sim_folder)
    shutil.move(dirname + 'files/', sim_folder)

    #open the last csv file
    df = pd.read_csv(sim_folder + '/final_0_0.csv')

    #update the position of the balls and springs
    balls_df['x'] = df['x[0]']
    balls_df['y'] = df['x[1]']
    balls_df['z'] = df['x[2]']
    springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']]) #this line probably not required

    #dfToVtk(balls_df,springs_df,filename=vtk_filename + str(i+1) + '.vtk', add_polygons = True)

    #save teh csv file if last file of stage
    if (i+1)%nb_iterations == 0:
        temp_stages_df = stages_df[stages_df["stage"] == int(i/nb_iterations)]
        stage_name = temp_stages_df.iloc[0]["stage_name"]
        balls_df.to_csv(dirname + 'sim_output/' + stage_name + '.csv')
        dfToVtk(balls_df,springs_df,filename=dirname + 'sim_output/stage_' + str(int(i/nb_iterations)+1) + '.vtk', add_polygons = True)





