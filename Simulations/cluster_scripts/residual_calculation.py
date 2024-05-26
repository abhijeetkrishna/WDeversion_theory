#from methods import *
from methods import *

map_index_dest= sys.argv[1] #"map_index_0.csv"
task_id= int(sys.argv[2]) #0

#functions
def compute_F(vertex_ID):
    #choose a vertex for which to evaluate the deformation gradient tensor
    #vertex_ID = 11
    #initial position of the vertex
    vertex_init = balls_init.loc[balls_init['ID'] == vertex_ID][["x", "y", "z"]].values.flatten()
    #final position of the vertex
    vertex_final = balls_init.loc[balls_final['ID'] == vertex_ID][["x_final", "y_final", "z_final"]].values.flatten()
    #find the neighbours of the vertex
    neighbours = springs_init.loc[springs_init['ball1'] == vertex_ID, 'ball2'].tolist() + springs_init.loc[springs_init['ball2'] == vertex_ID, 'ball1'].tolist()
    #make a dataframe of the neighbours
    neighbours_df = balls_init.loc[balls_init['ID'].isin(neighbours)][["ID", "x"]]
    #get the initial positions of the neighbours
    neighbours_df["x"] = neighbours_df.apply(lambda row: balls_init.loc[balls_init['ID'] == row["ID"]][["x", "y", "z"]].values.flatten(), axis = 1)
    #get the final positions of the neighbours
    neighbours_df["final_x"] = neighbours_df.apply(lambda row: balls_final.loc[balls_final['ID'] == row["ID"]][["x", "y", "z"]].values.flatten(), axis = 1)
    #rename
    neighbours_df = neighbours_df.reset_index(drop = True).rename(columns={'x': 'init_x', 'final_x': 'final_x'})
    #get edge vectors
    neighbours_df["init_edge"] = neighbours_df.apply(lambda row: (row["init_x"] - vertex_init).flatten(), axis = 1)
    neighbours_df["final_edge"] = neighbours_df.apply(lambda row: (row["final_x"] - vertex_final).flatten(), axis = 1)
    #find a least squares solution for a deformation gradient tensor
    X = np.stack(neighbours_df["init_edge"].values)
    x = np.stack(neighbours_df["final_edge"].values)
    lstsq_output = np.linalg.lstsq(X,x,rcond = None)
    F = lstsq_output[0].T #careful here, the output is transposed
    #if symmetrise : F = 0.5*(F + F.T)
    residuals = lstsq_output[1]
    #get fit x - matrix multiplication of F with X
    x_fit = np.dot(F,X.T).T #np.dot(X,F)
    #add fit to dataframe
    neighbours_df["final_edge_fit"] = x_fit.tolist()
    neighbours_df["final_x_fit"] = neighbours_df.apply(lambda row: (row["final_edge_fit"] + vertex_final).flatten(), axis = 1)
    #give the dataframe a name
    neighbours_df.vertex_ID = vertex_ID
    neighbours_df.vertex_init = vertex_init
    neighbours_df.vertex_final = vertex_final
    neighbours_df.F = F
    neighbours_df.residuals = residuals

    return neighbours_df

def iterate_compute_F(row):

    vertex_ID = row.ID
    #get neighbours and compute F in cartesian coordinates
    neighbours_df = compute_F(vertex_ID)
    #Get F in the local coordinate system
    F = neighbours_df.F
    #get the local coordinate system
    e_R = row[["e_R_x", "e_R_y", "e_R_z"]].values.flatten()
    e_phi = row[["e_phi_x", "e_phi_y", "e_phi_z"]].values.flatten()
    e_h = row[["e_h_x", "e_h_y", "e_h_z"]].values.flatten()
    e_local = [e_R, e_phi, e_h]
    #compute components along e_local
    F_aligned = np.zeros((3,3))
    for l in range(3):
        for m in range(3):
            e_local_1 = e_local[l]
            e_local_2 = e_local[m]
            for i in range(3):
                for j in range(3):
                    F_aligned[l,m] += e_local_1[i]*e_local_2[j]*F[i,j]
            #F_aligned[l,m] = np.einsum('i,j,ij', e_local[l], e_local[m], F) #this is giving an error in cluster
    neighbours_df.F_aligned = F_aligned
    #save to row
    row["F_rr"]  = F_aligned[0,0]
    row["F_phiphi"]  = F_aligned[1,1]
    row["F_hh"]  = F_aligned[2,2]
    row["F_rphi"]  = F_aligned[0,1]
    row["F_phir"]  = F_aligned[1,0]
    row["F_rh"]  = F_aligned[0,2]
    row["F_hr"]  = F_aligned[2,0]
    row["F_phih"]  = F_aligned[1,2]
    row["F_hphi"]  = F_aligned[2,1]
    for i in range(3):
        for j in range(3):
            row["F_" + str(i) + str(j)] = F[i,j]

    return(row)

#reading the map_index csv
print('reading the map_index csv')
map_index=pd.read_csv(map_index_dest)
input_var=list(map_index.columns) #variables for which values are being imported
for i in range(len(input_var)):
    if input_var[i] == 'folder_name':
        dirname = map_index.loc[task_id,str(input_var[i])]
    if input_var[i] == 'thickness':
        thickness = map_index.loc[task_id,str(input_var[i])]

#set the stage names
stage_init = "wL3"
stage_final_list = ["0hAPF", "2hAPF", "4hAPF"]
balls_w_deformation_gradient = pd.DataFrame()

for stage_final in stage_final_list:

    #stage_final = "4hAPF"

    #choose the initial timepoint
    balls_init = pd.read_csv(dirname + "sim_output/init_balls.csv")
    springs_init = pd.read_csv(dirname + "sim_output/init_springs.csv")
    #choose the final timepoint
    balls_final = balls_init.copy(deep = True)
    springs_final = springs_init.copy(deep = True)
    #final_pos = pd.read_csv(dirname + "sim_output/wL3 to 4hAPF.csv")
    final_pos = pd.read_csv(dirname + "sim_output/" + stage_init + " to " + stage_final + ".csv")
    balls_final[['x', 'y', 'z']] =  final_pos[['x', 'y', 'z']]
    springs_final = update_springs(springs_final, balls_final[['x', 'y', 'z']])
    #extract middle mesh
    #[balls_init, springs_init] = extract_thin_mesh(balls_init, springs_init, which = "middle", reindex = True)
    #[balls_final, springs_final] = extract_thin_mesh(balls_final, springs_final, which = "middle", reindex = True)
    balls_init[["x_final", "y_final", "z_final"]] = balls_final[["x", "y", "z"]]

    #compute the deformation gradient tensor for all vertices
    balls_init = balls_init.apply(iterate_compute_F, axis = 1)
    balls_init["F_iso"] = np.sqrt(balls_init["F_rr"]*balls_init["F_phiphi"])
    balls_init["F_aniso"] = np.sqrt(balls_init["F_rr"]/balls_init["F_phiphi"])

    #save
    #balls_init.to_csv(dirname + "lattice_csv/balls_with_fit_lambda.csv", index = False)

    # read
    #balls_init = pd.read_csv(dirname + "lattice_csv/balls_with_fit_lambda.csv")

    #read the lambda input values
    lambda_input = pd.read_pickle("input_lambda_df.pkl")

    lambda_iso_DV = np.poly1d(lambda_input.loc[(lambda_input["prop"] == "inDV_lambda_isotropic_coeffs") & (lambda_input["stage_final"] == stage_final), "value"].values[0])
    lambda_aniso_DV = np.poly1d(lambda_input.loc[(lambda_input["prop"] == "inDV_lambda_anisotropic_coeffs") & (lambda_input["stage_final"] == stage_final), "value"].values[0]) #np.poly1d([1])#
    lambda_iso_outDV = np.poly1d(lambda_input.loc[(lambda_input["prop"] == "lambda_isotropic_coeffs") & (lambda_input["stage_final"] == stage_final), "value"].values[0])
    lambda_aniso_outDV = np.poly1d(lambda_input.loc[(lambda_input["prop"] == "lambda_anisotropic_coeffs") & (lambda_input["stage_final"] == stage_final), "value"].values[0]) #np.poly1d([1])#

    def get_input_lambda_values(row):

        if row["DV_bool"] == 1:
            row["lambda_iso"] = lambda_iso_DV(row["pathlength_scaled"])
            row["lambda_aniso"] = lambda_aniso_DV(row["pathlength_scaled"])
        else:  
            row["lambda_iso"] = lambda_iso_outDV(row["pathlength_scaled"])
            row["lambda_aniso"] = lambda_aniso_outDV(row["pathlength_scaled"])
        row["lambda_rr"] = row["lambda_iso"]*row["lambda_aniso"]
        row["lambda_phiphi"] = row["lambda_iso"]/row["lambda_aniso"]
        return(row)

    balls_init = balls_init.apply(get_input_lambda_values, axis = 1)

    #compute the lambda residual isotropic and anisotropic
    balls_init["lambda_res_rr"] = balls_init["F_rr"]/balls_init["lambda_rr"]
    balls_init["lambda_res_phiphi"] = balls_init["F_phiphi"]/balls_init["lambda_phiphi"]
    balls_init["lambda_res_iso"] = np.sqrt(balls_init["lambda_res_rr"]*balls_init["lambda_res_phiphi"])
    balls_init["lambda_res_aniso"] = np.sqrt(balls_init["lambda_res_rr"]/balls_init["lambda_res_phiphi"])

    #add stage name
    balls_init["stage_final"] = stage_final

    #append balls_init to balls_w_deformation_gradient
    balls_w_deformation_gradient = pd.concat([balls_w_deformation_gradient, balls_init], axis = 0, ignore_index = True)

#drop some columns that are not needed 
cols_to_drop = ["neighbours", "spring1", "spring2",
                "e_h_x", "e_h_y", "e_h_z", "e_R_x", "e_R_y", "e_R_z", "e_phi_x", "e_phi_y", "e_phi_z",
                ]
balls_w_deformation_gradient = balls_w_deformation_gradient.drop(columns = cols_to_drop)

#save
balls_w_deformation_gradient.to_csv(dirname + "sim_output/balls_with_deformation_gradient.csv", index = False)



