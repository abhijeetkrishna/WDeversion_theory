# each simulation has a file called sim_output/balls_with_deformation_gradient.csv
# we want to combine all of them into one file

from methods import *
from sklearn.linear_model import LinearRegression


def lin_fit_w_sd(x_hat, x_data, y_data, sd_window=0.1, return_coef=False):
    # x_hat is the array of x values for which to evaluate the fit
    # x_data is the array of x values of the data
    # y_data is the array of y values of the data
    # sd_window is the window size for computing the standard deviation
    # get the standard deviation of the data in a sliding window
    sd = np.zeros(len(x_hat))
    for i in range(len(x_hat)):
        sd[i] = np.std(y_data[np.abs(x_data - x_hat[i]) < sd_window])
    # get the fit
    reg = LinearRegression().fit(x_data, y_data)
    # get the fit values
    y_hat = reg.predict(x_hat)
    if return_coef:
        coefs = np.append(reg.coef_, reg.intercept_)
        return (y_hat, sd, coefs)
    else:
        return (y_hat, sd)

#get list of folders for each thickness

param_comb_df = pd.read_csv(glob.glob('map_index_[0-9]*.csv')[0])
stage_final_array = ["0hAPF", "2hAPF", "4hAPF"]
thicknesses = np.sort(np.unique(param_comb_df["thickness"]))
roi_dict = {0: "outDV", 1: "DV"}

#genotype = 'ecadGFPnbG4'

#columns to fit
cols_to_fit = ['F_rr', 'F_phiphi', 'F_iso', 'F_aniso',
               'lambda_iso', 'lambda_aniso', 'lambda_rr', 'lambda_phiphi',
               'lambda_res_rr', 'lambda_res_phiphi', 'lambda_res_iso', 'lambda_res_aniso'
               ]

for thickness in thicknesses:

    combined_df = pd.DataFrame()
    combined_fit_df = pd.DataFrame()
    lambda_fit_coefs = pd.read_pickle("input_lambda_df.pkl")

    thickness_param_comb_df = param_comb_df.query("thickness == @thickness")
    folders = thickness_param_comb_df["folder_name"].values
    for folder in folders:
        temp_df = pd.read_csv(folder + "sim_output/balls_with_deformation_gradient.csv")
        combined_df = pd.concat([combined_df, temp_df], axis = 0, ignore_index = True)

    #we won't save the combined_df for now - since it must be about 300 MB
    #combined_df.to_csv("combined_df_" + str(thickness) + ".csv", index = False)

    #now we fit the lambdas for each stage
    #fit computed values to a straight line

    for stage in stage_final_array:

        #make a dataframe to store the fit
        balls_fit = pd.DataFrame(columns = ["pathlength_scaled", "DV_bool"] + cols_to_fit + [col + "_sd" for col in cols_to_fit])
        pathlength_scaled = np.linspace(0, 1, 100)
        balls_fit["DV_bool"] = np.append(np.zeros(len(pathlength_scaled)), np.ones(len(pathlength_scaled)))
        balls_fit["pathlength_scaled"] = np.append(pathlength_scaled, pathlength_scaled)
        balls_fit["stage_final"] = stage

        #fit columns one by one with a sliding window to compute the standard deviation
        sd_window = 0.1
        for col in cols_to_fit:
            #print("col", col)
            for DV_bool in [0,1]:
                #get the data
                balls_data = combined_df[(combined_df["DV_bool"] == DV_bool) & (combined_df["stage_final"] == stage)][["pathlength_scaled", col]]
                #remove NaNs
                balls_data = balls_data.dropna()
                #get the fit
                y_hat, sd, coefs = lin_fit_w_sd(pathlength_scaled.reshape(-1,1), balls_data["pathlength_scaled"].values.reshape(-1, 1), balls_data[col].values.reshape(-1, 1), sd_window = sd_window, return_coef = True)
                #save the fit
                balls_fit.loc[(balls_fit["DV_bool"] == DV_bool), col] = y_hat.flatten()
                balls_fit.loc[(balls_fit["DV_bool"] == DV_bool), col + "_sd"] = sd.flatten()
                #save the fit coefficients
                #concat the last row of the dataframe to itself
                lambda_fit_coefs = pd.concat([lambda_fit_coefs, lambda_fit_coefs.iloc[[-1]]], ignore_index=True)
                #lambda_fit_coefs = lambda_fit_coefs.append(lambda_fit_coefs.iloc[-1])
                #for the last row, update prop, roi and value columns
                lambda_fit_coefs.at[len(lambda_fit_coefs)-1, "prop"] = col
                lambda_fit_coefs.at[len(lambda_fit_coefs)-1, "value"] = coefs
                lambda_fit_coefs.at[len(lambda_fit_coefs)-1, "roi"] = roi_dict[DV_bool]
                lambda_fit_coefs.at[len(lambda_fit_coefs)-1, "stage_final"] = stage
                lambda_fit_coefs.at[len(lambda_fit_coefs)-1, "stage_name"] = "wL3 to " + stage
                lambda_fit_coefs.at[len(lambda_fit_coefs)-1, "stage"] = stage_final_array.index(stage)

        combined_fit_df = pd.concat([combined_fit_df, balls_fit], axis = 0, ignore_index = True)

    #save
    #balls_fit.to_csv(dirname + "lattice_csv/balls_with_linear_fit_lambda.csv", index = False)
    #lambda_fit_coefs.to_pickle(dirname + "analysis/def_grad_residual_fit_coefs.pkl")
    os.makedirs("analysis/", exist_ok = True)
    combined_fit_df.to_csv("analysis/balls_with_def_grad_residual_fit_thickness_" + str(thickness) + ".csv", index = False)
    lambda_fit_coefs.to_pickle("analysis/def_grad_residual_fit_coefs_thickness_" + str(thickness) + ".pkl")
