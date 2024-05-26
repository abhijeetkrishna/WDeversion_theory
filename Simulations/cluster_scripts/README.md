# Running the simulations

In order to run the simulations, you first need to have [OpenFPM](http://openfpm.mpi-cbg.de/) installed. 

To start running the simulation, you need to run the bash script [array_wd.sh](array_wd.sh). This is currently made for running an array job on the MPI-CBG cluster but can be easily changed to run locally or other clusters. 

The bash script first runs the python script [map_index_wd.py](map_index_wd.py) which makes a csv with a name like [map_index_49421574.csv](map_index_49421574.csv). This consists of a list of all simulations to run with all the different parameters in a single row. 

Then the bash script runs the python script [array_wd.py](array_wd.py). This script takes the csv file [map_index_49421574.csv](map_index_49421574.csv), sets the parameters. Then it uses a pre-existing spherical mesh pickled [file](../meshes/icosaspherical_cap/pickle/IcoSph_mesh_refine_factor_30.pkl) to prepare a spherical cap mesh. On this mesh, the spontaneous strains are input using the lambda functions specified in the [csv file](input_lambda_df.csv).

After this, the OpenFPM simulation is run using the [C++ script](main.cpp).

After running the simulation, the bash file runs python scripts [postprocess.py](postprocess.py) and [residual_calculation.py](residual_calculation.py) which are used to measure the output shape and calculate the residual strains respectively.

After all the simulations for different parameter combinations have been run, the bash script [sweep_analysis.sh](sweep_analysis.sh) runs the python script [sweep_analysis.py](sweep_analysis.py) which averages the outputs by combining simulations for different random seeds and parameter combinations. The csv files that are created by this script are used to plot the results in the notebook [readModelResults.ipynb](../notebooks/readModelResults.ipynb).