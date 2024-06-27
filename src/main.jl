include("VsDop.jl")

#"Default" vehicle parameters, you can also make your own with VehicleParameters struct
vehicle_params = VsDop.get_cessna172_params()
#Get our instance parameters, graph and graph parameters.
#Number of discrete speeds and heading angles can be specified in optional graph_params parameter.
dop = VsDop.get_graph("../instances/random.op", vehicle_params);
#Execute the algorithm. Max iterations and verbosity can be specified by the last two arguments.
path, config, score, time = VsDop.vs_dop(dop, 100, true)
#Plot resulting path/speeds
VsDop.Visual.plot_full_path(dop, path)
VsDop.Visual.plot_full_speeds(dop, path)