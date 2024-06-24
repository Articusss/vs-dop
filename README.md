# vs-dop
Implementation of the Variable-Speed Dubins Orienteering Problem

## Basic usage

The project can be used by importing into your REPL or into your Julia script by running `using VsDop`, then used as follows:

```julia
using VsDop
#"Default" vehicle parameters, you can also make your own with VehicleParameters struct
vehicle_params = VsDop.get_cessna172_params()
#Get our instance parameters, graph and graph parameters.
#Number of discrete speeds and heading angles can be specified in optional graph_params parameter.
dop = VsDop.get_graph("./instances/random.op", vehicle_params);
#Execute the algorithm. Max iterations and verbosity can be specified by the last two arguments.
path, config, score, time = VsDop.vs_dop(dop, max_iterations=1000, true)
#Plot resulting path/speeds
VsDop.Visual.plot_full_path(dop, path)
VsDop.Visual.plot_full_speeds(dop, path)
```

When running the project through the REPL, make sure to activate the project by pressing `]` to switch to the package manager and running 

```julia
activate .
```