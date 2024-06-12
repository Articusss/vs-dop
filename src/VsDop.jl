module VsDop
    import IterTools as itr
    import Plots as plt
    using Random
    include("AcceleratedDubins.jl")
    include("Visual.jl")
    include("Helper.jl")
    include("Vns.jl")
    include("FastReject.jl")

    struct VehicleParameters
        v_min::Float64
        v_max::Float64
        a_min::Float64
        a_max::Float64
        r_min::Float64
        r_max::Float64
    end

    struct DOPGraph
        graph::Array{Float64, 6}
        vehicle_params::VehicleParameters
        num_speeds::Int64
        num_headings::Int64
        radii_samples::Int64
        speeds::Vector{Float64}
        headings::Vector{Float64}
        radii::Vector{Float64}
        function DOPGraph(num_locations, v_params)
            num_speeds = 4
            num_headings = 8
            radii_samples = 8

            speeds = collect(range(v_params.v_min, AcceleratedDubins.speed_by_radius(v_params.r_max), num_speeds))
            headings = collect(range(0, 2 * pi, num_headings))
            radii = AcceleratedDubins.radii_samples_exp(v_params.r_min, v_params.r_max, radii_samples)
            graph = fill(-1, (num_locations, num_locations, num_speeds, num_speeds, num_headings, num_headings))

            return new(graph, v_params, num_speeds,num_headings,radii_samples,speeds,headings,radii)
        end

        function DOPGraph(num_speeds, num_headings, radii_samples, num_locations, v_params)
            speeds = collect(range(v_params.v_min, AcceleratedDubins.speed_by_radius(v_params.r_max), num_speeds))
            headings = collect(range(0, 2 * pi, num_headings))
            radii = AcceleratedDubins.radii_samples_exp(v_params.r_min, v_params.r_max, radii_samples)
            graph = fill(-1, (num_locations, num_locations, num_speeds, num_speeds, num_headings, num_headings))

            return new(graph, v_params, num_speeds, num_headings, radii_samples, speeds, headings, radii)
        end
    end

    mutable struct OpParameters
        graph::DOPGraph
        coordinates::Array{Tuple{Float64, Float64}, 1}
        scores::Vector{Float64}
        depots::Vector{Int64}
        tmax::Float64
    end

    function get_graph(instance_path::String, vehicle_params, graph_params = nothing)
        points, scores, depots, tmax = Helper.read_op_file(instance_path)
        
        if graph_params === nothing
            graph_params = DOPGraph(length(points), vehicle_params)
        end

        return OpParameters(compute_trajectories(points, graph_params), points, scores, depots, tmax)
    end

    function get_cessna172_params()
        return VehicleParameters(30., 67., -3., 2., 65.7, 264.2)
    end

    function vns_tsp(graph::Array{Float64, 6}, instance_path::String, vehicle_params, graph_params, max_iterations::Int64 = 2000, verbose::Bool = false)
        _, initial_seq, _ = greedy_solution(graph)
        vns_seq, time = variable_neighborhood_search_fast(graph, initial_seq, max_iterations, num_windows, verbose)

        points = Helper.read_tsp_file(instance_path)
        _, config = Helper.shortest_configuration_by_sequence(graph, vns_seq)
        path, _ = Helper.retrieve_path(points, config, params, speeds, headings, radii)

        return path, points, config
    end

    #TODO not calculate matrix diagonals (distance from the point to itself) can speed up
    #Also maybe take advantage of paths being borderline simmetric (only issue is radii and acceleration)
    function compute_trajectories(locations::Array{Tuple{Float64, Float64}, 1}, graph_params)
        #Build a 6-dimensional graph (starting node, ending node, starting speed, ending speed, starting heading angle, ending heading angle)
        num_locations = length(locations)
        num_speeds = graph_params.num_speeds
        num_headings = graph_params.num_headings
        
        graph = graph_params.graph

        speeds = graph_params.speeds
        headings = graph_params.headings
        radii = graph_params.radii

        # v_min v_max a_max -a_min
        vehicle_params = graph_params.vehicle_params
        params = [vehicle_params.v_min, vehicle_params.v_max, vehicle_params.a_max, -vehicle_params.a_min]
        
        for (node_i, node_f) in itr.product(1:num_locations, 1:num_locations)
            for (v_i, v_f) in itr.product(1:num_speeds, 1:num_speeds)
                for (h_i, h_f) in itr.product(1:num_headings, 1:num_headings)
                    start::Vector{Float64} = [locations[node_i][1], locations[node_i][2], headings[h_i]]
                    stop::Vector{Float64} = [locations[node_f][1], locations[node_f][2], headings[h_f]]
                    path, time, _ = AcceleratedDubins.fastest_path(start, stop, radii, params, [speeds[v_i], speeds[v_f]])
                    #Set to edge, note that we can't take advantage of simmetry because acceleration max/min is not necessarily the same
                    graph[node_i,node_f,v_i,v_f,h_i,h_f] = path === nothing ? Inf : time
                end
            end
        end

        return graph_params
    end

    function greedy_solution(op_params)
        num_headings = op_params.graph.num_headings
        num_speeds = op_params.graph.num_speeds
        
        sequence = [op_params.depots[1], op_params.depots[1]]
        seq_time = Helper.shortest_time_by_sequence(op_params, sequence)
        to_add = Set{Int64}(1:length(op_params.coordinates))
        delete!(to_add, op_params.depots[1])

        added = true
        while added
            added = false
            best_val = Inf
            best_position = -1
            best_point = -1
            new_time = -1

            #Find cheapest insertion -> Try to insert each node between every point in solution, select the best position and insert
            for candidate in to_add
                for pos in 2:length(sequence)
                    cand_sequence = insert!(deepcopy(sequence), pos, candidate)
                    time = Helper.shortest_time_by_sequence(op_params, cand_sequence)
                    difference = (time - seq_time) / op_params.scores[candidate]
                    if difference < best_val && time <= op_params.tmax
                        added = true
                        new_time = time
                        best_val = difference
                        best_point = candidate
                        best_position = pos
                    end
                end
            end

            if added
                #Insert best insertion
                delete!(to_add, best_point)
                insert!(sequence, best_position, best_point)
                seq_time = new_time
            end
        end

        #Randomly add rest of nodes
        sequence = vcat(sequence, shuffle(collect(to_add)))

        return sequence, seq_time
    end

    #BASE VNS -> no fast reject
    function variable_neighborhood_search(graph::Array{Float64, 6}, initial_sequence::Vector{Int64}, max_iterations::Int64, verbose::Bool = false)
        best_time = Helper.shortest_time_by_sequence(graph, initial_sequence)
        best_sequence = deepcopy(initial_sequence)
        len = length(initial_sequence)

        for i in 1:max_iterations
            if verbose
                println((i, best_time))
            end
            #Shake
            local_sequence = rand([Vns.path_exchange, Vns.path_move])(deepcopy(best_sequence))
            local_time = Helper.shortest_time_by_sequence(graph, local_sequence)

            #Search
            for j in 1:len^2
                search = Vns.rand([Vns.one_point_move, Vns.open_point_exchange])(deepcopy(local_sequence))
                search_time = Helper.shortest_time_by_sequence(graph, search)

                if search_time < local_time
                    local_time = search_time
                    local_sequence = search
                end
            end

            if local_time < best_time
                best_time = local_time
                best_sequence = local_sequence
            end
        end


        return best_sequence, best_time
    end

    #VNS with fast reject
    function variable_neighborhood_search_fast(graph::Array{Float64, 6}, initial_sequence::Vector{Int64}, max_iterations::Int64, num_windows::Int64, verbose::Bool = false)
        best_time = Helper.shortest_time_by_sequence(graph, initial_sequence)
        best_sequence = deepcopy(initial_sequence)
        len = length(initial_sequence)
        stored = fill(-1., (len,len))

        for i in 1:max_iterations
            if verbose
                println((i, best_time))
            end
            #Shake
            local_sequence = rand([Vns.path_exchange, Vns.path_move])(deepcopy(best_sequence))
            local_costs::Vector{Float64} = fill(-1., num_windows)
            
            #Search
            for j in 1:len^2
                search = Vns.rand([Vns.one_point_move, Vns.open_point_exchange])(deepcopy(local_sequence))

                reject, costs = FastReject.fast_reject(graph, search, local_sequence, local_costs, num_windows, stored, best_time)
                if !reject
                    local_sequence = search
                    local_costs = costs
                end
            end

            local_time = Helper.shortest_time_by_sequence(graph, local_sequence)
            if local_time < best_time
                best_time = local_time
                best_sequence = local_sequence
            end
        end


        return best_sequence, best_time
    end
end # module VsTsp