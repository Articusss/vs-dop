module VsDop
    import IterTools as itr
    using Random
    include("AcceleratedDubins.jl")
    include("Visual.jl")
    include("Helper.jl")
    include("Vns.jl")

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

            max_speed_discrete = min(v_params.v_max, AcceleratedDubins.speed_by_radius(v_params.r_max))
            speeds = collect(range(v_params.v_min, max_speed_discrete, num_speeds))
            headings = collect(range(0, 2 * pi, num_headings))
            radii = AcceleratedDubins.radii_samples_exp(v_params.r_min, v_params.r_max, radii_samples)
            graph = fill(-1, (num_locations, num_locations, num_speeds, num_speeds, num_headings, num_headings))

            return new(graph, v_params, num_speeds,num_headings,radii_samples,speeds,headings,radii)
        end

        function DOPGraph(num_speeds, num_headings, radii_samples, num_locations, v_params)
            max_speed_discrete = min(v_params.v_max, AcceleratedDubins.speed_by_radius(v_params.r_max))
            speeds = collect(range(v_params.v_min, max_speed_discrete, num_speeds))
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

    mutable struct Results
        scores::Vector{Float64}
        travel_times::Vector{Float64}
        timestamps::Vector{Float64}
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

    function vs_dop(op_params, max_iterations::Int64 = 2000, verbose::Bool = false)
        initial_seq, ini_time, ini_score = greedy_solution(op_params)
        if verbose
            println("Initial Score: ", ini_score, ". Initial time: ", ini_time)
        end

        vns_seq, time, score = variable_neighborhood_search(op_params, initial_seq, max_iterations, verbose)

        config = Helper.shortest_configuration_by_sequence(op_params, vns_seq)
        path = Helper.retrieve_path(op_params, config)

        return path, config, score, time
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
        
        @Threads.threads for node_i in 1:num_locations
            for node_f in 1:num_locations
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
        end

        return graph_params
    end

    function greedy_solution(op_params)
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

        return sequence, seq_time, Helper.get_score(op_params, sequence)
    end

    #BASE VNS
    function variable_neighborhood_search(op_params, initial_sequence::Vector{Int64}, max_iterations::Int64, verbose::Bool = false)
        best_time = Helper.shortest_time_by_sequence(op_params, initial_sequence)
        best_sequence = deepcopy(initial_sequence)
        best_score = Helper.get_score(op_params, initial_sequence)
        len = length(initial_sequence)

        #For dynamic programming
        distances = fill((0., false), (length(op_params.coordinates) + 1, op_params.graph.num_speeds, op_params.graph.num_headings))
        distances_cand = fill((0., false), (length(op_params.coordinates) + 1, op_params.graph.num_speeds, op_params.graph.num_headings))
        
        for i in 1:max_iterations
            if verbose
                println((i, best_time, best_score))
            end

            l = 1
            while l <= 2 
                #Shake
                local_sequence = Vns.shake(deepcopy(best_sequence), l)
                local_time, depot_pos = Helper.shortest_time_by_sequence(op_params, local_sequence, distances)
                local_score = Helper.get_score(op_params, local_sequence)

                #Search
                for j in 1:len^2
                    search, change_pos = Vns.search(deepcopy(local_sequence), l)

                    #If change is done above depot, nothing has happened and can be skipped
                    if (change_pos > depot_pos) 
                        j -= 1 #Run again
                        continue
                    end
                    
                    search_score = Helper.get_score(op_params, search)

                    #Only check time if solutions score higher OR local_sequence is not valid
                    if (search_score >= local_score) || (local_time > op_params.tmax) 
                        #"Borrow" distances from local sequence
                        distances_cand[1:change_pos-1, :, :], distances[1:change_pos-1, :, :] = distances[1:change_pos-1, :, :], distances_cand[1:change_pos-1, :, :]
                        search_time, cand_depot_pos = Helper.shortest_time_by_sequence(op_params, search, distances_cand, change_pos - 1)
                        
                        if (search_score > local_score && search_time <= op_params.tmax) || (search_score == local_score && search_time < local_time) 
                            local_time = search_time
                            local_sequence = search
                            local_score = search_score
                            depot_pos = cand_depot_pos
                            
                            #Fully swap
                            distances, distances_cand = distances_cand, distances
                        else
                            #Return distances
                            distances_cand[1:change_pos-1, :, :], distances[1:change_pos-1, :, :] = distances[1:change_pos-1, :, :], distances_cand[1:change_pos-1, :, :]
                        end
                    end
                end

                #Higher score is found, only thing that matters is TMAX
                if local_time <= op_params.tmax && local_score > best_score
                    best_time = local_time
                    best_sequence = local_sequence
                    best_score = local_score
                    l = 1
                #Equal score is found, swap if time is better
                elseif local_score == best_score && local_time < best_time
                    best_time = local_time
                    best_sequence = local_sequence
                    best_score = local_score
                    l = 1
                else
                    l += 1
                end
            end
        end

        return best_sequence, best_time, best_score
    end

    #VNS WITH BENCHMARKING TOOLS
    function variable_neighborhood_search_benchmark(op_params, initial_sequence::Vector{Int64}, max_time::Int64, results, t0)
        best_time = Helper.shortest_time_by_sequence(op_params, initial_sequence)
        best_sequence = deepcopy(initial_sequence)
        best_score = Helper.get_score(op_params, initial_sequence)
        len = length(initial_sequence)

        #For dynamic programming
        distances = fill((0., false), (length(op_params.coordinates) + 1, op_params.graph.num_speeds, op_params.graph.num_headings))
        distances_cand = fill((0., false), (length(op_params.coordinates) + 1, op_params.graph.num_speeds, op_params.graph.num_headings))
        
        initial_time = time()
        while time() - initial_time <= max_time
            #Store info
            if (time() - t0) - results.timestamps[end] > 5
                push!(results.scores, best_score)
                push!(results.travel_times, best_time)
                push!(results.timestamps, time() - t0)
            end

            l = 1
            while l <= 2 
                #Shake
                local_sequence = Vns.shake(deepcopy(best_sequence), l)
                local_time, depot_pos = Helper.shortest_time_by_sequence(op_params, local_sequence, distances)
                local_score = Helper.get_score(op_params, local_sequence)

                #Search
                for j in 1:len^2
                    search, change_pos = Vns.search(deepcopy(local_sequence), l)

                    #If change is done above depot, nothing has happened and can be skipped
                    if (change_pos > depot_pos) 
                        j -= 1 #Run again
                        continue
                    end
                    
                    search_score = Helper.get_score(op_params, search)

                    #Only check time if solutions score higher OR local_sequence is not valid
                    if (search_score >= local_score) || (local_time > op_params.tmax) 
                        #"Borrow" distances from local sequence
                        distances_cand[1:change_pos-1, :, :], distances[1:change_pos-1, :, :] = distances[1:change_pos-1, :, :], distances_cand[1:change_pos-1, :, :]
                        search_time, cand_depot_pos = Helper.shortest_time_by_sequence(op_params, search, distances_cand, change_pos - 1)
                        
                        if (search_score > local_score && search_time <= op_params.tmax) || (search_score == local_score && search_time < local_time) 
                            local_time = search_time
                            local_sequence = search
                            local_score = search_score
                            depot_pos = cand_depot_pos
                            
                            #Fully swap
                            distances, distances_cand = distances_cand, distances
                        else
                            #Return distances
                            distances_cand[1:change_pos-1, :, :], distances[1:change_pos-1, :, :] = distances[1:change_pos-1, :, :], distances_cand[1:change_pos-1, :, :]
                        end
                    end
                end

                #Higher score is found, only thing that matters is TMAX
                if local_time <= op_params.tmax && local_score > best_score
                    best_time = local_time
                    best_sequence = local_sequence
                    best_score = local_score
                    l = 1
                #Equal score is found, swap if time is better
                elseif local_score == best_score && local_time < best_time
                    best_time = local_time
                    best_sequence = local_sequence
                    best_score = local_score
                    l = 1
                else
                    l += 1
                end
            end
        end

        return best_sequence, best_time, best_score
    end

    function benchmark(op, num_instances, exec_time)
        r = []
        seqs = []
        best_scores = Vector{Float64}()
        best_times = Vector{Float64}()
        for _ in 1:num_instances
            push!(r, Results([], [], []))
            push!(seqs, [])
            push!(best_scores, -1)
            push!(best_times, -1)
        end
        Threads.@threads for i in 1:num_instances
            #Run algo
            t0 = time()
            ini_seq, t = greedy_solution(op)
            push!(r[i].scores, VsDop.Helper.get_score(op, ini_seq))
            push!(r[i].travel_times, t)
            push!(r[i].timestamps, time() - t0)

            #Run vns
            seq, t, score = variable_neighborhood_search_benchmark(op, ini_seq, exec_time, r[i], t0)
            seqs[i] = seq
            best_scores[i] = score
            best_times[i] = t
            
            push!(r[i].scores, score)
            push!(r[i].travel_times, t)
            push!(r[i].timestamps, time() - t0)
        end

        return r, seqs, best_scores, best_times, findmax(best_scores)[2]
    end

    function bulk_benchmark(configurations, run_per_config, base_path, instance_list, exec_time)
        #configurations is a vector of (num_speeds, num_heading angles)
        #TODO remove temp code
        vehicle_params = [VehicleParameters(30., 30., -3., 2., 65.7, 264.2), VehicleParameters(66.98, 66.98, -3., 2., 65.7, 264.2)]
        for inst in instance_list
            points, scores, depots, tmax = Helper.read_op_file("$base_path/$inst")
            for v_params in vehicle_params
                graph_params = DOPGraph(2, 8, 8, length(points), v_params)
                op = OpParameters(compute_trajectories(points, graph_params), points, scores, depots, tmax)
                res, best_seqs, best_scores, best_times, best_idx = benchmark(op, run_per_config, exec_time)

                dir = "results/CONST$(v_params.v_max)v"
                isdir(dir) || mkdir(dir)
                #Write to file
                
                f = open("$dir/$inst", "w")
                println(f, "BEST_SCORE = $(maximum(best_scores))")
                println(f, "BEST_TIME = $(best_times[best_idx])")
                println(f, "BEST_SEQ = $(best_seqs[best_idx])")
                println(f, "SCORES = $best_scores")
                println(f, "TIMES = $best_times")
                println(f, "R_BEST_SCORES = $(res[best_idx].scores)")
                println(f, "R_BEST_TRAVELTIMES = $(res[best_idx].travel_times)")
                println(f, "R_BEST_TIMESTAMPS = $(res[best_idx].timestamps)")
                close(f)
                
            end
        end
    end
end # module VsTsp