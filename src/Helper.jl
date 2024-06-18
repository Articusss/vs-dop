module Helper
    import IterTools as itr
    using ..AcceleratedDubins

    function read_op_file(filename::String)
        coordinates::Array{Tuple{Float64, Float64}, 1} = []
        values::Vector{Float64} = []
        depots::Vector{Int64} = []
        tmax = -1.

        open(filename) do file
            in_coordinates = false
            in_scores = false
            in_depots = false
            
            for line in eachline(file)
                tokens = split(line)
                if line == "NODE_COORD_SECTION"
                    in_coordinates = true
                elseif line == "NODE_SCORE_SECTION"
                    in_coordinates = false
                    in_scores = true
                elseif line == "DEPOT_SECTION"
                    in_scores = false
                    in_depots = true
                elseif tokens[1] == "COST_LIMIT"
                    tmax = parse(Float64, tokens[3])
                elseif line == "EOF"
                    break
                elseif in_coordinates
                    x_coord = parse(Float64, tokens[2])
                    y_coord = parse(Float64, tokens[3])
                    push!(coordinates, (x_coord, y_coord))
                elseif in_scores
                    score = parse(Float64, tokens[2])
                    push!(values, score)
                elseif in_depots
                    depot = parse(Int64, tokens[1])
                    if depot == -1
                        in_depots = false
                    else
                        push!(depots, depot)
                    end                  
                end
            end
        end

        return coordinates, values, depots, tmax
    end

    #Get path (DubinsPathR2) vector with configurations vector
    function retrieve_path(op_params, configurations::Vector{Tuple{Int64, Int64, Int64}})
        #(dubinspath, starting_speed, ending_speed)
        full_path::Vector{Tuple{AcceleratedDubins.DubinsPathR2, Float64, Float64}} = []

        locations = op_params.coordinates
        headings = op_params.graph.headings
        speeds = op_params.graph.speeds
        radii = op_params.graph.radii

        vehicle_params = op_params.graph.vehicle_params
        params = [vehicle_params.v_min, vehicle_params.v_max, vehicle_params.a_max, -vehicle_params.a_min]

        for i in 1:length(configurations)-1
            next = i + 1

            n_i, v_i, h_i = configurations[i]
            n_f, v_f, h_f = configurations[next]

            starting::Vector{Float64} = [locations[n_i][1], locations[n_i][2], headings[h_i]]
            ending::Vector{Float64} = [locations[n_f][1], locations[n_f][2], headings[h_f]]

            path, time, _ = AcceleratedDubins.fastest_path(starting, ending, radii, params, [speeds[v_i], speeds[v_f]])

            push!(full_path, (path, speeds[configurations[i][2]], speeds[configurations[next][2]]))
        end

        return full_path
    end

    #Configuration vector time
    function configuration_time(graph::Array{Float64,6}, configurations::Vector{Tuple{Int64, Int64, Int64}})
        total_time = 0
        for i in 1:length(configurations)-1
            next = i + 1
            n_i, v_i, h_i = configurations[i]
            n_f, v_f, h_f = configurations[next]
            total_time += graph[n_i, n_f, v_i, v_f, h_i, h_f]
        end
        return total_time
    end

    function get_score(op_params, sequence::Vector{Int64})
        depot = op_params.depots[1]
        total_score = 0
        #Depot does not have score
        i = 2
        while sequence[i] != depot
            total_score += op_params.scores[sequence[1]]
            i += 1
        end

        return total_score
    end

    #Forward graph search
    function shortest_time_by_sequence(op_params, sequence::Vector{Int64})
        num_headings = op_params.graph.num_headings
        num_speeds = op_params.graph.num_speeds
        graph = op_params.graph.graph
        depot = op_params.depots[1]
    
        #Holds best path of arbitrary position considering (speed, headingAngle)
        prev::Array{Tuple{Float64, Bool}, 2} = fill((0., false), (num_speeds, num_headings))
        curr::Array{Tuple{Float64, Bool}, 2} = fill((0., false), (num_speeds, num_headings))

        #Validates value, to remove necessity of creating new matrix everytime this runs
        valid = true
        depot_found = false
        pos = 1
        #Update best path for pos + 1, starting from second position
        while !depot_found
            curr_node = sequence[pos + 1]
            prev_node = sequence[pos]
            depot_found = curr_node == depot
            for (prev_speed, curr_speed) in itr.product(1:num_speeds, 1:num_speeds)
                for (prev_heading, curr_heading) in itr.product(1:num_headings, 1:num_headings)
                    #Not valid, assign first value for future comparisons
                    if curr[curr_speed, curr_heading][2] != valid
                        curr[curr_speed, curr_heading] = (prev[prev_speed,prev_heading][1] + graph[prev_node, curr_node, prev_speed, curr_speed, prev_heading, curr_heading], valid)
                    else
                        #Valid, compare with previously calculated value
                        val = prev[prev_speed,prev_heading][1] + graph[prev_node, curr_node, prev_speed, curr_speed, prev_heading, curr_heading]
                        if val < curr[curr_speed, curr_heading][1]
                            curr[curr_speed, curr_heading] = (val, valid)
                        end
                    end
                end
            end
            
            #Swap, swap validity if necessary
            curr, prev = prev, curr
            valid = pos % 2 == 0 ? !valid : valid
            pos += 1
        end

        return minimum(prev)[1]
    end

    function shortest_configuration_by_sequence(op_params, sequence::Vector{Int64})
        num_headings = op_params.graph.num_headings
        num_speeds = op_params.graph.num_speeds
        graph = op_params.graph.graph
        depot = op_params.depots[1]
    
        #Holds best path of arbitrary position considering (speed, headingAngle)
        prev::Array{Tuple{Float64, Bool}, 2} = fill((0., false), (num_speeds, num_headings))
        curr::Array{Tuple{Float64, Bool}, 2} = fill((0., false), (num_speeds, num_headings))

        prev_configs::Matrix{Vector{Tuple{Int64, Int64, Int64}}} = [[(sequence[1], i, j)] for i in 1:num_speeds, j in 1:num_headings]
        curr_configs = Matrix{Vector{Tuple{Int64, Int64, Int64}}}(undef, num_speeds, num_headings)

        #Validates value, to remove necessity of creating new matrix everytime this runs
        valid = true
        depot_found = false
        pos = 1
        #Update best path for pos + 1, starting from second position
        while !depot_found
            curr_node = sequence[pos + 1]
            prev_node = sequence[pos]
            depot_found = curr_node == depot
            for (prev_speed, curr_speed) in itr.product(1:num_speeds, 1:num_speeds)
                for (prev_heading, curr_heading) in itr.product(1:num_headings, 1:num_headings)
                    #Not valid, assign first value for future comparisons
                    if curr[curr_speed, curr_heading][2] != valid
                        curr[curr_speed, curr_heading] = (prev[prev_speed,prev_heading][1] + graph[prev_node, curr_node, prev_speed, curr_speed, prev_heading, curr_heading], valid)
                    
                        curr_configs[curr_speed, curr_heading] = vcat(prev_configs[prev_speed, prev_heading], [(curr_node, curr_speed, curr_heading)])
                    else
                        #Valid, compare with previously calculated value
                        val = prev[prev_speed,prev_heading][1] + graph[prev_node, curr_node, prev_speed, curr_speed, prev_heading, curr_heading]
                        if val < curr[curr_speed, curr_heading][1]
                            curr[curr_speed, curr_heading] = (val, valid)
                            curr_configs[curr_speed, curr_heading] = vcat(prev_configs[prev_speed, prev_heading], [(curr_node, curr_speed, curr_heading)])
                        end
                    end
                end
            end
            
            #Swap, swap validity if necessary
            curr, prev = prev, curr
            curr_configs, prev_configs = prev_configs, curr_configs
            valid = pos % 2 == 0 ? !valid : valid
            pos += 1
        end

        return prev_configs[findmin(prev)[2]]
    end
end