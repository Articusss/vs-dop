module Visual
    import PyPlot as plt
    using PyCall
    using ..AcceleratedDubins

    function filter_values(v, times)
        #Easiest example
        if length(v) == 5
            return v[2:4], times[2:4]
        elseif length(v) == 4
            return v[2:4], times[2:4] 
        end
        res_speed, res_time = [], []
        starting = v[1]
        last_time = times[1]
        i = 2
        while v[i] == starting
            i += 1
        end
        
        #Found segment start, append values
        push!(res_speed, v[i-1])
        push!(res_speed, v[i])
        push!(res_time, times[i-1])
        push!(res_time, times[i])

        # Check if next speed is different, if its not were done (peak)
        if v[i] != v[i+1]
            push!(res_speed, v[i+1])
            push!(res_time, times[i+1])
            return res_speed, res_time
        end

        #Find end
        starting = v[i]
        while v[i] == starting
            i+= 1
        end

        push!(res_speed, v[i-1])
        push!(res_speed, v[i])
        push!(res_time, times[i-1])
        push!(res_time, times[i])

        return res_speed, res_time
    end

    function sample_line(op, times, speeds, x_endpoints, y_endpoints, num_points = 1000)
        np = pyimport("numpy")
        x1, x2 = x_endpoints
        y1, y2 = y_endpoints

        t = np.linspace(0, 1, num_points)

        x_points = [x1 + x * (x2 - x1) for x in t]
        y_points = [y1 + x * (y2 - y1) for x in t]

        #Now calculate speed for each point
        starting_time = times[1]
        peak_time = times[2]
        starting_speed = speeds[1]
        peak_speed = speeds[2]

        #Discover peak distance/start
        peak_distance = starting_speed * (peak_time - starting_time) + (op.graph.vehicle_params.a_max * (peak_time - starting_time)^2) / 2
        mantains_peak = length(speeds) > 3
        peak_stop = peak_distance + peak_speed * (times[3] - peak_time)
        speed_at_conf = []
        for i in eachindex(x_points)
            x,y =  x_points[i], y_points[i]
            s = sqrt((x - x1)^2 + (y - y1)^2)
            #discover time
            if s <= peak_distance #accelerating
                v = np.sqrt(starting_speed^2 + 2 * op.graph.vehicle_params.a_max * s)
            elseif mantains_peak && s <= peak_stop #Mantaining
                v = peak_speed
            elseif mantains_peak && s > peak_stop
                v = np.sqrt(peak_speed^2 + 2 * op.graph.vehicle_params.a_min * (s - peak_stop))
            else #Deaccelerating
                v = np.sqrt(peak_speed^2 + 2 * op.graph.vehicle_params.a_min * (s - peak_distance))
            end
            push!(speed_at_conf, v)
        end

        return x_points, y_points, speed_at_conf
    end

    function plot_full_path(op, paths)
        parameters = op.graph.vehicle_params
        params = [parameters.v_min, parameters.v_max, parameters.a_max, -parameters.a_min]

        FIG_W = 14
        FIG_H = FIG_W * 9/16
        
        fig = plt.figure(figsize=(FIG_W,FIG_H), dpi=100)
        ax = fig.add_subplot(111, aspect="equal", facecolor="#EAEAF2FF")
        ax.tick_params(axis="both", which="major", labelsize=16)
        # ax.tick_params(axis='both', which='minor', labelsize=8)
        plt.grid(color="w", linestyle="solid")

        ax.set_xlabel("x [m]", fontsize=16)
        ax.set_ylabel("y [m]", fontsize=16)

        #fig, ax = plt.subplots()
        LineCollection = plt.matplotlib[:collections][:LineCollection]
        np = pyimport("numpy")

        speeds = op.graph.speeds
        vmin, vmax = minimum(speeds), maximum(speeds)
        norm = plt.matplotlib[:pyplot][:Normalize](vmin=vmin, vmax=vmax)
        #cmap = plt.get_cmap("viridis")
        #cmap_path = plt.get_cmap("YlOrRd")
        cmap_path = plt.get_cmap("autumn_r")
        cmap_locs = plt.get_cmap("cool")

        scatter = plt.scatter([x[1] for x in op.coordinates], [x[2] for x in op.coordinates], c=op.scores, cmap=cmap_locs, zorder=10, s=150, marker="s")
        plt.scatter([op.coordinates[1][1]], [op.coordinates[1][2]],  c="r", s=75)
        colorbar = plt.colorbar(scatter, pad=-0.025)
        #colorbar.set_label("Reward")
        colorbar.ax.set_title("Reward", fontsize=16)  
        colorbar.ax.tick_params(labelsize=16)
        #colorbar.set_ticks(collect(range(0, 100, 10)))   

        x, y = [], []
        for (path, v_i, v_f) in paths
            confx, confy = AcceleratedDubins.sample_path(path)
            times, speeds = AcceleratedDubins.speed_profile(path, params, [v_i, v_f])
            speeds, times = filter_values(speeds, times)
            #Plot first curve
            ax.plot(confx[1], confy[1], color=cmap_path(norm(v_i)), linewidth=3)

            #Plot straight segment
            straightx, straighty, speeds_at_conf = sample_line(op, times, speeds, confx[2], confy[2])
            points = np.reshape(np.transpose(np.array([straightx, straighty])), (-1,1,2))
            segments = np.concatenate([points[1:end-1, :, :], points[2:end, :, :]], axis=1)
            lc = LineCollection(segments, cmap=cmap_path, norm=norm, linewidth=3)
            lc.set_array(speeds_at_conf)
            ax.add_collection(lc)

            #Plot last curve
            ax.plot(confx[3], confy[3], color=cmap_path(norm(v_f)), linewidth=3)
        end 

        sm = plt.cm.ScalarMappable(cmap=cmap_path, norm=norm)
        sm.set_array(speeds)
        cbar = plt.colorbar(sm, ax=ax)
        #cbar.set_label("Speed")
        cbar.ax.set_title("Speed", fontsize=16)
        cbar.ax.tick_params(labelsize=16)
        #cbar.ax.set_clim(25, 70)
        #cbar.set_ticks(collect(range(30, 70, 8)))

        fig.savefig("path.pdf", bbox_inches="tight")
    end

    function plot_full_speeds(op, paths)
        parameters = op.graph.vehicle_params
        params = [parameters.v_min, parameters.v_max, parameters.a_max, -parameters.a_min]

        x, y = [], []
        base_time = 0
        for (path, v_i, v_f) in paths
            times, speeds = AcceleratedDubins.speed_profile(path, params, [v_i, v_f])
            x = vcat(x, [t + base_time for t in times])
            y = vcat(y, speeds)
            
            base_time += last(times)
        end

        plt.plot(x, y)
    end
end