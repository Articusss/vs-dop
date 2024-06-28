module Visual
    import PyPlot as plt
    using PyCall
    using ..AcceleratedDubins

    function sample_line(x_endpoints, y_endpoints, num_points = 1000)
        np = pyimport("numpy")
        x1, x2 = x_endpoints
        y1, y2 = y_endpoints

        t = np.linspace(0, 1, num_points)

        x_points = [x1 + x * (x2 - x1) for x in t]
        y_points = [y1 + x * (y2 - y1) for x in t]

        return x_points, y_points
    end

    function calculate_section_params(path, speed, times, timestamp)
        if timestamp == length(speed)
            return 0, Inf #Keep speed constant
        end

        delta_time = times[timestamp + 1] - times[timestamp]
        delta_speed = speed[timestamp + 1] - speed[timestamp]
        accel = delta_speed / delta_time
        #s = s0 + vt + 1/2at^2
        total_len = (speed[timestamp] * delta_time) + (accel * delta_time^2) / 2

        return accel, total_len
    end

    function sample_speeds(path, speed, times, resolution::Float64 = 0.05)
        points_x, points_y, speed_at_point = [], [], []
        # Sample first arc
        curr_timestamp = 1
        prev_len = 0
        curr_acc, section_len = calculate_section_params(path, speed, times, curr_timestamp)
        
        for l in 0.:resolution:sum(path.lengths)
            x, y, _ = AcceleratedDubins.get_configuration(path, l)
            push!(points_x, x)
            push!(points_y, y)
            #Calculate distance to reach a shift in acceleration
            if l - prev_len >= section_len
                curr_timestamp += 1
                prev_len += section_len
                curr_acc, section_len = calculate_section_params(path, speed, times, curr_timestamp)
            end
            #Calculate speed at point
            v = sqrt(speed[curr_timestamp]^2 + 2 * curr_acc * (l - prev_len))
            push!(speed_at_point, v)
        end
        return points_x, points_y, speed_at_point
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

        scatter = plt.scatter([x[1] for x in op.coordinates], [x[2] for x in op.coordinates], c=op.scores, cmap=cmap_locs, zorder=2, s=150, marker="s")
        plt.scatter([op.coordinates[1][1]], [op.coordinates[1][2]],  c="g", zorder=5, s=75)
        colorbar = plt.colorbar(scatter, pad=-0.025)
        #colorbar.set_label("Reward")
        colorbar.ax.set_title("Reward", fontsize=16)  
        colorbar.ax.tick_params(labelsize=16)
        #colorbar.set_ticks(collect(range(0, 100, 10)))   

        x, y = [], []
        for (path, v_i, v_f) in paths
            confx, confy = AcceleratedDubins.sample_path(path)
            times, speeds = AcceleratedDubins.speed_profile(path, params, [v_i, v_f])

            confx, confy, speeds_at_conf = sample_speeds(path, speeds, times)

            points = np.reshape(np.transpose(np.array([confx, confy])), (-1,1,2))
            segments = np.concatenate([points[1:end-1, :, :], points[2:end, :, :]], axis=1)
            lc = LineCollection(segments, cmap=cmap_path, norm=norm, linewidth=3)
            lc.set_array(speeds_at_conf)
            ax.add_collection(lc)
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