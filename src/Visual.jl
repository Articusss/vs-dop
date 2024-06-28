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

        fig, ax = plt.subplots()
        LineCollection = plt.matplotlib[:collections][:LineCollection]
        np = pyimport("numpy")

        speeds = op.graph.speeds
        vmin, vmax = minimum(speeds), maximum(speeds)
        norm = plt.matplotlib[:pyplot][:Normalize](vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("viridis")

        x, y = [], []
        for (path, v_i, v_f) in paths
            confx, confy = AcceleratedDubins.sample_path(path)
            times, speeds = AcceleratedDubins.speed_profile(path, params, [v_i, v_f])

            #Plot first curve
            confx, confy, speeds_at_conf = sample_speeds(path, speeds, times)

            #Plot straight segment
            points = np.reshape(np.transpose(np.array([confx, confy])), (-1,1,2))
            segments = np.concatenate([points[1:end-1, :, :], points[2:end, :, :]], axis=1)
            lc = LineCollection(segments, cmap=cmap, norm=norm)
            lc.set_array(speeds_at_conf)
            ax.add_collection(lc)
        end

        plt.axis("equal")

        scatter = plt.scatter([x[1] for x in op.coordinates], [x[2] for x in op.coordinates], c=op.scores, cmap="viridis", zorder=10)
        plt.scatter([op.coordinates[1][1]], [op.coordinates[1][2]], zorder=15, c="r", s=75)
        colorbar = plt.colorbar(scatter)
        colorbar.set_label("Reward")

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array(speeds)
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label("Speed")
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