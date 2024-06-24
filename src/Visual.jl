module Visual
    import PyPlot as plt
    using ..AcceleratedDubins

    function plot_full_path(op, paths)
        fig, ax = plt.subplots()
        speeds = op.graph.speeds
        cmap = plt.get_cmap("viridis")
        vmin, vmax = minimum(speeds), maximum(speeds)

        x, y = [], []
        for (path, v_i, v_f) in paths
            confx, confy = AcceleratedDubins.sample_path(path)
            in_color = (v_i - vmin) / (vmax - vmin)
            out_color = (v_f- vmin) / (vmax - vmin)
            
            #Plot first curve
            ax.plot(confx[1], confy[1], color=cmap(in_color))

            #Plot straight segment
            ax.plot(confx[2], confy[2], color="r")

            #Plot last curve
            ax.plot(confx[3], confy[3], color=cmap(out_color))
        end

        plt.axis("equal")

        scatter = plt.scatter([x[1] for x in op.coordinates], [x[2] for x in op.coordinates], c=op.scores, cmap="viridis", zorder=10)
        colorbar = plt.colorbar(scatter)
        colorbar.set_label("Reward")

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.matplotlib[:pyplot][:Normalize](vmin=vmin, vmax=vmax))
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