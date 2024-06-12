module Visual
    import Plots as plt
    using ..AcceleratedDubins

    function plot_full_path(positions, paths)
        p = plt.scatter([x[1] for x in positions], [x[2] for x in positions])

        x, y = [], []
        for (path, v_i, v_f) in paths
            confx, confy = AcceleratedDubins.sample_path(path)
            x = vcat(x, confx)
            y = vcat(y, confy)
        end
        plt.plot!(p, x, y)
        plt.display(p)
    end

    function plot_full_speeds(paths, parameters)
        p = plt.plot()
        params = [parameters.v_min, parameters.v_max, parameters.a_max, -parameters.a_min]

        x, y = [], []
        base_time = 0
        for (path, v_i, v_f) in paths
            times, speeds = AcceleratedDubins.speed_profile(path, params, [v_i, v_f])
            x = vcat(x, [t + base_time for t in times])
            y = vcat(y, speeds)
            
            base_time += last(times)
        end

        plt.plot!(p, x, y)
        plt.display(p)
    end
end