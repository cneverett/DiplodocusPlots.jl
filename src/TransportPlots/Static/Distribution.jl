"""
    MomentumDistributionPlot(sol,species,PhaseSpace;step=1,uDis=false,logt=false,plot_limits=(nothing,nothing),theme=DiplodocusDark())

Plots the angle averaged distribution function of of a given particle `species` as a function of time given by the `sol` based on the conditions held in `PhaseSpace`. 
"""
function MomentumDistributionPlot(sol,species::String,PhaseSpace::PhaseSpaceStruct;step=1,uDis=false,logt=false,plot_limits=(nothing,nothing),theme=DiplodocusDark())

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    fig = Figure()
    ax = Axis(fig[1,1],xlabel=L"$\log_{10}p$ $[m_\text{Ele}c]$",ylabel=L"$\log_{10}\left(p^2\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\right)$ $[\text{m}^{-3}]$",aspect=DataAspect())
    ax.limits = plot_limits

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    species_index = findfirst(x->x==species,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    dp = Grids.dpx_list[species_index]
    du = Grids.dpy_list[species_index]
    meanp = Grids.mpx_list[species_index]

    d = zeros(Float32,p_num,u_num)
    dj = zeros(Float32,p_num,2)

    t_save = length(sol.t)

    t_plot = ceil(Int64,t_save/step)

    color_counter = 1

    if logt
        log_t = log10.(sol.t)
        values = findall(x->x%1==0,log_t)
        my_colors = [cgrad(:rainbow)[z] for z ∈ range(0.0, 1.0, length = length(values)+1)]
    else
        values = (1:t_save)*step
        my_colors = [cgrad(:rainbow)[z] for z ∈ range(0.0, 1.0, length = t_plot+1)]
    end

    for i in 1:t_save

    if i in values || i == 1

        t = sol.t[i]

        d .= reshape(sol.f[i].x[species_index],(p_num,u_num))

        @. d = d*(d!=Inf)

        if uDis
            dj .= log10.(d .* dp)[:,4:3:8] # print just aligned and antialigned with field

            scatterlines!(ax,log10.(meanp),dj[:,1],linewidth=2.0,color = my_colors[color_counter],markersize=1.0)
            scatterlines!(ax,log10.(meanp),dj[:,2],linewidth=2.0,color = my_colors[color_counter],markersize=1.0,linestyle=:dash)
        else
            dj[:,1] .= log10.(d*du .* dp)
            scatterlines!(ax,log10.(meanp),dj[:,1],linewidth=2.0,color = my_colors[color_counter],markersize=1.0)
        end

        color_counter += 1

    end

    end

    Colorbar(fig[1,2],colormap = my_colors,limits=(log10(sol.t[2]),log10(sol.t[end])),label=L"$\log_{10}(t)$ $[\text{s} * \sigma_{T}c]$")
    
    return fig

    end # with_theme

end