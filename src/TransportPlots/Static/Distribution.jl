"""
    MomentumDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::PlotType;step=1,order=1,uDis=false,logt=false,plot_limits=(nothing,nothing),theme=DiplodocusDark())

Plots the angle averaged distribution function of of a given vector of particle `species` as a function of time given by the `sol` based on the conditions held in `PhaseSpace`. 
    
The plot can be either static, animated or interactive depending on the `type` argument. `Static` and `Animated` plots are generated using CairoMakie and are best for publications and presentations, while `Interactive` plots are generated using GLMakie and allow for user interaction with the plot.

Common arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.
- `TimeUnits`: a function that converts the time given in code units to the desired units for plotting
- `plot_limits`: the limits of the x and y axes, default is `(nothing,nothing)` which sets the limits automatically based on the data.
- `wide`: if `true`, the plot is generated in a wide format (double column 8:3 aspect ratio), default is `false` (single column 4:3 aspect ratio).
- `legend`: if `true`, a legend is added to the plot, default is `true`.
- `thermal`: default is `false`. If `true` the expected thermal distribution for each species is plotted based on the final time step of the simulation.

Static arguments:
- `step`: the step size in time to plot, default is 1.

Animated arguments:
- `framerate`: the frame rate of the animation, default is 12 fps.
- `filename`: the name of the file to save the animation to, default is "MomentumDistribution.mp4".
- `figure`: default is `nothing`, which creates a new figure. If a figure is provided, the plot is added to that figure instead of creating a new one, this should be of the form of a tuple of `figure` and `time_idx` from the main plot.
"""
function MomentumDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Static;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,thermal=false,plot_limits=(nothing,nothing),wide=false,legend=true,step=1)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    xlab = L"$\log_{10}\left(p [m_\text{Ele}c]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}]\right)$"
    elseif order != 1
        ylab = L"$\log_{10}\left(p^{%$(order)}\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_\text{Ele}c\right)^{%$(order-1)}]\right)$"
    end
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    linestyles = [:solid,(:dash,:dense),(:dot,:dense),(:dashdot,:dense),(:dashdotdot,:dense)]
    max_total = -Inf32
    p_min = Inf32
    p_max = -Inf32

    legend_elements = []
    line_labels = []

    t_save = length(sol.t)
    t_plot = ceil(Int64,t_save/step)

    values = (1:t_save)*step .+ 2 # add 2 to skip initial and kernel steps

    for (species_idx, species_name) in enumerate(species) 

    species_index = findfirst(x->x==species_name,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    dp = Grids.dpx_list[species_index]
    du = Grids.dpy_list[species_index]
    meanp = Grids.mpx_list[species_index]
    meanu = Grids.mpy_list[species_index]
    meanh = Grids.mpz_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    mass = Grids.mass_list[species_index]

    f3D = zeros(Float32,p_num,u_num,h_num)
    f2D = zeros(Float32,p_num,u_num)
    f1D = zeros(Float32,p_num)

    p_min = min(p_min,p_r[1])
    p_max = max(p_max,p_r[end])

    for i in 1:t_save

        if (i in values || i == 1 || i == 2) # plot first step for initial conds, second for kernel 

            t = sol.t[i]
            #println("t=$(CodeToSIUnitsTime(t))")
            if Time.t_grid == "u"
                color = theme.colormap[][(t - sol.t[1]) / (sol.t[end] - sol.t[1])]
            elseif Time.t_grid == "l"
                color = theme.colormap[][(log10(t) - log10(sol.t[1])) / (log10(sol.t[end]) - log10(sol.t[1]))]
            end

            f3D .= reshape(sol.f[i].x[species_index],(p_num,u_num,h_num))

            @. f3D = f3D*(f3D!=Inf)
            # scale by order
            # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
            for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
            end

            # sum along u and h directions
            pdNdp = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
            if sum(@. !isnan(pdNdp) * !isinf(pdNdp) * !iszero(pdNdp)) == 1 # there is only one valid position so scatterlines doesn't work
                idx = findfirst(!iszero,pdNdp)
                lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp[idx])],linewidth=2.0,color = color,linestyle=linestyles[species_idx])
            else
                scatterlines!(ax,log10.(meanp),log10.(pdNdp),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])
            end

            max_f = maximum(x for x in pdNdp if !isnan(x))
            max_total = max(max_f,max_total)

        end

    end

    if thermal

        # expected thermal spectrum based on final time step
        f = sol.f[end].x[species_index]
        Nᵃ = DiplodocusTransport.FourFlow(f,p_num,u_num,p_r,u_r,mass)
        Uₐ = [-1.0,0.0,0.0,0.0] # static observer
        num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
        Δab = DiplodocusTransport.ProjectionTensor(Uₐ)
        Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f,p_num,u_num,p_r,u_r,mass)
        Pressure = DiplodocusTransport.ScalarPressure(Tᵃᵇ,Δab)
        Temperature = DiplodocusTransport.ScalarTemperature(Pressure,num)

        MJ = DiplodocusTransport.MaxwellJuttner_Distribution(PhaseSpace,species[species_idx],Temperature;n=num)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        @. MJ *= (meanp^(order)) / dp

        scatterlines!(ax,log10.(meanp),log10.(MJ),linewidth=1.0,color = theme.textcolor[],markersize=0.0,label="Maxwell-Juttner")

    end

    push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
    push!(line_labels,species_name)

    end # species loop 

    t_unit_string = TimeUnits()

    if Time.t_grid == "u"
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(TimeUnits(sol.t[1]),TimeUnits(sol.t[end])),label=L"$t %$t_unit_string$")
    elseif Time.t_grid == "l"
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end]))),label=L"$\log_{10}\left(t %$t_unit_string \right)$")
    end

    if legend
        axislegend(ax,legend_elements,line_labels,position = :lt)
    end

    if plot_limits == (nothing,nothing)
        xlims!(ax,(log10(p_min)-1.0,log10(p_max)+1.0))
        ylims!(ax,(log10(max_total)-9.0,log10(max_total)+1.0)) 
    end
    #println("$((log10(p_min)-1.0,log10(p_max)+1.0))")
    #println("$((log10(max_total)-9.0,log10(max_total)+1.0))")

    return fig

    end # with_theme

end

function MomentumDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,thermal=false,plot_limits=(nothing,nothing),wide=false,legend=true,framerate=12,filename="MomentumDistribution.mp4",initial=true,figure=nothing)

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    xlab = L"$\log_{10}\left(p [m_\text{Ele}c]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}]\right)$"
    elseif order != 1
        ylab = L"$\log_{10}\left(p^{%$(order)}\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_\text{Ele}c\right)^{%$(order-1)}]\right)$"
    end

    if isnothing(figure)
        time_idx = Observable(1) # index of the current time step
        t = @lift(sol.t[$time_idx])
        fig = Figure(size =(3.25inch,3.25inch)) # 1:1 aspect ratio
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    else
        fig, time_idx = figure # use the provided figure and time index
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    end

    ax.limits = plot_limits

    line_labels = []
    legend_elements = []    

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    for (species_idx, species_name) in enumerate(species) 

    color=theme.palette.color[][mod(2*species_idx-1,7)+1]
    push!(legend_elements,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
    push!(line_labels,species_name)

    species_index = findfirst(x->x==species_name,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    dp = Grids.dpx_list[species_index]
    meanp = Grids.mpx_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    mass = Grids.mass_list[species_index]

    pdNdp = @lift begin
        f3D = zeros(Float32,p_num,u_num,h_num)
        f3D .= reshape(sol.f[$time_idx].x[species_index],(p_num,u_num,h_num))
        @. f3D = f3D*(f3D!=Inf)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        for px in 1:p_num, py in 1:u_num, pz in 1:h_num
            f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
        end
        # sum along u and h directions
        log10.(dropdims(sum(f3D, dims=(2,3)),dims=(2,3)))
    end

    scatterlines!(ax,log10.(meanp),pdNdp,linewidth=2.0,color = color,markersize=0.0,linestyle=:solid)


    if initial

        pdNdp_initial = begin
            f3D_initial = zeros(Float32,p_num,u_num,h_num)
            f3D_initial .= reshape(sol.f[1].x[species_index],(p_num,u_num,h_num))
            @. f3D_initial = f3D_initial*(f3D_initial!=Inf)
            # scale by order
            # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
            for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                f3D_initial[px,py,pz] = f3D_initial[px,py,pz] * (meanp[px]^(order)) / dp[px]
            end
            # sum along u and h directions
            log10.(dropdims(sum(f3D_initial, dims=(2,3)),dims=(2,3)))
        end

        if sum(@. !isnan(pdNdp_initial) * !isinf(pdNdp_initial) * !iszero(pdNdp_initial)) == 1 # there is only one valid position so scatterlines doesn't work
            idx = findfirst(@. !iszero(pdNdp_initial) & !isnan(pdNdp_initial) & !isinf(pdNdp_initial))
            lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, pdNdp_initial[idx]],linewidth=2.0,color = color,linestyle=(:dot))
        else
            scatterlines!(ax,log10.(meanp),log10.(pdNdp_initial),linewidth=2.0,color = color,markersize=0.0,linestyle=(:dot))
        end

    end

    if thermal

        # expected thermal spectrum based on final time step
        f = sol.f[end].x[species_index]
        Nᵃ = DiplodocusTransport.FourFlow(f,p_num,u_num,p_r,u_r,mass)
        Uₐ = [-1.0,0.0,0.0,0.0] # static observer
        num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
        Δab = DiplodocusTransport.ProjectionTensor(Uₐ)
        Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f,p_num,u_num,p_r,u_r,mass)
        Pressure = DiplodocusTransport.ScalarPressure(Tᵃᵇ,Δab)
        Temperature = DiplodocusTransport.ScalarTemperature(Pressure,num)

        MJ = DiplodocusTransport.MaxwellJuttner_Distribution(PhaseSpace,species[species_idx],Temperature;n=num)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        @. MJ *= (meanp^(order)) / dp

        scatterlines!(ax,log10.(meanp),log10.(MJ),linewidth=2.0,color = theme.textcolor[],markersize=0.0,label="Maxwell-Juttner",linestyle=(:dash,:dense))

    end

    end # species loop 

    if legend
        axislegend(ax,legend_elements,line_labels,position = :lt)
    end

    
    if !isnothing(filename)
        # recording the animation
        time_idxs = 1:length(sol.t)
        record(fig,filename,time_idxs,framerate=framerate,backend=CairoMakie) do frame
            println("$frame")
            time_idx[] = frame
        end
    end

    end # with_theme

end

"""
    MomentumAndPolarAngleDistributionPlot(sol,species::String,PhaseSpace::PhaseSpaceStruct,type::PlotType)

Plots the distribution function of of a given particle `species` as a function of momentum ``p`` and polar angle ``u`` as a function of time given by the `sol` based on the conditions held in `PhaseSpace`. 

The plot can be either static, animated or interactive depending on the `type` argument. `Static` and `Animated` plots are generated using CairoMakie and are best for publications and presentations, while `Interactive` plots are generated using GLMakie and allow for user interaction with the plot.

Common arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.

Static arguments:
- `timevalues`: a NOT OPTIONAL `tuple` of 3 either `Int64` or `Float64` time values to be plotted. `Int64` values are taken to be the time index in `sol.t`, whereas `Float64` values are taken to be actural time valus in code units that are then converted to the closest index in `sol.t`.

Animated arguments:
- `framerate`: the frame rate of the animation, default is 12 fps.
- `filename`: the name of the file to save the animation to, default is "MomentumAndPolarAngleDistribution.mp4".

"""
function MomentumAndPolarAngleDistributionPlot(sol,species::String,PhaseSpace::PhaseSpaceStruct,type::Static,timevalues::T;theme=DiplodocusDark(),order::Int64=1) where T <: Union{Tuple{Float64,Float64,Float64},Tuple{Int64,Int64,Int64}}

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    species_index = findfirst(x->x==species,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    dp = Grids.dpx_list[species_index]
    du = Grids.dpy_list[species_index]
    meanp = Grids.mpx_list[species_index]
    meanu = Grids.mpy_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    mass = Grids.mass_list[species_index]

    t_idx = zeros(Int64,3)
    t = zeros(Float64,3)

    if typeof(timevalues) == Tuple{Int64,Int64,Int64}
        for i in 1:3
        t_idx[i] = timevalues[i]
        t[i] = sol.t[t_idx[i]]
        end
    elseif typeof(timevalues) == Tuple{Float64,Float64,Float64}
        for i in 1:3
        t[i] = timevalues[i]
        t_idx[i] = findmin(abs.(sol.t .- t[i]))[2] #findfirst(x->x==t[i],sol.t)
        end
    end

    fig = Figure(size=(576,276)) # 8:3 aspect ratio

    dis1 = dropdims(sum(reshape(sol.f[t_idx[1]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    replace!(dis1,0.0 => NaN) # replace Inf with NaN for plotting
    dis2 = dropdims(sum(reshape(sol.f[t_idx[2]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    replace!(dis2,0.0 => NaN) # replace Inf with NaN for plotting
    dis3 = dropdims(sum(reshape(sol.f[t_idx[3]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    replace!(dis3,0.0 => NaN) # replace Inf with NaN for plotting

    # scale by order
    # f = dN/dpdudh * dpdudh therefore dN/dpdu = f / dpdu and p^order * dN/dpdu = f * mp^order / dpdu
    for px in 1:p_num, py in 1:u_num
        dis1[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
        dis2[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
        dis3[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
    end

    max_dis = maximum(x for x in [dis1; dis2; dis3] if !isnan(x))
    min_dis = minimum(x for x in [dis1; dis2; dis3] if !isnan(x))
    col_range = (log10(max_dis)-20.0,log10(max_dis))

    ax1 = PolarAxis(fig[1,1+1],theta_0=-pi/2,direction=-1,width=176)
    ax1.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax1,0,pi)

    ax2 = PolarAxis(fig[1,2+1],theta_0=-pi/2,direction=-1,width=176)
    ax2.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax2,0,pi)

    ax3 = PolarAxis(fig[1,3+1],theta_0=-pi/2,direction=-1,width=176)
    ax3.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax3,0,pi)

    #hm1 = heatmap!(ax1,acos.(u_r),log10.(p_r),log10.(dis1'),colormap=theme.colormap,colorrange=col_range)
    #hm2 = heatmap!(ax2,acos.(u_r),log10.(p_r),log10.(dis2'),colormap=theme.colormap,colorrange=col_range)
    #hm3 = heatmap!(ax3,acos.(u_r),log10.(p_r),log10.(dis3'),colormap=theme.colormap,colorrange=col_range)

    u_as_theta_grid = zeros(Float64,length(u_r))
    u_as_theta_grid_tick_values = Vector{Float64}(-1:0.5:1)
    u_as_theta_grid_tick_values_string = string.(u_as_theta_grid_tick_values)
    u_as_theta_grid_tick_values_string[end] = "u=1.0" 
    u_as_theta_grid_tick_locations = zeros(Float64,length(u_as_theta_grid_tick_values))
    @. u_as_theta_grid = pi - pi * (u_r+1)/2 # convert u grid to a set of theta values such that u can be plotted as polar angle
    @. u_as_theta_grid_tick_locations = pi - pi * (u_as_theta_grid_tick_values+1)/2 # convert u grid ticks to a set of theta values such that u can be plotted as polar angle

    hm1 = heatmap!(ax1,u_as_theta_grid,log10.(p_r),log10.(dis1'),colormap=theme.colormap_var,colorrange=col_range)
    hm2 = heatmap!(ax2,u_as_theta_grid,log10.(p_r),log10.(dis2'),colormap=theme.colormap_var,colorrange=col_range)
    hm3 = heatmap!(ax3,u_as_theta_grid,log10.(p_r),log10.(dis3'),colormap=theme.colormap_var,colorrange=col_range)

    rlims!(ax1,log10(p_r[1]),log10(p_r[end])+1.0)
    rlims!(ax2,log10(p_r[1]),log10(p_r[end])+1.0)
    rlims!(ax3,log10(p_r[1]),log10(p_r[end])+1.0)
    ax1.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    ax2.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    ax3.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    #hidethetadecorations!(ax1, grid=false)
    #hidethetadecorations!(ax2, grid=false)
    #hidethetadecorations!(ax3, grid=false)

    if order == 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}[\text{m}^{-3}]\right)$",flipaxis=false,height=176,tellheight=false)
    elseif order != 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p^{%$order}\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}[\text{m}^{-3}\left(m_\text{Ele}c\right)^{%$order-1}]\right)$ $$",flipaxis=false,height=176,tellheight=false)
    end

    pt = 4/3
    text!(ax1,L"$\log_{10}\left(p[m_\text{Ele}c]\right)$",position=(-3.05,log10(p_r[end])),rotation=pi/2,fontsize=9pt)
    text!(ax1,L"$t=%$(t[1])$",position=(2.8,log10(p_r[end])+4.5),fontsize=10pt)
    text!(ax2,L"$t=%$(t[2])$",position=(2.8,log10(p_r[end])+4.5),fontsize=10pt)
    text!(ax3,L"$t=%$(t[3])$",position=(2.8,log10(p_r[end])+4.5),fontsize=10pt)

    colsize!(fig.layout,1,Relative(0.1))
    colsize!(fig.layout,2,Relative(0.3))
    colsize!(fig.layout,3,Relative(0.3))
    colsize!(fig.layout,4,Relative(0.3))
    
    return fig

    end # with_theme

end

function MomentumAndPolarAngleDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,framerate=12,filename="MomentumAndPolarAngleDistribution.mp4",figure=nothing)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    if isnothing(figure)
        time_idx = Observable(1) # index of the current time step
        t = @lift(sol.t[$time_idx])
        fig = Figure(size = (3.25inch,3.25inch)) # 1:1 aspect ratio
    else
        fig, time_idx = figure # use the provided figure and time index
    end


    num_species = length(species)

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    for (species_idx, species_name) in enumerate(species)

    species_index = findfirst(x->x==species_name,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    dp = Grids.dpx_list[species_index]
    du = Grids.dpy_list[species_index]
    meanp = Grids.mpx_list[species_index]
    meanu = Grids.mpy_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    mass = Grids.mass_list[species_index]


    dis = @lift begin
        f2D = dropdims(sum(reshape(sol.f[$time_idx].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
        replace!(f2D,0.0 => NaN) # replace Inf with NaN for plotting
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dpdu = f / dpdu and p^order * dN/dpdu = f * mp^order / dpdu
        for px in 1:p_num, py in 1:u_num
            f2D[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
        end
        log10.(f2D)'
    end

    max_dis = @lift(maximum(x for x in $dis if !isnan(x)))
    #min_dis = @lift(minimum(x for x in [dis] if !isnan(x)))
    col_range = @lift(($max_dis-20.0,$max_dis))

    ax = PolarAxis(fig[1,1+species_idx],theta_0=-pi/2,direction=-1,width=Relative(1.2))
    ax.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax,0,pi)

    #hm1 = heatmap!(ax1,acos.(u_r),log10.(p_r),log10.(dis1'),colormap=theme.colormap,colorrange=col_range)
    #hm2 = heatmap!(ax2,acos.(u_r),log10.(p_r),log10.(dis2'),colormap=theme.colormap,colorrange=col_range)
    #hm3 = heatmap!(ax3,acos.(u_r),log10.(p_r),log10.(dis3'),colormap=theme.colormap,colorrange=col_range)

    u_as_theta_grid = zeros(Float64,length(u_r))
    u_as_theta_grid_tick_values = Vector{Float64}(-1:0.5:1)
    u_as_theta_grid_tick_values_string = string.(u_as_theta_grid_tick_values)
    u_as_theta_grid_tick_values_string[end] = "u=1.0" 
    u_as_theta_grid_tick_locations = zeros(Float64,length(u_as_theta_grid_tick_values))
    @. u_as_theta_grid = pi - pi * (u_r+1)/2 # convert u grid to a set of theta values such that u can be plotted as polar angle
    @. u_as_theta_grid_tick_locations = pi - pi * (u_as_theta_grid_tick_values+1)/2 # convert u grid ticks to a set of theta values such that u can be plotted as polar angle

    hm = heatmap!(ax,u_as_theta_grid,log10.(p_r),dis,colormap=theme.colormap_var,colorrange=col_range)

    rlims!(ax,log10(p_r[1]),log10(p_r[end])+1.0)
    ax.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    #hidethetadecorations!(ax1, grid=false)
    #hidethetadecorations!(ax2, grid=false)
    #hidethetadecorations!(ax3, grid=false)

    pt = 4/3
    text!(ax,L"$\log_{10}\left(p[m_\text{Ele}c]\right)$",position=(-3.00,log10(p_r[end])+1.0),rotation=pi/2,fontsize=9pt)
    if num_species != 1
        text!(ax,L"$species_name",position=(2.6,log10(p_r[end])+3.2),fontsize=10pt)
    end

    if species_idx == length(species)
    if order == 1
        Colorbar(fig[1,1],hm,label=L"$\log_{10}\left(p\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}[\text{m}^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    elseif order != 1
        Colorbar(fig[1,1],hm,label=L"$\log_{10}\left(p^{%$order}\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}[\text{m}^{-3}\left(m_\text{Ele}c\right)^{%$order-1}]\right)$ $$",flipaxis=false,height=Relative(0.75),tellheight=false)
    end
    end

    end # species loop

    colsize!(fig.layout,1,Relative(0.1))
    colsize!(fig.layout,2,Relative(0.9))
    #colgap!(fig.layout,1,Relative(-0.1))
    #colsize!(fig.layout,3,Relative(0.3))
    #colsize!(fig.layout,4,Relative(0.3))

    
    if !isnothing(filename)
        # recording the animation
        time_idxs = 1:length(sol.t)
        record(fig,filename,time_idxs,framerate=framerate,backend=CairoMakie) do frame
            println("$frame")
            time_idx[] = frame
        end
    end

    end # with_theme

end

"""
    MomentumComboAnimation(sol,species::String,PhaseSpace::PhaseSpaceStruct)

Animates the angle averaged and angle dependent particle distribution function for `species` as a function of time given by the `sol` based on the conditions held in `PhaseSpace`. 

The plot can be either static, animated or interactive depending on the `type` argument. `Static` and `Animated` plots are generated using CairoMakie and are best for publications and presentations, while `Interactive` plots are generated using GLMakie and allow for user interaction with the plot.

Arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.
- `framerate`: the frame rate of the animation, default is 12 fps.
- `filename`: the name of the file to save the animation to, default is "MomentumComboAnimation.mp4".


"""
function MomentumComboAnimation(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct;theme=DiplodocusDark(),order::Int64=1,framerate=12,filename="MomentumComboAnimation.mp4",plot_limits_momentum=(nothing,nothing))

    CairoMakie.activate!(inline=true) 

    with_theme(theme) do

        fig = Figure(size=(5.416inch,3.25inch)) # 5:3 aspect ratio

        time_idx = Observable(1) # index of the current time step
        t = @lift(sol.t[$time_idx])

        MomentumDistributionPlot(sol,species,PhaseSpace,Animated();theme=theme,order=order,TimeUnits=CodeToCodeUnitsTime,thermal=true,plot_limits=plot_limits_momentum,wide=false,legend=false,framerate=framerate,filename=nothing,initial=true,figure=(fig[2,1],time_idx))
        MomentumAndPolarAngleDistributionPlot(sol,species,PhaseSpace,Animated();order=order,theme=theme,framerate=framerate,filename=nothing,figure=(fig[1:3,2],time_idx))

        grid = fig[3,1] = GridLayout()
               
        Label(grid[1,1],@lift("t = $(round($t, digits = 1))"),fontsize=20pt)

        rowsize!(fig.layout,1,Relative(0.02))
        rowsize!(fig.layout,2,Relative(0.8))
        rowsize!(fig.layout,3,Relative(0.18))
        colsize!(fig.layout,1,Relative(0.4))
        colsize!(fig.layout,2,Relative(0.6))

        time_idxs = 1:length(sol.t)

        record(fig,filename,time_idxs,framerate=framerate,backend=CairoMakie,compression=1) do frame
            println("$frame")
            time_idx[] = frame
        end
        
    end

end




# ============== AM3 Test Plots ============== #

function AM3_MomentumDistributionPlot(filePath,t_max,t_min,t_grid;plot_limits=(nothing,nothing),theme=DiplodocusDark(),wide=false)

    fileExist = isfile(filePath)

    if fileExist
        f = DC.jldopen(filePath,"r+");

        meanp_ele = f["meanp_ele"];
        f_ele = f["f_ele"];
        t_ele = f["t_ele"];

        meanp_pho = f["meanp_pho"];
        f_pho = f["f_pho"];
        t_pho = f["t_pho"];

        DC.close(f)  
    else
        error("no file at $filePath found")
    end

    # unit conversion 
    eV_to_mElec2 = 1.60217e-19 / 9.109e-31 / 2.9979e8^2
    cm3_to_m3 = 1e6 

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    fig = Figure()
    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    xlab = L"$\log_{10}\left(p [m_\text{Ele}c]\right)$"
    ylab = L"$\log_{10}\left(p^2\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_\text{Ele}c\right)]\right)$"
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    linestyles = [:solid,:dash,:dot,:dashdot,:dashdotdot]
    legend_elements = []
    line_labels = []
    
    for i in 1:length(t_pho)

        t = t_pho[i]
        println("t=$t")
        if log10(t) % 1 == 0.0 # 10^n timesteps
            if t_grid == "u"
                color = theme.colormap[][(t - t_min) / (t_max - t_min)]
            elseif t_grid == "l"
                color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
            end

            # sum along u and h directions
            pdNdp = f_pho[i,:][1]
            #println("$pdNdp")
            scatterlines!(ax,log10.(meanp_pho .* eV_to_mElec2),log10.(meanp_pho .* pdNdp .* cm3_to_m3 .* eV_to_mElec2),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])
        end

    end

    push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
    push!(line_labels,"Pho")

    for i in 1:length(t_ele)

            t = t_ele[i]
            if log10(t) % 1 == 0.0 # 10^n timesteps
                if t_grid == "u"
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                elseif t_grid == "l"
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                end

                # sum along u and h directions
                pdNdp = f_ele[i,:][1]
                scatterlines!(ax,log10.(sqrt.((meanp_ele .* eV_to_mElec2).^2 .-1)),log10.(sqrt.((meanp_ele .* eV_to_mElec2).^2 .-1) .* pdNdp .* cm3_to_m3),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])
            end

    end

    push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[2],linewidth = 2.0))
    push!(line_labels,"Ele")

    if t_grid == "u"
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(TimeUnits(t_min),TimeUnits(t_max)),label=L"$t$ $[\text{s} * \sigma_{T}c]$")
    elseif t_grid == "l"
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(t_min),log10(t_max)),label=L"$\log_{10}\left(t [\text{s}]\right)$")
    end

    axislegend(ax,legend_elements,line_labels,position = :lt)
   
    return fig

    end # with_theme

end