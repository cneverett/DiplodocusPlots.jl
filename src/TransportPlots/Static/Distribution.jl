"""
    MomentumDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::PlotType;step=1,order=1,legend=true,thermal=false,paraperp=false,plot_limits=(nothing,nothing),TimeUnits=CodeToCodeUnitsTime,theme=DiplodocusDark())

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
- `paraperp`: default is `false`. If `true` the first and center `u` bins will be plotted to represent the distribution parallel to the axis and perpendicular.

Static arguments:
- `step`: the step size in time to plot, default is 1.

Animated arguments:
- `framerate`: the frame rate of the animation, default is 12 fps.
- `filename`: the name of the file to save the animation to, default is "MomentumDistribution.mp4".
- `figure`: default is `nothing`, which creates a new figure. If a figure is provided, the plot is added to that figure instead of creating a new one, this should be of the form of a tuple of `figure` and `time_idx` from the main plot.
- `initial`: default is `false`, if `true` causes the initial distribution to remain on the plot
"""
function MomentumDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Static;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,thermal=false,plot_limits=(nothing,nothing),wide=false,legend=true,paraperp=false,step=1)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}]\right)$"
    elseif order == 2
        ylab = L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)]\right)$"
    elseif order == -2
        ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$"
    else
        ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
    end
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]
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
    h_r = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]

    f3D = zeros(Float32,p_num,u_num,h_num)
    f1D = zeros(Float32,p_num*u_num*h_num)

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

            f1D .= copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=species_index))

            f3D .= reshape(f1D,(p_num,u_num,h_num))

            @. f3D = f3D*(f3D!=Inf)
            # scale by order
            # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
            for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
            end

            if paraperp == false
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
            elseif paraperp == true # assumes only a single particle species
                pdNdp_para = dropdims(sum(f3D, dims=(3)),dims=(3))[:,end]
                pdNdp_perp = dropdims(sum(f3D, dims=(3)),dims=(3))[:,ceil(Int64,u_num/2)]

                if sum(@. !isnan(pdNdp_para) * !isinf(pdNdp_para) * !iszero(pdNdp_para)) == 1 # there is only one valid position so scatterlines doesn't work
                    idx = findfirst(!iszero,pdNdp_para)
                    lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp_para[idx])],linewidth=2.0,color = color,linestyle=linestyles[1])
                else
                    scatterlines!(ax,log10.(meanp),log10.(pdNdp_para),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])
                end

                max_f = maximum(x for x in pdNdp_para if !isnan(x))
                max_total = max(max_f,max_total)

                if sum(@. !isnan(pdNdp_perp) * !isinf(pdNdp_perp) * !iszero(pdNdp_perp)) == 1 # there is only one valid position so scatterlines doesn't work
                    idx = findfirst(!iszero,pdNdp_perp)
                    lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp_perp[idx])],linewidth=2.0,color = color,linestyle=linestyles[2])
                else
                    scatterlines!(ax,log10.(meanp),log10.(pdNdp_perp),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])
                end

                max_f = maximum(x for x in pdNdp_perp if !isnan(x))
                max_total = max(max_f,max_total)

            end

            

        end

    end

    if thermal

        # expected thermal spectrum based on final time step
        f = copy(Location_Species_To_StateVector(sol.f[end],PhaseSpace,species_index=species_index))
        Nᵃ = DiplodocusTransport.FourFlow(f,p_num,u_num,h_num,p_r,u_r,h_r,mass)
        Uₐ = [-1.0,0.0,0.0,0.0] # static observer
        num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
        Δab = DiplodocusTransport.ProjectionTensor(Uₐ)
        Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f,p_num,u_num,h_num,p_r,u_r,h_r,mass)
        Pressure = DiplodocusTransport.ScalarPressure(Tᵃᵇ,Δab)
        Temperature = DiplodocusTransport.ScalarTemperature(Pressure,num)

        MJ = DiplodocusTransport.MaxwellJuttner_Distribution(PhaseSpace,species[species_idx],Temperature;n=num)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        @. MJ *= (meanp^(order)) / dp

        scatterlines!(ax,log10.(meanp),log10.(MJ),linewidth=1.0,color = theme.textcolor[],markersize=0.0,label="Maxwell-Juttner")

    end

    if paraperp == false
        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
        push!(line_labels,species_name)
    end

    end # species loop 

    if paraperp == true
        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[2],linewidth = 2.0))
        push!(line_labels,L"\parallel")
        push!(line_labels,L"\perp")
    end

    t_unit_string = TimeUnits()

    if Time.t_grid == "u"
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(TimeUnits(sol.t[1]),TimeUnits(sol.t[end])),label=L"$t\,$ $%$t_unit_string$")
    elseif Time.t_grid == "l"
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(round(TimeUnits(sol.t[1]),sigdigits=5)),log10(round(TimeUnits(sol.t[end]),sigdigits=5))),label=L"$\log_{10}\left(t\,%$t_unit_string \right)$")
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

function MomentumDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,thermal=false,plot_limits=(nothing,nothing),wide=false,legend=true,framerate=12,filename="MomentumDistribution.mp4",initial=true,paraperp=false,figure=nothing)

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}]\right)$"
    elseif order == 2
        ylab = L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)]\right)$"
    elseif order == -2
        ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$"
    else 
        ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
    end

    if isnothing(figure)
        time_idx = Observable(1) # index of the current time step
        t = @lift(sol.t[$time_idx])
        if wide
            fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
        else
            fig = Figure() # default single column 4:3 aspect ratio
        end
        #fig = Figure(size =(3.25inch,3.25inch)) # 1:1 aspect ratio
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    else
        fig, time_idx = figure # use the provided figure and time index
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    end

    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot,:dense),(:dashdotdot,:dense)]
    line_labels = []
    legend_elements = []    

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    for (species_idx, species_name) in enumerate(species) 

    color=theme.palette.color[][mod(2*species_idx-1,7)+1]
    if paraperp == false
        push!(legend_elements,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
        push!(line_labels,species_name)
    end

    species_index = findfirst(x->x==species_name,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    dp = Grids.dpx_list[species_index]
    meanp = Grids.mpx_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    h_r = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]

    if paraperp == false 
        pdNdp = @lift begin
        f1D = zeros(Float32,p_num*u_num*h_num)
        f1D .= copy(Location_Species_To_StateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index))
        f3D = zeros(Float32,p_num,u_num,h_num)
        f3D .= reshape(f1D,(p_num,u_num,h_num))
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

    elseif paraperp == true

        pdNdp_para = @lift begin
        f1D = zeros(Float32,p_num*u_num*h_num)
        f1D .= copy(Location_Species_To_StateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index))
        f3D = zeros(Float32,p_num,u_num,h_num)
        f3D .= reshape(f1D,(p_num,u_num,h_num))
        @. f3D = f3D*(f3D!=Inf)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        for px in 1:p_num, py in 1:u_num, pz in 1:h_num
            f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
        end
        # sum along u and h directions
        log10.(dropdims(sum(f3D, dims=(3)),dims=(3)))[:,1]
        end

        pdNdp_perp = @lift begin
        f1D = zeros(Float32,p_num*u_num*h_num)
        f1D .= copy(Location_Species_To_StateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index))
        f3D = zeros(Float32,p_num,u_num,h_num)
        f3D .= reshape(f1D,(p_num,u_num,h_num))
        @. f3D = f3D*(f3D!=Inf)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        for px in 1:p_num, py in 1:u_num, pz in 1:h_num
            f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
        end
        # sum along u and h directions
        log10.(dropdims(sum(f3D, dims=(3)),dims=(3)))[:,round(Int64,u_num/2)]
        end

        scatterlines!(ax,log10.(meanp),pdNdp_para,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])
        scatterlines!(ax,log10.(meanp),pdNdp_perp,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])

    end

    if initial

        pdNdp_initial = begin
            f1D_initial = zeros(Float32,p_num*u_num*h_num)
            f1D_initial .= copy(Location_Species_To_StateVector(sol.f[1],PhaseSpace,species_index=species_index))
            f3D_initial = zeros(Float32,p_num,u_num,h_num)
            f3D_initial .= reshape(f1D_initial,(p_num,u_num,h_num))
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
            scatterlines!(ax,log10.(meanp),pdNdp_initial,linewidth=2.0,color = color,markersize=0.0,linestyle=(:dot))
        end

    end

    if thermal

        # expected thermal spectrum based on final time step
        f = copy(Location_Species_To_StateVector(sol.f[end],PhaseSpace,species_index=species_index))
        Nᵃ = DiplodocusTransport.FourFlow(f,p_num,u_num,h_num,p_r,u_r,h_r,mass)
        Uₐ = [-1.0,0.0,0.0,0.0] # static observer
        num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
        Δab = DiplodocusTransport.ProjectionTensor(Uₐ)
        Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f,p_num,u_num,h_num,p_r,u_r,h_r,mass)
        Pressure = DiplodocusTransport.ScalarPressure(Tᵃᵇ,Δab)
        Temperature = DiplodocusTransport.ScalarTemperature(Pressure,num)

        MJ = DiplodocusTransport.MaxwellJuttner_Distribution(PhaseSpace,species[species_idx],Temperature;n=num)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        @. MJ *= (meanp^(order)) / dp

        scatterlines!(ax,log10.(meanp),log10.(MJ),linewidth=2.0,color = theme.textcolor[],markersize=0.0,label="Maxwell-Juttner",linestyle=(:dash,:dense))

    end

    end # species loop 

    if paraperp == true
        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[2],linewidth = 2.0))
        push!(line_labels,L"\parallel")
        push!(line_labels,L"\perp")
    end

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
    AngleDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::PlotType;angle_step=1,order=1,plot_limits=(nothing,nothing),theme=DiplodocusDark())

Plots the angle averaged distribution function of of a given vector of particle `species` as a function of time given by the `sol` based on the conditions held in `PhaseSpace`. 
    
The plot can be either static, animated or interactive depending on the `type` argument. `Static` and `Animated` plots are generated using CairoMakie and are best for publications and presentations, while `Interactive` plots are generated using GLMakie and allow for user interaction with the plot.

Common arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.
- `TimeUnits`: a function that converts the time given in code units to the desired units for plotting
- `plot_limits`: the limits of the x and y axes, default is `(nothing,nothing)` which sets the limits automatically based on the data.
- `wide`: if `true`, the plot is generated in a wide format (double column 8:3 aspect ratio), default is `false` (single column 4:3 aspect ratio).
- `legend`: if `true`, a legend is added to the plot, default is `true`.
- `angle_step`: the angular grid step to reduce plotting lines.

Static arguments:
- `time_idx`: the time index to plot at.

Animated arguments:
- `framerate`: the frame rate of the animation, default is 12 fps.
- `filename`: the name of the file to save the animation to, default is "AngleDistribution.mp4".
- `figure`: default is `nothing`, which creates a new figure. If a figure is provided, the plot is added to that figure instead of creating a new one, this should be of the form of a tuple of `figure` and `time_idx` from the main plot.
"""
function AngleDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Static,time_idx::Int64;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,plot_limits=(nothing,nothing),wide=false,legend=true,angle_step::Int64=1)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    xlab = L"$\log_{10}\left(p [m_ec]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}]\right)$"
    elseif order != 1
        ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
    end
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot,:dense),(:dashdotdot,:dense)]
    max_total = -Inf32
    p_min = Inf32
    p_max = -Inf32

    legend_elements = []
    line_labels = []
    legend_elements_angle = []
    line_labels_angle = []

    counter = 1

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
    h_r = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]

    f3D = zeros(Float32,p_num,u_num,h_num)
    f1D = zeros(Float32,p_num*u_num*h_num)

    p_min = min(p_min,p_r[1])
    p_max = max(p_max,p_r[end])

    t = sol.t[time_idx]

    f1D .= copy(Location_Species_To_StateVector(sol.f[time_idx],PhaseSpace,species_index=species_index))

    f3D .= reshape(f1D,(p_num,u_num,h_num))
    @. f3D = f3D*(f3D!=Inf)
    # scale by order
    # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
    for px in 1:p_num, py in 1:u_num, pz in 1:h_num
        f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
    end

    for i in ceil(Int64,u_num/2):angle_step:u_num

        u_val = meanu[i]

        color = theme.colormap[][(u_val - u_r[ceil(Int64,u_num/2)]) / (u_r[end] - u_r[ceil(Int64,u_num/2)])]

        pdNdp= dropdims(sum(f3D, dims=(3)),dims=(3))[:,i]

        if sum(@. !isnan(pdNdp) * !isinf(pdNdp) * !iszero(pdNdp)) == 1 # there is only one valid position so scatterlines doesn't work
            idx = findfirst(!iszero,pdNdp)
            lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp[idx])],linewidth=2.0,color = color,linestyle=linestyles[species_idx])
        else
            scatterlines!(ax,log10.(meanp),log10.(pdNdp),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])
        end
        max_f = maximum(x for x in pdNdp if !isnan(x))
        max_total = max(max_f,max_total)

        if legend==true && counter == 1 
        u_val = meanu[i]
        color = theme.colormap[][(u_val - u_r[ceil(Int64,u_num/2)]) / (u_r[end] - u_r[ceil(Int64,u_num/2)])]
        push!(legend_elements_angle,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
        push!(line_labels_angle,L"$u=%$(round(u_val,sigdigits=2))$")
        end

    end # angle loop

        counter += 1

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
        push!(line_labels,species_name)

    end # species loop 


    if legend
        axislegend(ax,legend_elements,line_labels,position = :lt)
        axislegend(ax,legend_elements_angle,line_labels_angle,position = :lb)
    end

    if plot_limits == (nothing,nothing)
        xlims!(ax,(log10(p_min)-1.0,log10(p_max)+1.0))
        ylims!(ax,(log10(max_total)-9.0,log10(max_total)+1.0)) 
    end

    return fig

    end # with_theme

end

function AngleDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,plot_limits=(nothing,nothing),wide=false,legend=true,angle_step::Int64=1,framerate=12,filename="AngleDistribution.mp4",figure=nothing)

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    xlab = L"$\log_{10}\left(p [m_ec]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}]\right)$"
    elseif order == -2
        ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$"
    elseif order != 1
        ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
    end

    if isnothing(figure)
        time_idx = Observable(1) # index of the current time step
        t = @lift(TimeUnits(sol.t[$time_idx]))
        if wide
            fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
        else
            fig = Figure() # default single column 4:3 aspect ratio
        end
        #fig = Figure(size =(3.25inch,3.25inch)) # 1:1 aspect ratio
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    else
        fig, time_idx = figure # use the provided figure and time index
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    end

    ax.limits = plot_limits

    t_unit_string = TimeUnits()
               
    text!(ax,@lift("t=$(round($(t), sigdigits = 3))  "),fontsize=12pt,align=(:right,:center),space=:relative,offset=(375.0,135.0))
    text!(ax,L"%$t_unit_string",fontsize=12pt,align=(:left,:center),space=:relative,offset=(375.0,135.0))

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot,:dense),(:dashdotdot,:dense)]
    line_labels = []
    legend_elements = []  
    legend_elements_angle = []
    line_labels_angle = []  

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    counter = 1;

    for (species_idx, species_name) in enumerate(species) 

    species_index = findfirst(x->x==species_name,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    dp = Grids.dpx_list[species_index]
    meanp = Grids.mpx_list[species_index]
    meanu = Grids.mpy_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    mass = Grids.mass_list[species_index]

    for i in ceil(Int64,u_num/2):angle_step:u_num

        u_val = meanu[i]
        color = theme.colormap[][(u_val - u_r[ceil(Int64,u_num/2)]) / (u_r[end] - u_r[ceil(Int64,u_num/2)])]

        pdNdp = @lift begin
        f1D = zeros(Float32,p_num*u_num*h_num)
        f1D .= copy(Location_Species_To_StateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index))
        f3D = zeros(Float32,p_num,u_num,h_num)
        f3D .= reshape(f1D,(p_num,u_num,h_num))
        @. f3D = f3D*(f3D!=Inf)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        for px in 1:p_num, py in 1:u_num, pz in 1:h_num
            f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
        end
        # sum along u and h directions
        log10.(dropdims(sum(f3D, dims=(3)),dims=(3)))[:,i]
        end

        scatterlines!(ax,log10.(meanp),pdNdp,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])

        if legend==true && counter == 1 
        u_val = meanu[i]
        color = theme.colormap[][(u_val - u_r[ceil(Int64,u_num/2)]) / (u_r[end] - u_r[ceil(Int64,u_num/2)])]
        push!(legend_elements_angle,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
        push!(line_labels_angle,L"$\theta=%$(round(acos(u_val)/pi,sigdigits=2)) \pi$")
        end

    end # angle loop 

        counter += 1

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
        push!(line_labels,species_name)

    end # species loop 

    if legend
        axislegend(ax,legend_elements,line_labels,position = :lt)
        axislegend(ax,legend_elements_angle,line_labels_angle,position = :lb)
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
function MomentumAndPolarAngleDistributionPlot(sol,species::String,PhaseSpace::PhaseSpaceStruct,type::Static,timevalues::T;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime) where T <: Union{Tuple{Float64,Float64,Float64},Tuple{Int64,Int64,Int64}}

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

    f1 = copy(Location_Species_To_StateVector(sol.f[t_idx[1]],PhaseSpace,species_index=species_index))
    f2 = copy(Location_Species_To_StateVector(sol.f[t_idx[2]],PhaseSpace,species_index=species_index))
    f3 = copy(Location_Species_To_StateVector(sol.f[t_idx[3]],PhaseSpace,species_index=species_index))

    #dis1 = dropdims(sum(reshape(sol.f[t_idx[1]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    #dis2 = dropdims(sum(reshape(sol.f[t_idx[2]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    #dis3 = dropdims(sum(reshape(sol.f[t_idx[3]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)

    dis1 = dropdims(sum(reshape(f1,(p_num,u_num,h_num)),dims=3),dims=3)
    dis2 = dropdims(sum(reshape(f2,(p_num,u_num,h_num)),dims=3),dims=3)
    dis3 = dropdims(sum(reshape(f3,(p_num,u_num,h_num)),dims=3),dims=3)

    # scale by order
    # f = dN/dpdudh * dpdudh therefore dN/dpdu = f / dpdu and p^order * dN/dpdu = f * mp^order / dpdu
    for px in 1:p_num, py in 1:u_num
        dis1[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
        dis2[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
        dis3[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
    end
    replace!(dis1,0.0 => NaN) # replace Inf with NaN for plotting
    replace!(dis2,0.0 => NaN) # replace Inf with NaN for plotting
    replace!(dis3,0.0 => NaN) # replace Inf with NaN for plotting

    max_dis = maximum(x for x in [dis1; dis2; dis3] if !isnan(x))
    min_dis = minimum(x for x in [dis1; dis2; dis3] if !isnan(x))
    col_range = (log10(max_dis)-16.0,log10(max_dis))

    ax1 = PolarAxis(fig[1,1+1],theta_0=-pi/2,direction=-1,width=176)
    ax1.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax1,0,pi)

    ax2 = PolarAxis(fig[1,2+1],theta_0=-pi/2,direction=-1,width=176)
    ax2.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax2,0,pi)

    ax3 = PolarAxis(fig[1,3+1],theta_0=-pi/2,direction=-1,width=176)
    ax3.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax3,0,pi)

    u_as_theta_grid = zeros(Float64,length(u_r))
    u_as_theta_grid_tick_values = Vector{Float64}(-1:0.5:1)
    u_as_theta_grid_tick_values_string = string.(u_as_theta_grid_tick_values)
    u_as_theta_grid_tick_values_string[end] = "u=1.0" 
    u_as_theta_grid_tick_locations = zeros(Float64,length(u_as_theta_grid_tick_values))
    @. u_as_theta_grid = pi - pi * (u_r+1)/2 # convert u grid to a set of theta values such that u can be plotted as polar angle
    @. u_as_theta_grid_tick_locations = pi - pi * (u_as_theta_grid_tick_values+1)/2 # convert u grid ticks to a set of theta values such that u can be plotted as polar angle

    hm1 = heatmap!(ax1,u_as_theta_grid,log10.(p_r),log10.(dis1'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
    hm2 = heatmap!(ax2,u_as_theta_grid,log10.(p_r),log10.(dis2'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
    hm3 = heatmap!(ax3,u_as_theta_grid,log10.(p_r),log10.(dis3'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))

    translate!(hm1,0,0,-100)
    translate!(hm2,0,0,-100)
    translate!(hm3,0,0,-100)

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
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    elseif order == -2
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    elseif order != 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    end

    t_unit_string = TimeUnits()

    pt = 4/3
    text!(ax1,L"$\log_{10}\left(p\,[m_ec]\right)$",position=(-3.05,log10(p_r[end])),rotation=pi/2,fontsize=9pt)
    ax1_label=Axis(fig[2,2])
    ax2_label=Axis(fig[2,3])
    ax3_label=Axis(fig[2,4])
    hidedecorations!(ax1_label)
    hidedecorations!(ax2_label)
    hidedecorations!(ax3_label)
    hidespines!(ax1_label)
    hidespines!(ax2_label)
    hidespines!(ax3_label)
    text!(ax1_label,L"$t=%$(round(TimeUnits(t[1]),sigdigits=3))$ $%$t_unit_string$",space=:relative,position=(0.5,0.5),fontsize=10pt,align=(:center,:center))
    text!(ax2_label,L"$t=%$(round(TimeUnits(t[2]),sigdigits=3))$ $%$t_unit_string$",space=:relative,position=(0.5,0.5),fontsize=10pt,align=(:center,:center))
    text!(ax3_label,L"$t=%$(round(TimeUnits(t[3]),sigdigits=3))$ $%$t_unit_string$",space=:relative,position=(0.5,0.5),fontsize=10pt,align=(:center,:center))

    colsize!(fig.layout,1,Relative(0.1))
    colsize!(fig.layout,2,Relative(0.3))
    colsize!(fig.layout,3,Relative(0.3))
    colsize!(fig.layout,4,Relative(0.3))
    rowsize!(fig.layout,2,Relative(0.06))
    rowgap!(fig.layout,1,0.0)
    
    return fig

    end # with_theme

end

function MomentumAndPolarAngleDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,framerate=12,filename="MomentumAndPolarAngleDistribution.mp4",figure=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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
        f1D = copy(Location_Species_To_StateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index))
        f2D = dropdims(sum(reshape(f1D,(p_num,u_num,h_num)),dims=3),dims=3)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dpdu = f / dpdu and p^order * dN/dpdu = f * mp^order / dpdu
        for px in 1:p_num, py in 1:u_num
            f2D[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
        end
        replace!(f2D,0.0 => NaN) # replace Inf with NaN for plotting
        log10.(f2D)'
    end

    max_dis = @lift(maximum(x for x in $dis if !isnan(x)))
    #min_dis = @lift(minimum(x for x in [dis] if !isnan(x)))
    col_range = @lift(($max_dis-16.0,$max_dis))

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

    hm = heatmap!(ax,u_as_theta_grid,log10.(p_r),dis,colormap=theme.colormap_var,colorrange=col_range,colorscale=@lift(x->asinh(x-$max_dis)))

    translate!(hm,0,0,-100)

    rlims!(ax,log10(p_r[1]),log10(p_r[end])+1.0)
    ax.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    #hidethetadecorations!(ax1, grid=false)
    #hidethetadecorations!(ax2, grid=false)
    #hidethetadecorations!(ax3, grid=false)

    pt = 4/3
    text!(ax,L"$\log_{10}\left(p\,[m_ec]\right)$",position=(-3.00,log10(p_r[end])+1.0),rotation=pi/2,fontsize=9pt)
    if num_species != 1
        text!(ax,L"$species_name",position=(2.6,log10(p_r[end])+3.2),fontsize=10pt)
    end

    if species_idx == length(species)
    if order == 1
        Colorbar(fig[1,1],hm,label=L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    elseif order == -2
        Colorbar(fig[1,1],hm,label=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    elseif order != 1
        Colorbar(fig[1,1],hm,label=L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
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
    MomentumAndAzimuthalAngleDistributionPlot(sol,species::String,PhaseSpace::PhaseSpaceStruct,type::PlotType)

Plots the distribution function of of a given particle `species` as a function of momentum ``p`` and azimuthal angle ``h`` as a function of time given by the `sol` based on the conditions held in `PhaseSpace`. 

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
function MomentumAndAzimuthalAngleDistributionPlot(sol,species::String,PhaseSpace::PhaseSpaceStruct,type::Static,timevalues::T;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime) where T <: Union{Tuple{Float64,Float64,Float64},Tuple{Int64,Int64,Int64}}

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
    dh = Grids.dpz_list[species_index]
    meanp = Grids.mpx_list[species_index]
    meanu = Grids.mpy_list[species_index]
    meanh = Grids.mpz_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    h_r = Grids.pzr_list[species_index]
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

    f1 = copy(Location_Species_To_StateVector(sol.f[t_idx[1]],PhaseSpace,species_index=species_index))
    f2 = copy(Location_Species_To_StateVector(sol.f[t_idx[2]],PhaseSpace,species_index=species_index))
    f3 = copy(Location_Species_To_StateVector(sol.f[t_idx[3]],PhaseSpace,species_index=species_index))

    #dis1 = dropdims(sum(reshape(sol.f[t_idx[1]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    #dis2 = dropdims(sum(reshape(sol.f[t_idx[2]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    #dis3 = dropdims(sum(reshape(sol.f[t_idx[3]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)

    # sum over u angle
    dis1 = dropdims(sum(reshape(f1,(p_num,u_num,h_num)),dims=2),dims=2)
    dis2 = dropdims(sum(reshape(f2,(p_num,u_num,h_num)),dims=2),dims=2)
    dis3 = dropdims(sum(reshape(f3,(p_num,u_num,h_num)),dims=2),dims=2)

    # scale by order
    # f = dN/dpdudh * dpdudh therefore dN/dpdu = f / dpdu and p^order * dN/dpdu = f * mp^order / dpdu
    for px in 1:p_num, pz in 1:h_num
        dis1[px,pz] *= (meanp[px]^(order)) / dp[px] / dh[pz]
        dis2[px,pz] *= (meanp[px]^(order)) / dp[px] / dh[pz]
        dis3[px,pz] *= (meanp[px]^(order)) / dp[px] / dh[pz]
    end
    replace!(dis1,0.0 => NaN) # replace Inf with NaN for plotting
    replace!(dis2,0.0 => NaN) # replace Inf with NaN for plotting
    replace!(dis3,0.0 => NaN) # replace Inf with NaN for plotting

    max_dis = maximum(x for x in [dis1; dis2; dis3] if !isnan(x))
    min_dis = minimum(x for x in [dis1; dis2; dis3] if !isnan(x))
    col_range = (log10(max_dis)-24.0,log10(max_dis))

    ax1 = PolarAxis(fig[1,1+1],width=176)
    thetalims!(ax1,0,2pi)

    ax2 = PolarAxis(fig[1,2+1],width=176)
    thetalims!(ax2,0,2pi)

    ax3 = PolarAxis(fig[1,3+1],width=176)
    thetalims!(ax3,0,2pi)

    #u_as_theta_grid = zeros(Float64,length(u_r))
    #u_as_theta_grid_tick_values = Vector{Float64}(-1:0.5:1)
    #u_as_theta_grid_tick_values_string = string.(u_as_theta_grid_tick_values)
    #u_as_theta_grid_tick_values_string[end] = "u=1.0" 
    #u_as_theta_grid_tick_locations = zeros(Float64,length(u_as_theta_grid_tick_values))
    #@. u_as_theta_grid = pi - pi * (u_r+1)/2 # convert u grid to a set of theta values such that u can be plotted as polar angle
    #@. u_as_theta_grid_tick_locations = pi - pi * (u_as_theta_grid_tick_values+1)/2 # convert u grid ticks to a set of theta values such that u can be plotted as polar angle

    if PhaseSpace.Momentum.px_grid_list[species_index] == "l"
        hm1 = heatmap!(ax1,h_r,log10.(p_r),log10.(dis1'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
        hm2 = heatmap!(ax2,h_r,log10.(p_r),log10.(dis2'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
        hm3 = heatmap!(ax3,h_r,log10.(p_r),log10.(dis3'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
        rlims!(ax1,log10(p_r[1]),log10(p_r[end])+1.0)
        rlims!(ax2,log10(p_r[1]),log10(p_r[end])+1.0)
        rlims!(ax3,log10(p_r[1]),log10(p_r[end])+1.0)
        ax1.radius_at_origin = log10(p_r[1])-1.0
        ax2.radius_at_origin = log10(p_r[1])-1.0
        ax3.radius_at_origin = log10(p_r[1])-1.0
    elseif PhaseSpace.Momentum.px_grid_list[species_index] == "u"
        hm1 = heatmap!(ax1,h_r,p_r,log10.(dis1'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
        hm2 = heatmap!(ax2,h_r,p_r,log10.(dis2'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
        hm3 = heatmap!(ax3,h_r,p_r,log10.(dis3'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
        rlims!(ax1,0.0,p_r[end])
        rlims!(ax2,0.0,p_r[end])
        rlims!(ax3,0.0,p_r[end])
    end
    #ax1.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    #ax2.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    #ax3.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    #hidethetadecorations!(ax1, grid=false)
    #hidethetadecorations!(ax2, grid=false)
    #hidethetadecorations!(ax3, grid=false)
    ax1.thetagridcolor=(:grey45,0.5)
    ax1.rgridcolor=(:grey45,0.5)
    ax2.thetagridcolor=(:grey45,0.5)
    ax2.rgridcolor=(:grey45,0.5)
    ax3.thetagridcolor=(:grey45,0.5)
    ax3.rgridcolor=(:grey45,0.5)


    if order == 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}\phi\mathrm{d}V}[\text{m}^{-3}]\right)$",flipaxis=false,height=176,tellheight=false)
    elseif order != 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}\phi\mathrm{d}V}[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$",flipaxis=false,height=176,tellheight=false)
    end

    t_unit_string = TimeUnits()

    pt = 4/3
    #text!(ax1,L"$\log_{10}\left(p[m_ec]\right)$",position=(-3.05,log10(p_r[end])),rotation=pi/2,fontsize=9pt)
    #text!(ax1,L"$t=%$(round(TimeUnits(t[1]),sigdigits=3))$ $%$t_unit_string$",position=(2.8,log10(p_r[end])+4.5),fontsize=10pt)
    #text!(ax2,L"$t=%$(round(TimeUnits(t[2]),sigdigits=3))$ $%$t_unit_string$",position=(2.8,log10(p_r[end])+4.5),fontsize=10pt)
    #text!(ax3,L"$t=%$(round(TimeUnits(t[3]),sigdigits=3))$ $%$t_unit_string$",position=(2.8,log10(p_r[end])+4.5),fontsize=10pt)

    colsize!(fig.layout,1,Relative(0.1))
    colsize!(fig.layout,2,Relative(0.3))
    colsize!(fig.layout,3,Relative(0.3))
    colsize!(fig.layout,4,Relative(0.3))
    
    return fig

    end # with_theme

end

function MomentumAndAzimuthalAngleDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,framerate=12,filename="MomentumAndPolarAngleDistribution.mp4",figure=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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
    dh = Grids.dpz_list[species_index]
    meanp = Grids.mpx_list[species_index]
    meanu = Grids.mpy_list[species_index]
    meanh = Grids.mpz_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    h_r = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]


    dis = @lift begin
        f1D = copy(Location_Species_To_StateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index))
        # sum over u angles
        f2D = dropdims(sum(reshape(f1D,(p_num,u_num,h_num)),dims=2),dims=2)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dpdu = f / dpdu and p^order * dN/dpdu = f * mp^order / dpdu
        for px in 1:p_num, pz in 1:h_num
            f2D[px,pz] *= (meanp[px]^(order)) / dp[px] / dh[pz]
        end
        replace!(f2D,0.0 => NaN) # replace Inf with NaN for plotting
        log10.(f2D)'
    end

    max_dis = @lift(maximum(x for x in $dis if !isnan(x)))
    #min_dis = @lift(minimum(x for x in [dis] if !isnan(x)))
    col_range = @lift(($max_dis-24.0,$max_dis))

    ax = PolarAxis(fig[1,1+species_idx],width=Relative(1.2))
    ax.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax,0,2pi)

    #hm1 = heatmap!(ax1,acos.(u_r),log10.(p_r),log10.(dis1'),colormap=theme.colormap,colorrange=col_range)
    #hm2 = heatmap!(ax2,acos.(u_r),log10.(p_r),log10.(dis2'),colormap=theme.colormap,colorrange=col_range)
    #hm3 = heatmap!(ax3,acos.(u_r),log10.(p_r),log10.(dis3'),colormap=theme.colormap,colorrange=col_range)

    #u_as_theta_grid = zeros(Float64,length(u_r))
    #u_as_theta_grid_tick_values = Vector{Float64}(-1:0.5:1)
    #u_as_theta_grid_tick_values_string = string.(u_as_theta_grid_tick_values)
    #u_as_theta_grid_tick_values_string[end] = "u=1.0" 
    #u_as_theta_grid_tick_locations = zeros(Float64,length(u_as_theta_grid_tick_values))
    #@. u_as_theta_grid = pi - pi * (u_r+1)/2 # convert u grid to a set of theta values such that u can be plotted as polar angle
    #@. u_as_theta_grid_tick_locations = pi - pi * (u_as_theta_grid_tick_values+1)/2 # convert u grid ticks to a set of theta values such that u can be plotted as polar angle

    hm = heatmap!(ax,h_r,log10.(p_r),dis,colormap=theme.colormap_var,colorrange=col_range,colorscale=@lift(x->asinh(x-$max_dis)))

    rlims!(ax,log10(p_r[1]),log10(p_r[end])+1.0)
    #ax.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    #hidethetadecorations!(ax1, grid=false)
    #hidethetadecorations!(ax2, grid=false)
    #hidethetadecorations!(ax3, grid=false)

    pt = 4/3
    text!(ax,L"$\log_{10}\left(p[m_ec]\right)$",position=(-3.00,log10(p_r[end])+1.0),rotation=pi/2,fontsize=9pt)
    if num_species != 1
        text!(ax,L"$species_name",position=(2.6,log10(p_r[end])+3.2),fontsize=10pt)
    end

    if species_idx == length(species)
    if order == 1
        Colorbar(fig[1,1],hm,label=L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}\phi\mathrm{d}V}[\text{m}^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    elseif order != 1
        Colorbar(fig[1,1],hm,label=L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}\phi\mathrm{d}V}[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
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
- `plot_limits_momentum`: the limits of the momentum plot, default is `(nothing,nothing)`.
- `thermal`: whether to plot the expected thermal spectrum based on the final time step, default is `false`.
- `paraperp`: default is `false`. If `true` the first and center `u` bins will be plotted to represent the distribution parallel to the axis and perpendicular.
- `initial`: default is `false`, if `true` causes the initial distribution to remain on the plot


"""
function MomentumComboAnimation(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct;theme=DiplodocusDark(),order::Int64=1,thermal=false,paraperp=false,legend=false,initial=false,framerate=12,filename="MomentumComboAnimation.mp4",plot_limits_momentum=(nothing,nothing),TimeUnits::Function=CodeToCodeUnitsTime)

    CairoMakie.activate!(inline=true) 

    with_theme(theme) do

        fig = Figure(size=(5.416inch,3.25inch)) # 5:3 aspect ratio

        time_idx = Observable(1) # index of the current time step
        t = @lift(TimeUnits(sol.t[$time_idx]))

        MomentumDistributionPlot(sol,species,PhaseSpace,Animated();theme=theme,order=order,TimeUnits=CodeToCodeUnitsTime,thermal=thermal,plot_limits=plot_limits_momentum,wide=false,legend=legend,framerate=framerate,filename=nothing,initial=initial,paraperp=paraperp,figure=(fig[2,1],time_idx))
        MomentumAndPolarAngleDistributionPlot(sol,species,PhaseSpace,Animated();order=order,theme=theme,framerate=framerate,filename=nothing,figure=(fig[1:3,2],time_idx))

        grid = fig[3,1] = GridLayout()

        t_unit_string = TimeUnits()
               
        Label(grid[1,1],@lift("t=$(round($(t), sigdigits = 3))"),fontsize=18pt)
        Label(grid[1,2],L"%$t_unit_string",fontsize=18pt)

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


############################ MISC Undocumented Plots for papers ###########################################
###########################################################################################################
###########################################################################################################

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
    xlab = L"$\log_{10}\left(p [m_ec]\right)$"
    ylab = L"$\log_{10}\left(p^2\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)]\right)$"
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot,:dense),(:dashdotdot,:dense)]
    legend_elements = []
    line_labels = []
    
    for i in 1:length(t_pho)

        t = t_pho[i]
        println("t=$t")
        if log10(t) % 1 == 0.0 && t <= t_max # 10^n timesteps
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
            if log10(t) % 1 == 0.0 && t <= t_max # 10^n timesteps 
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

function AM3_DIP_Combo_MomentumDistributionPlot(filePath_AM3,sol_DIP,PhaseSpace_DIP,t_max,t_min,t_grid;plot_limits=(nothing,nothing),theme=DiplodocusDark(),ele_err=true,pos_err=true)

    name_list = PhaseSpace_DIP.name_list
    Momentum = PhaseSpace_DIP.Momentum
    Grids = PhaseSpace_DIP.Grids
    Time = PhaseSpace_DIP.Time

    # load AM3 Data
    fileExist = isfile(filePath_AM3)
    if fileExist
        f = DC.jldopen(filePath_AM3,"r+");

        meanp_ele = f["meanp_ele"];
        f_ele = f["f_ele"];
        t_ele = f["t_ele"];

        meanp_pos = f["meanp_pos"];
        f_pos = f["f_pos"];
        t_pos = f["t_pos"];

        meanp_pho = f["meanp_pho"];
        f_pho = f["f_pho"];
        t_pho = f["t_pho"];

        DC.close(f)  
    else
        error("no file at $filePath_AM3 found")
    end

    # unit conversion 
    eV_to_mElec2 = 1.60217e-19 / 9.109e-31 / 2.9979e8^2
    cm3_to_m3 = 1e6 

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    fig = Figure(size=(576,460)) # double column 8:6 aspect ratio

    #g1 = fig[1,1] = GridLayout()
    g2 = fig[1,1]  = GridLayout()
    #g3 = fig[1,3] = GridLayout()

    #colsize!(g2.layout,0,Relative(0.05))
    #colsize!(g2.layout,1,Relative(0.9))
    #colsize!(g2.layout,2,Relative(0.05))

    xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
    ylab = L"$\log_{10}\left(p^2\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)]\right)$"
    ylab_err = L"$\text{Frac. Diff.}$"
    
    ax_DIP = Axis(g2[1,1],aspect=DataAspect(),ylabel=ylab)
    hidexdecorations!(ax_DIP,grid=false)
    ax_AM3 = Axis(g2[2,1],aspect=DataAspect())
    hidexdecorations!(ax_AM3,grid=false)
    ax_err = Axis(g2[3,1],xlabel=xlab,aspect=DataAspect(),ylabel=ylab_err)

    #pt = 4/3
    #Label(g2[1:2,0],ylab,rotation = pi/2,fontsize=9pt)
    #Label(g2[3,0],ylab_err,rotation = pi/2,fontsize=9pt)

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]
    legend_elements = []
    line_labels = []

    order = 2

    ######## Photon Lines ########
    ##############################

    species_name = "Pho"
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
    h_r = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]

    f3D = zeros(Float32,p_num,u_num,h_num)
    f1D = zeros(Float32,p_num*u_num*h_num)
    
    for i in 1:length(t_pho)

        # AM3
        t = t_pho[i]
        #println("t=$t")
        if log10(t) % 1 == 0.0 && t <= t_max # 10^n timesteps
            t_plot = t
            println("t=$t")

            if t_grid == "u"
                color = theme.colormap[][(t - t_min) / (t_max - t_min)]
            elseif t_grid == "l"
                color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
            end

            meanp_AM3 = meanp_pho .* eV_to_mElec2

            # AM3
            pdNdp_AM3 = meanp_AM3 .* f_pho[i,:][1] .* cm3_to_m3
            #println("$pdNdp")
            scatterlines!(ax_AM3,log10.(meanp_AM3),log10.(pdNdp_AM3),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])

            # DIP
            t_idx = findfirst(x->round(CodeToSIUnitsTime(x),sigdigits=4)==t_plot,sol_DIP.t)
            println(t_idx)
            f1D .= copy(Location_Species_To_StateVector(sol_DIP.f[t_idx],PhaseSpace_DIP,species_index=species_index))
            f3D .= reshape(f1D,(p_num,u_num,h_num))
            @. f3D = f3D*(f3D!=Inf)
            # scale by order
            # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
            for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
            end
            pdNdp_DIP = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
            if sum(@. !isnan(pdNdp_DIP) * !isinf(pdNdp_DIP) * !iszero(pdNdp_DIP)) == 1 # there is only one valid position so scatterlines doesn't work
                idx = findfirst(!iszero,pdNdp)
                lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp_DIP[idx])],linewidth=2.0,color = color,linestyle=linestyles[1])
            else
                scatterlines!(ax_DIP,log10.(meanp),log10.(pdNdp_DIP),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])
            end

            # error plot
            # regrid AM3 results (assuming it has finest grid) to grid size of DIP
            pdNdp_AM3_DIP_Grid = zeros(Float64,size(meanp))
            for p in eachindex(meanp)
                p_idx = find_closest(meanp_AM3,meanp[p])
                pdNdp_AM3_DIP_Grid[p] = pdNdp_AM3[p_idx]
            end
            err = (pdNdp_AM3_DIP_Grid.-pdNdp_DIP)./pdNdp_DIP
            replace!(err,-1.0=>NaN) # remove values for which AM3 does not have values
            replace!(err,Inf=>NaN) # remove values for which DIP does not have values

            scatterlines!(ax_err,log10.(meanp),err,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])


        end

    end

    push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
    push!(line_labels,"Pho")

    ######## Electron Lines ########
    ################################

    species_name = "Ele"
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
    h_r = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]

    f3D = zeros(Float32,p_num,u_num,h_num)
    f1D = zeros(Float32,p_num*u_num*h_num)

    for i in 1:length(t_ele)

        t = t_ele[i]
        if log10(t) % 1 == 0.0 && t <= t_max # 10^n timesteps 
            t_plot = t
            if t_grid == "u"
                color = theme.colormap[][(t - t_min) / (t_max - t_min)]
            elseif t_grid == "l"
                color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
            end
            meanp_AM3 = sqrt.((meanp_ele .* eV_to_mElec2).^2 .-1)

            # AM3
            pdNdp_AM3 = meanp_AM3 .* f_ele[i,:][1] .* cm3_to_m3
            scatterlines!(ax_AM3,log10.(meanp_AM3),log10.(pdNdp_AM3),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])

            # DIP
            t_idx = findfirst(x->round(CodeToSIUnitsTime(x),sigdigits=3)==t_plot,sol_DIP.t)
            f1D .= copy(Location_Species_To_StateVector(sol_DIP.f[t_idx],PhaseSpace_DIP,species_index=species_index))
            f3D .= reshape(f1D,(p_num,u_num,h_num))
            @. f3D = f3D*(f3D!=Inf)
            # scale by order
            # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
            for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
            end
            pdNdp_DIP = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
            if sum(@. !isnan(pdNdp_DIP) * !isinf(pdNdp_DIP) * !iszero(pdNdp_DIP)) == 1 # there is only one valid position so scatterlines doesn't work
                idx = findfirst(!iszero,pdNdp)
                lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp_DIP[idx])],linewidth=2.0,color = color,linestyle=linestyles[2])
            else
                scatterlines!(ax_DIP,log10.(meanp),log10.(pdNdp_DIP),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])
            end

            if ele_err
                # error plot
                # regrid AM3 results (assuming it has finest grid) to grid size of DIP
                pdNdp_AM3_DIP_Grid = zeros(Float64,size(meanp))
                for p in eachindex(meanp)
                    p_idx = find_closest(meanp_AM3,meanp[p])
                    pdNdp_AM3_DIP_Grid[p] = pdNdp_AM3[p_idx]
                end
                err = (pdNdp_AM3_DIP_Grid.-pdNdp_DIP)./pdNdp_DIP
                replace!(err,-1.0=>NaN) # remove values for which AM3 does not have values
                replace!(err,Inf=>NaN) # remove values for which DIP does not have values

                scatterlines!(ax_err,log10.(meanp),err,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])
            end

        end

    end

    push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[2],linewidth = 2.0))
    push!(line_labels,"Ele")

    ######## Positron Lines ########
    ################################

    species_name = "Pos"
    species_index = findfirst(x->x==species_name,name_list)
    if species_index != nothing
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
        h_r = Grids.pzr_list[species_index]
        mass = Grids.mass_list[species_index]

        f3D = zeros(Float32,p_num,u_num,h_num)
        f1D = zeros(Float32,p_num*u_num*h_num)

        for i in 1:length(t_pos)

            t = t_pos[i]
            if log10(t) % 1 == 0.0 && t <= t_max # 10^n timesteps 
                t_plot = t
                if t_grid == "u"
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                elseif t_grid == "l"
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                end
                meanp_AM3 = sqrt.((meanp_pos .* eV_to_mElec2).^2 .-1)

                # AM3
                pdNdp_AM3 = meanp_AM3 .* f_pos[i,:][1] .* cm3_to_m3
                scatterlines!(ax_AM3,log10.(meanp_AM3),log10.(pdNdp_AM3),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[3])

                # DIP
                t_idx = findfirst(x->round(CodeToSIUnitsTime(x),sigdigits=3)==t_plot,sol_DIP.t)
                f1D .= copy(Location_Species_To_StateVector(sol_DIP.f[t_idx],PhaseSpace_DIP,species_index=species_index))
                f3D .= reshape(f1D,(p_num,u_num,h_num))
                @. f3D = f3D*(f3D!=Inf)
                # scale by order
                # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
                for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                    f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
                end
                pdNdp_DIP = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
                if sum(@. !isnan(pdNdp_DIP) * !isinf(pdNdp_DIP) * !iszero(pdNdp_DIP)) == 1 # there is only one valid position so scatterlines doesn't work
                    idx = findfirst(!iszero,pdNdp)
                    lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp_DIP[idx])],linewidth=2.0,color = color,linestyle=linestyles[3])
                else
                    scatterlines!(ax_DIP,log10.(meanp),log10.(pdNdp_DIP),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[3])
                end

                if pos_err
                # error plot
                # regrid AM3 results (assuming it has finest grid) to grid size of DIP
                pdNdp_AM3_DIP_Grid = zeros(Float64,size(meanp))
                for p in eachindex(meanp)
                    p_idx = find_closest(meanp_AM3,meanp[p])
                    pdNdp_AM3_DIP_Grid[p] = pdNdp_AM3[p_idx]
                end
                err = (pdNdp_AM3_DIP_Grid.-pdNdp_DIP)./pdNdp_DIP
                replace!(err,-1.0=>NaN) # remove values for which AM3 does not have values
                replace!(err,Inf=>NaN) # remove values for which DIP does not have values

                scatterlines!(ax_err,log10.(meanp),err,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[3])
            end

            end
        end
        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[3],linewidth = 2.0))
        push!(line_labels,"Pos")
    end

    if t_grid == "u"
        Colorbar(g2[1:2,2],colormap = theme.colormap,limits=(TimeUnits(t_min),TimeUnits(t_max)),label=L"$t\,$ $[\text{s} * \sigma_{T}c]$",height=Relative(0.66))
    elseif t_grid == "l"
        Colorbar(g2[1:2,2],colormap = theme.colormap,limits=(log10(t_min),log10(t_max)),label=L"$\log_{10}\left(t \,[\text{s}]\right)$",height=Relative(0.66))
    end


    axislegend(ax_DIP,legend_elements,line_labels,position = :lt)

    linkxaxes!(ax_DIP,ax_err)
    linkxaxes!(ax_DIP,ax_AM3)
    #ax_DIP.limits = plot_limits
    xlims!(ax_DIP,plot_limits[1][1],plot_limits[1][2]) 
    ylims!(ax_DIP,plot_limits[2][1],plot_limits[2][2]) 
    #ax_AM3.limits = plot_limits
    xlims!(ax_AM3,plot_limits[1][1],plot_limits[1][2])
    ylims!(ax_AM3,plot_limits[2][1],plot_limits[2][2]) 

    xlims!(ax_err,plot_limits[1][1],plot_limits[1][2])
    ylims!(ax_err,-1.0,1.0)
    
    rowsize!(g2,1,Relative(0.4))
    rowsize!(g2,2,Relative(0.4))
    rowsize!(g2,3,Relative(0.2))

    rowgap!(g2,1,0.0)
    rowgap!(g2,2,0.0)
   
    return fig

    end # with_theme

end


function TwoSolAngleDistributionPlot(twosol::Tuple{OutputStruct,OutputStruct},species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,plot_limits=(nothing,nothing),wide=false,legend=true,framerate=12,filename="TwoSolAngleDistribution.mp4",figure=nothing)

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    xlab = L"$\log_{10}\left(p [m_ec]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}]\right)$"
    elseif order == 2
        ylab = L"$\log_{10}\left(p^{2}\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)]\right)$"
    else 
        ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
    end

    if isnothing(figure)
        time_idx = Observable(1) # index of the current time step
        t = @lift(TimeUnits(twosol[1].t[$time_idx]))
        if wide
            fig = Figure(size=(500,216)) # double column 8:3 aspect ratio
        else
            fig = Figure() # default single column 4:3 aspect ratio
        end
        #fig = Figure(size =(3.25inch,3.25inch)) # 1:1 aspect ratio
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    else
        fig, time_idx = figure # use the provided figure and time index
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    end

    ax.limits = plot_limits

    t_unit_string = TimeUnits()
               
    text!(ax,@lift("t=$(round($(t), sigdigits = 3))  "),fontsize=12pt,align=(:right,:center),space=:relative,offset=(375.0,135.0))
    text!(ax,L"%$t_unit_string",fontsize=12pt,align=(:left,:center),space=:relative,offset=(375.0,135.0))

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot,:dense),(:dashdotdot,:dense)]
    line_labels = []
    legend_elements = []  
    legend_elements_angle = []
    line_labels_angle = []  

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    counter = 1;

    for (species_idx, species_name) in enumerate(species) 

        species_index = findfirst(x->x==species_name,name_list)

        p_num = Momentum.px_num_list[species_index]  
        u_num = Momentum.py_num_list[species_index]
        h_num = Momentum.pz_num_list[species_index]
        dp = Grids.dpx_list[species_index]
        meanp = Grids.mpx_list[species_index]
        meanu = Grids.mpy_list[species_index]
        p_r = Grids.pxr_list[species_index]
        u_r = Grids.pyr_list[species_index]
        mass = Grids.mass_list[species_index]

        for i in 1:2

            color = theme.palette.color[][mod(2*i-1,7)+1]

            sol = twosol[i]

            pdNdp = @lift begin
                f1D = zeros(Float32,p_num*u_num*h_num)
                f1D .= copy(Location_Species_To_StateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index))
                f3D = zeros(Float32,p_num,u_num,h_num)
                f3D .= reshape(f1D,(p_num,u_num,h_num))
                @. f3D = f3D*(f3D!=Inf)
                # scale by order
                # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
                for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                    f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
                end
                # sum along u and h directions
                log10.(dropdims(sum(f3D, dims=(3)),dims=(3)))[:,2]
            end

            scatterlines!(ax,log10.(meanp),pdNdp,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])

            if legend==true && counter == 1 
                push!(legend_elements_angle,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
                if i == 1
                    push!(line_labels_angle,L"\text{Iso}")
                end
                if i == 2 
                    push!(line_labels_angle,L"\text{Ani}")
                end
            end

        end # sol loop 

        counter += 1

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
        push!(line_labels,species_name)

    end # species loop 

    if legend
        axislegend(ax,legend_elements,line_labels,position = :lt)
        axislegend(ax,legend_elements_angle,line_labels_angle,position = :lb)
    end
    
    if !isnothing(filename)
        # recording the animation
        time_idxs = 1:length(twosol[1].t)
        record(fig,filename,time_idxs,framerate=framerate,backend=CairoMakie) do frame
            println("$frame")
            time_idx[] = frame
        end
    end

    end # with_theme

end

function BFieldObserverPlot(sols::Vector{OutputStruct},PhaseSpaces::Vector{PhaseSpaceStruct},time_idx::Int64,ObserverAngle::Float64;ObserverDistance::Float64=1.0,theme=DiplodocusDark(),plot_limits=(nothing,nothing),TimeUnits::Function=CodeToCodeUnitsTime,title=nothing)

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

        fig = Figure(size=(432,216)) # 6:3 ratio
        t = round(TimeUnits(sols[1].t[time_idx]),sigdigits=3)
        t_unit_string = TimeUnits()
        xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
        ylab = L"$\log_{10}\left(pF_{p}\,[\text{m}^{-3}m_ec]\right)$"
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())

        if !isnothing(title)
            titlestr = L"t=%$(t)\, %$t_unit_string, \theta_\text{Obs}=%$ObserverAngle \pi"
            ax.title = titlestr
        end

        ax.limits = plot_limits

        max_total = -Inf32

        name_list = PhaseSpaces[1].name_list
        Space = PhaseSpaces[1].Space
        Momentum = PhaseSpaces[1].Momentum
        Grids = PhaseSpaces[1].Grids
        Time = PhaseSpaces[1].Time

        photon_index = findfirst(x->x=="Pho",name_list)

        ur = Grids.pyr_list[photon_index]
        pr = Grids.pxr_list[photon_index]
        mp = Grids.mpx_list[photon_index]
        dp = Grids.dpx_list[photon_index]

        legend_elements = []
        
        for sol in eachindex(sols)

            Fν = ObserverFlux(PhaseSpaces[sol],sols[sol],[ObserverAngle],ObserverDistance)
            νFν = log10.(mp .* Fν[time_idx,1,:])
            max_f = maximum(x for x in νFν if !isnan(x))
            max_total = max(max_total,max_f)

            if sol != length(sols)
                scatterlines!(ax,log10.(mp),νFν,color=theme.palette.color[][mod(sol,7)+1],markersize=0.0)
                push!(legend_elements,LineElement(color = theme.palette.color[][mod(sol,7)+1], linestyle = :solid,linewidth = 2.0))
            else
                scatterlines!(ax,log10.(mp),νFν,color=theme.textcolor[],markersize=0.0,linestyle=(:dash,:dense))
                push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = (:dash,:dense),linewidth = 2.0))
            end

            

        end

        if plot_limits == (nothing,nothing)
            xlims!(ax,(log10(pr[1]),log10(pr[end])))
            ylims!(ax,(max_total-9.0,max_total+1.0)) 
        end

        #line_labels=[L"\theta_\text{B}=0.0\pi",L"\theta_\text{B}=0.1\pi",L"\theta_\text{B}=0.2\pi",L"\theta_\text{B}=0.3\pi",L"\theta_\text{B}=0.4\pi",L"\theta_\text{B}=0.5\pi",L"\text{Iso}"]

        line_labels=[L"\theta_B=0",L"\theta_B=\pi/6",L"\theta_B=\pi/3",L"\theta_B=\pi/2",L"\text{Iso}"]

        axislegend(ax,legend_elements,line_labels,position = :lt)

        return fig

    end # with theme

end


function find_closest(A::AbstractArray{T}, b::T) where {T<:Real}
    if length(A) <= 1
        return firstindex(A)
    end

    i = searchsortedfirst(A, b)

    if i == firstindex(A)
        return i
    elseif i > lastindex(A)
        return lastindex(A)
    else
        prev_dist = b - A[i-1]
        next_dist = A[i] - b

        if prev_dist < next_dist
            return i - 1
        else
            return i
        end
    end
end