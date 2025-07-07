"""
    MomentumDistributionPlot(sol,species,PhaseSpace;step=1,order=1,uDis=false,logt=false,plot_limits=(nothing,nothing),theme=DiplodocusDark())

Plots the angle averaged distribution function of of a given particle `species` as a function of time given by the `sol` based on the conditions held in `PhaseSpace`. 

Optional arguments:
- `step`: the step size in time to plot, default is 1.
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.
"""
function MomentumDistributionPlot(sol,species::String,PhaseSpace::PhaseSpaceStruct;step=1,order::Int64=1,thermal=false,uDis=false,logt=false,plot_limits=(nothing,nothing),theme=DiplodocusDark())

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    fig = Figure()
    xlab = L"$\log_{10}p$ $[m_\text{Ele}c]$"
    if order == 1
        ylab = L"$\log_{10}\left(p\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\right)$ $[\text{m}^{-3}]$"
    elseif order != 1
        ylab = L"$\log_{10}\left(p^{%$(order)}\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\right)$ $[\text{m}^{-3}]$"
    end
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
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
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    mass = Grids.mass_list[species_index]

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

    max_total = -Inf32

    for i in 1:t_save

        if i in values || i == 1

            t = sol.t[i]

            d .= reshape(sol.f[i].x[species_index],(p_num,u_num))

            max_d = maximum(x for x in d if !isnan(x))
            max_total = max(max_d,max_total)

            @. d = d*(d!=Inf)

            if uDis
                dj .= log10.(d .* dp)[:,4:3:8] # print just aligned and antialigned with field

                scatterlines!(ax,log10.(meanp),dj[:,1],linewidth=2.0,color = my_colors[color_counter],markersize=1.0)
                scatterlines!(ax,log10.(meanp),dj[:,2],linewidth=2.0,color = my_colors[color_counter],markersize=1.0,linestyle=:dash)
            else
                # sum along u direction
                pdNdp = dropdims(sum(d, dims=2),dims=2)
                # scale by order-1
                dj[:,1] .= log10.(pdNdp .* (meanp.^(order-1)))
                scatterlines!(ax,log10.(meanp),dj[:,1],linewidth=2.0,color = my_colors[color_counter],markersize=1.0)
            end

            color_counter += 1

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

        MJ = DiplodocusTransport.MaxwellJuttner_Distribution(PhaseSpace,"Sph",Temperature,mass;n=num)

        MJ = MJ .* (meanp.^(order-1))

        scatterlines!(ax,log10.(meanp),log10.(MJ),linewidth=2.0,color = :white,markersize=0.0,linestyle=:dash,label="Maxwell-Juttner")

    end

    if Time.t_grid == "u"
        Colorbar(fig[1,2],colormap = my_colors,limits=(sol.t[1],sol.t[end]),label=L"$t$ $[\text{s} * \sigma_{T}c]$")
    elseif Time.t_grid == "l"
        Colorbar(fig[1,2],colormap = my_colors,limits=(log10(sol.t[1]),log10(sol.t[end])),label=L"$\log_{10}(t)$ $[\text{s} * \sigma_{T}c]$")
    end

    if plot_limits == (nothing,nothing)
        xlims!(ax,(log10.(p_r[1]),log10.(p_r[end])))
        ylims!(ax,(max_total-11.0,max_total+1.0)) 
    end
    
    return fig

    end # with_theme

end


"""
    MomentumAndPolarAngleDistributionPlot(sol,species,PhaseSpace,time;order=1,theme=DiplodocusDark())

Plots the distribution function of of a given particle `species` as a function of momentum ``p`` and polar angle ``u`` at a given `time given by the `sol` based on the conditions held in `PhaseSpace`. 

Optional arguments:
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.
"""
function MomentumAndPolarAngleDistributionPlot(sol,species::String,PhaseSpace::PhaseSpaceStruct,time::T;order::Int64=1,theme=DiplodocusDark()) where T <: Union{Tuple{Float64,Float64,Float64},Tuple{Int64,Int64,Int64}}

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

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
    meanu = Grids.mpy_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    mass = Grids.mass_list[species_index]

    t_idx = zeros(Int64,3)
    t = zeros(Float64,3)

    if typeof(time) == Tuple{Int64,Int64,Int64}
        for i in 1:3
        t_idx[i] = time[i]
        t[i] = sol.t[t_idx[i]]
        end
    elseif typeof(time) == Tuple{Float64,Float64,Float64}
        for i in 1:3
        t[i] = time[i]
        t_idx[i] = findfirst(x->x==t[i],sol.t)
        end
    end

    fig = Figure(size=(800,300))

    dis1 = reshape(sol.f[t_idx[1]].x[species_index],(p_num,u_num))
    dis1 = dis1 .* (meanp.^(order-1))

    replace!(dis1,0.0 => NaN) # replace Inf with NaN for plotting
    dis2 = reshape(sol.f[t_idx[2]].x[species_index],(p_num,u_num))
    dis2 = dis2 .* (meanp.^(order-1))
    replace!(dis2,0.0 => NaN) # replace Inf with NaN for plotting
    dis3 = reshape(sol.f[t_idx[3]].x[species_index],(p_num,u_num))
    dis3 = dis3 .* (meanp.^(order-1))
    replace!(dis3,0.0 => NaN) # replace Inf with NaN for plotting

    max_dis = maximum(x for x in [dis1; dis2; dis3] if !isnan(x))
    min_dis = minimum(x for x in [dis1; dis2; dis3] if !isnan(x))
    col_range = (log10(max_dis)-20.0,log10(max_dis))

    ax1 = PolarAxis(fig[1,1+1],theta_0=-pi/2,direction=-1)
    ax1.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax1,0,pi)

    ax2 = PolarAxis(fig[1,2+1],theta_0=-pi/2,direction=-1)
    ax2.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax2,0,pi)

    ax3 = PolarAxis(fig[1,3+1],theta_0=-pi/2,direction=-1)
    ax3.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax3,0,pi)

    hm1 = heatmap!(ax1,acos.(u_r),log10.(p_r),log10.(dis1'),colorrange=col_range)
    hm2 = heatmap!(ax2,acos.(u_r),log10.(p_r),log10.(dis2'),colorrange=col_range)
    hm3 = heatmap!(ax3,acos.(u_r),log10.(p_r),log10.(dis3'),colorrange=col_range)

    rlims!(ax1,log10(p_r[1]),log10(p_r[end])+1.0)
    rlims!(ax2,log10(p_r[1]),log10(p_r[end])+1.0)
    rlims!(ax3,log10(p_r[1]),log10(p_r[end])+1.0)
    hidethetadecorations!(ax1, grid=false)
    hidethetadecorations!(ax2, grid=false)
    hidethetadecorations!(ax3, grid=false)

    if order == 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\right)$ $[\text{m}^{-3}]$",flipaxis=false,height=220,tellheight=false)
    elseif order != 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p^{%$order}\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\right)$ $[\text{m}^{-3}]$",flipaxis=false,height=220,tellheight=false)
    end

    text!(fig[1,2],L"$\log_{10}p$ $[m_\text{Ele}c]$",position=(-3.05,5.1),rotation=pi/2,fontsize=18.0f0)

    text!(fig[1,2],L"$t=%$(time[1])$",position=(2.6,7),fontsize=20.0f0)

    text!(fig[1,3],L"$t=%$(time[2])$",position=(2.6,7),fontsize=20.0f0)

    text!(fig[1,4],L"$t=%$(time[3])$",position=(2.6,7),fontsize=20.0f0)

    colgap!(fig.layout,1,Relative(-0.05))
    colgap!(fig.layout,2,Relative(-0.2))
    colgap!(fig.layout,3,Relative(-0.2))
    
    return fig

    end # with_theme

end