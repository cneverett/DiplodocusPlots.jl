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
    elseif order == 2
        ylab = L"$\log_{10}\left(p^2\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\right)$ $[\text{m}^{-3}]$"
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
                if order == 1
                    dj[:,1] .= log10.(pdNdp) 
                elseif order == 2
                    dj[:,1] .= log10.(pdNdp .* meanp)
                end
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

        if order == 2
            MJ = MJ .* meanp
        end

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