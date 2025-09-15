"""
    TimeScalePlot(PhaseSpace, scheme, state)

Function plots the time scales for interaction loss rates of a simulation given the system `state` at a time index `t_idx`.

Optional arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `p_timescale`: whether to plot the timescale for momentum magnitude loss or state vector loss, default is `false` i.e. plot state vector losses to aid in assessing time step limits for stability.
"""
function TimeScalePlot(method::DiplodocusTransport.SteppingMethodType,state::Vector{F},t_idx::Int64;wide=false,paraperp::Bool=false,p_timescale=false,legend=true,horz_lines=nothing,plot_limits=(nothing,nothing),theme=DiplodocusDark(),TimeUnits::Function=CodeToCodeUnitsTime) where F<:AbstractFloat

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    PhaseSpace=method.PhaseSpace
    Momentum = PhaseSpace.Momentum
    Time = PhaseSpace.Time
    Grids = PhaseSpace.Grids
    tr = PhaseSpace.Grids.tr
    p_num_list = Momentum.px_num_list 
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    dp_list = Grids.dpx_list
    mp_list = Grids.mpx_list
    name_list = PhaseSpace.name_list

    # only works for a single particle
    #state = reshape(state,(p_num_list[1],u_num_list[1],h_num_list[1]))
    #state = mp_list[1] .* state
    #state = reshape(state,p_num_list[1]*u_num_list[1]*h_num_list[1])

    dstate = zeros(eltype(state),size(state))
    timescale = zeros(eltype(state),size(state))

    dt0 = tr[2] - tr[1]
    dt = tr[t_idx+1] - tr[t_idx]
    t = tr[t_idx]

    method(dstate,state,dt0,dt,t)

    dstate = DiplodocusTransport.diag(method.temp) .* state

    #dstate = method.FluxM.Ap_Flux \ (DiplodocusTransport.diag(method.FluxM.I_Flux .+ method.FluxM.J_Flux) .* (dt / dt0)) .* state

    @. timescale =  -dt * state / dstate / 0.5868763104768393 # REMOVE THIS LATER

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    t_unit_string = TimeUnits()

    xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
    if p_timescale
        ylab = L"$\log_{10}\left(\frac{p}{\mathrm{d}p/\mathrm{d}t}\,%$t_unit_string \right)$"
    else
        ylab = L"$\log_{10}\left(Timescale\,%$t_unit_string\right)$"
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]
    legend_elements_species = []
    line_labels_species = []
    legend_elements_angle = []
    line_labels_angle = []

    if !isnothing(horz_lines)

        for i in eachindex(horz_lines)

            if Time.t_grid == "u"
                t_horz = horz_lines[i]
                color_val = (t_horz - TimeUnits(tr[1]))/(TimeUnits(tr[end])-TimeUnits(tr[1]))
                println("colval=$color_val")
                color = theme.colormap[][color_val]
            elseif Time.t_grid == "l"
                t_horz = horz_lines[i]
                color_val = (t_horz - log10(TimeUnits(tr[1])))/(log10(TimeUnits(tr[end]))-log10(TimeUnits(tr[1])))
                color = theme.colormap[][color_val]
            end

            hlines!(ax,log10.(horz_lines[i]),color=color,linewidth=1.5)

        end

    end

    for species in eachindex(name_list)

        name = name_list[species]
        p_num = p_num_list[species]
        u_num = u_num_list[species]
        h_num = h_num_list[species]
        dp = dp_list[species]
        mp = mp_list[species]
        mu = Grids.mpy_list[species]

        timescale1D = copy(Location_Species_To_StateVector(timescale,PhaseSpace,species_index=species))
        timescale3D = reshape(timescale1D,(p_num,u_num,h_num))
        @. timescale3D = timescale3D*(timescale3D!=Inf)
        @. timescale3D = timescale3D*(timescale3D!=-Inf)
        #@. timescale3D = timescale3D*(timescale3D<=0.0)

        timescale2D = dropdims(sum(timescale3D, dims=(3)),dims=(3))
        # timescale is timescale for particle losses not particle energy losses, to get that we need to scale by dp/mp (i.e. dp/p)
        if p_timescale
            timescale2D = mp ./ dp .* timescale2D
        end

        if paraperp==true

            println("tscale=$(TimeUnits.(Float64.(abs.(timescale2D[10,ceil(Int64,u_num/2)]))))")

            scatterlines!(ax,log10.(mp),log10.(TimeUnits.(Float64.(abs.(timescale2D[:,ceil(Int64,u_num/2)])))),linewidth=2.0,color = color=theme.textcolor[],markersize=0.0,linestyle=linestyles[1])
            scatterlines!(ax,log10.(mp),log10.(TimeUnits.(Float64.(abs.(timescale2D[:,end])))),linewidth=2.0,color = color=theme.textcolor[],markersize=0.0,linestyle=linestyles[2])

            push!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
            push!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = linestyles[2],linewidth = 2.0))
            push!(line_labels_angle,L"\parallel")
            push!(line_labels_angle,L"\perp")

        else
            for u in 1:u_num

                if u == 1 || meanu[u]==0.0
                    scatterlines!(ax,log10.(mp),log10.(TimeUnits.(Float64.(timescale2D[:,u]))),linewidth=2.0,color = theme.palette.color[][mod(2*u-1,7)+1],markersize=0.0,linestyle=linestyles[species])
                end

                if species == 1
                    push!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = :solid,linewidth = 2.0))
                    push!(line_labels_angle,L"%$(mu[u])")
                end

            end 
        end

        #push!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = :solid,linewidth = 2.0))
        #push!(line_labels_angle,L"%$(mu[u])")

        push!(legend_elements_species,LineElement(color = theme.textcolor[], linestyle = linestyles[species],linewidth = 2.0))
        push!(line_labels_species,name)

    end # species loop

        if paraperp
            axislegend(ax,legend_elements_angle,line_labels_angle,position = :lb)
        end
        if legend
            axislegend(ax,legend_elements_species,line_labels_species,position = :lb)
        end

    return fig

    end # with theme

end