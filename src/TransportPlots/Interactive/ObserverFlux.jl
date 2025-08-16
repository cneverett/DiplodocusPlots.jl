function ObserverFluxPlot(PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,ObserverAngles::Vector{Float64},ObserverDistance::Float64;plot_limits=(nothing,nothing),theme=DiplodocusDark(),title=nothing)

    with_theme(theme) do

    GLMakie.activate!(inline=false)

    name_list = PhaseSpace.name_list
    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    β = Space.space_coordinates.β # local frame angle

    photon_index = findfirst(x->x=="Pho",name_list)

    ur = Grids.pyr_list[photon_index]
    mp = Grids.mpx_list[photon_index]

    Fν = ObserverFlux(PhaseSpace,sol,ObserverAngles,ObserverDistance)

    fig = Figure()
    
    sg = SliderGrid(fig[2,1],
    (label = "t_idx", range = 1:length(sol.t), startvalue = 1, update_while_dragging = false),
    )

    t_idx = sg.sliders[1].value

    t = @lift(sol.t[$t_idx])

    t_v = t[]

    ax = Axis(fig[1,1],ylabel=L"$\log_{10}\left(pF_{p}\right)$ $[\text{m}^{-3}]$",aspect=DataAspect())
    ax.limits = plot_limits

    if !isnothing(title)
        titlestr = @lift("Observer Flux at distance $ObserverDistance, with B-Field at an angle of β =$(β)π, at t=$(sol.t[$t_idx])")
        ax.title = titlestr
    end

    for θ in 1:length(ObserverAngles)

        flux_val = @lift(log10.(mp .* Fν[$t_idx,θ,:]))

        scatterlines!(ax,log10.(mp),flux_val,color=theme.palette.color[][mod(2*θ-1,7)+1],markersize=0.0,label= "θ=$(ObserverAngles[θ])π")

    end

    axislegend(ax)

    end # theme

    return fig

end