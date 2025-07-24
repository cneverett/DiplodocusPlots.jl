"""
    ObserverFluxPlot(PhaseSpace,sol,time_idx,ObserverAngles,ObserverDistance;plot_limits=(nothing,nothing),theme=DiplodocusDark(),title=nothing)

"""
function ObserverFluxPlot(PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,time_idx::Int64,ObserverAngles::Vector{Float64},ObserverDistance::Float64;plot_limits=(nothing,nothing),theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

    CairoMakie.activate!(inline=false)

    with_theme(theme) do

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

    t = TimeUnits(sol.t[time_idx])
    t_unit_string = TimeUnits()

    ax = Axis(fig[1,1],ylabel=L"$\log_{10}\left(pF_{p}\right)$ $[\text{m}^{-3}]$",aspect=DataAspect())
    ax.limits = plot_limits

    if !isnothing(title)
        titlestr = L"Observer Flux at distance %$ObserverDistance [\text{m}], with B-Field at an angle of β =%$(β)\pi, at t=%$(t) %$t_unit_string"
        ax.title = titlestr
    end

    for θ in 1:length(ObserverAngles)

        flux_val = log10.(mp .* Fν[time_idx,θ,:])

        scatterlines!(ax,log10.(mp),flux_val,color=theme.palette.color[][mod(2*θ-1,7)+1],markersize=0.0,label= L"θ=%$(ObserverAngles[θ])\pi")

    end

    axislegend(ax)

    return fig

    end # theme

end