"""
    TimeScalePlot(PhaseSpace, scheme, state)

Function plots the time scales of a simulation given the system `state` at a time index `t_idx`.
"""
function TimeScalePlot(method::DiplodocusTransport.SteppingMethodType,state::Vector{F},t_idx::Int64;wide=false,plot_limits=(nothing,nothing),theme=DiplodocusDark()) where F<:AbstractFloat

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    PhaseSpace=method.PhaseSpace
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    tr = PhaseSpace.Grids.tr
    mp_list = Grids.mpx_list
    name_list = PhaseSpace.name_list

    dstate = zeros(eltype(state),size(state))
    timescale = zeros(eltype(state),size(state))

    dt0 = tr[2] - tr[1]
    dt = tr[t_idx+1] - tr[t_idx]
    t = tr[t_idx]

    method(dstate,state,dt0,dt,t)

    @. timescale = dstate / state / dt

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    xlab = L"$\log_{10}\left(p [m_ec]\right)$"
    ylab = L"$\tau$"

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot,:dense),(:dashdot,:dense),(:dashdotdot,:dense)]
    legend_elements_species = []
    line_labels_species = []
    legend_elements_angle = []
    line_labels_angle = []

    for species in eachindex(name_list)

        name = name_list[species]
        p_num = Momentum.px_num_list[species]  
        u_num = Momentum.py_num_list[species]
        h_num = Momentum.pz_num_list[species]
        meanp = Grids.mpx_list[species]
        meanu = Grids.mpy_list[species]

        timescale1D = copy(Location_Species_To_StateVector(timescale,PhaseSpace,species_index=species))
        timescale3D = reshape(timescale1D,(p_num,u_num,h_num))
        @. timescale3D = timescale3D*(timescale3D!=Inf)

        timescale2D = dropdims(sum(timescale3D, dims=(3)),dims=(3))

        for u in 1:u_num

            println(timescale2D[:,u])
        
            scatterlines!(ax,log10.(meanp),log10.(abs.(timescale2D[:,u])),linewidth=2.0,color = color=theme.palette.color[][mod(u,7)],markersize=0.0,linestyle=linestyles[species])

            if species == 1
                push!(legend_elements_angle,LineElement(color = theme.palette.color[][mod(u,7)], linestyle = :solid,linewidth = 2.0))
                push!(line_labels_angle,L"%$(meanu[u])")
            end

        end

        push!(legend_elements_species,LineElement(color = theme.textcolor[], linestyle = linestyles[species],linewidth = 2.0))
        push!(line_labels_species,name)

    end # species loop

    end # with theme

end