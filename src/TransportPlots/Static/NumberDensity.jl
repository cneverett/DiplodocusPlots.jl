"""
    NumberDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the number density of all species as a function of time.
"""
function NumberDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark())

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    t_grid = Time.t_grid

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list

    mass_list = Grids.mass_list

    num = zeros(Float64,length(sol.t))

    if isnothing(fig)
        fig = Figure(size=(600,300))
        ax = Axis(fig[1,1],title="Number Density",xlabel="Time",ylabel=L"$n$ $[\mathrm{m}^{-3}]$")
    else
        ax = Axis(fig,title="Number Density",xlabel="Time",ylabel=L"$n$ $[\mathrm{m}^{-3}]$")
    end

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            Nᵃ = DiplodocusTransport.FourFlow(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
            #Ua = HydroFourVelocity(Na)
            Uₐ = [-1,0,0,0] # static observer
            num[i] = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
        
        end

        if t_grid == "u"
            scatterlines!(ax,sol.t,num,marker = :circle,markersize=0.0,label=name_list[j])
            xlims!(ax,sol.t[1],sol.t[end])
        elseif t_grid == "l"
            scatterlines!(ax,log10.(sol.t),num,marker = :circle,markersize=0.0,label=name_list[j])
            xlims!(ax,log10(sol.t[1]),log10(sol.t[end]))
        end

    end

    fig[1,2] = Legend(fig,ax,"Particles") 

    end # with_theme

    return fig

end



"""
    FracNumberDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the fractional change in number density of all species between time setups as a function of time.
"""
function FracNumberDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark())

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    t_grid = Time.t_grid

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list

    mass_list = Grids.mass_list

    if isnothing(fig)
        fig = Figure(size=(600,300))
        ax = Axis(fig[1,1],title="Frac. Change in Number Density",xlabel="Time",ylabel="Frac. Change")
    else
        ax = Axis(fig,title="Frac. Change in Number Density",xlabel="Time",ylabel="Frac. Change")
    end

    frac_num = zeros(Float64,length(sol.t))
    num = 0.0

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            Nᵃ = DiplodocusTransport.FourFlow(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
            Uₐ = [-1,0,0,0]

            if i == 1
                frac_num[i] = 0.0 # initial value
            else
                frac_num[i] = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)/num - 1.0
            end

            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
        
        end

        if t_grid == "u"
            scatterlines!(ax,sol.t,frac_num,marker = :circle,markersize=0.0,label=name_list[j])
            xlims!(ax,sol.t[1],sol.t[end])
        elseif t_grid == "l"
            scatterlines!(ax,log10.(sol.t),frac_num,marker = :circle,markersize=0.0,label=name_list[j])
            xlims!(ax,log10(sol.t[1]),log10(sol.t[end]))
        end

    end

    fig[1,2] = Legend(fig,ax,"Particles")

    end # with_theme

    return fig
end

