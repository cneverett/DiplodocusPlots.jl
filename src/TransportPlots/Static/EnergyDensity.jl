"""
    EnergyDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the energy density of all species as a function of time.
"""
function EnergyDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark())

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

    eng = zeros(Float64,length(sol.t))

    if isnothing(fig)
        fig = Figure(size=(400,300))
        ax = Axis(fig[1,1],title="Energy Density",xlabel="Time",ylabel=L"$e/c$ $[\mathrm{J}\mathrm{m}^{-3}]$",xgridvisible=false,ygridvisible=false) # check units
    else
        ax = Axis(fig,title="Energy Density",xlabel="Time",ylabel=L"$e/c$ $[\mathrm{J}\mathrm{m}^{-3}]$",xgridvisible=false,ygridvisible=false) # check units
    end

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            Nᵃ = DiplodocusTransport.FourFlow(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
            Uₐ = [-1,0,0,0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)

            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j]) 

            eng[i] = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num)
        
        end

        if t_grid == "u"
            scatterlines!(ax,sol.t,eng,marker = :circle,markersize=0.0,label=name_list[j])
            xlims!(ax,sol.t[1],sol.t[end])
        elseif t_grid == "l"
            scatterlines!(ax,log10.(sol.t),eng,marker = :circle,markersize=0.0,label=name_list[j])
            xlims!(ax,log10(sol.t[1]),log10(sol.t[end]))
        end

    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,"Particles",position = :rb)

    end # with_theme

    return fig

end



"""
    FracEnergyDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the fractional change in energy density of all species between time setups as a function of time.
"""
function FracEnergyDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark())

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
        fig = Figure(size=(400,300))
        ax = Axis(fig[1,1],title="Frac. Change in Energy Density",xlabel="Time",ylabel="Frac. Change",xgridvisible=false,ygridvisible=false) # check units
    else
        ax = Axis(fig,title="Frac. Change in Energy Density",xlabel="Time",ylabel="Frac. Change",xgridvisible=false,ygridvisible=false) # check units
    end

    #j = findfirst(x->x==species,name_list)

    frac_eng = zeros(Float64,length(sol.t))
    eng = 0.0

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            Nᵃ = DiplodocusTransport.FourFlow(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
            Uₐ = [-1,0,0,0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j]) 

            if i == 1
                frac_eng[i] = 0.0 # initial value
            else
                frac_eng[i] = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num)/eng - 1.0
            end

            eng = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num)
        
        end

        if t_grid == "u"
            scatterlines!(ax,sol.t,frac_eng,marker = :circle,markersize=0.0,label=name_list[j])
            xlims!(ax,sol.t[1],sol.t[end])
        elseif t_grid == "l"
            scatterlines!(ax,log10.(sol.t),frac_eng,marker = :circle,markersize=0.0,label=name_list[j])
            xlims!(ax,log10(sol.t[1]),log10(sol.t[end]))
        end

    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,"Particles",position = :rb)

    end # with_theme

    return fig
end