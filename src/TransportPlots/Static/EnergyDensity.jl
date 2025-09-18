"""
    EnergyDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark(),title=nothing)

Returns a plot of the energy density of all species as a function of time.
"""
function EnergyDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,perparticle=false)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    t_grid = Time.t_grid

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list
    hr_list = Grids.pzr_list

    mass_list = Grids.mass_list

    eng = zeros(Float64,length(sol.t))
    eng_total = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits()

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
    end

    if perparticle
        ylab = L"$\bar{p^0}$ $[\mathrm{J}/c]$"
    else
        ylab = L"$e/c$ $[\mathrm{J}\mathrm{m}^{-3}]$"
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false) # check units
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false) # check units
    end

    if !isnothing(title)
        ax.title = title
    end

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f1D = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))

            Nᵃ = DiplodocusTransport.FourFlow(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j])
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)

            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j]) 

            eng[i] = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num,perparticle=perparticle)

            if species== "All"
                eng_total[i] += eng[i]
            end
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),eng,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),eng,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    if species == "All"
        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end
    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,position = :rc)

    end # with_theme

    return fig

end



"""
    FracEnergyDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark(),title=nothing,only_all=false)

Returns a plot of the fractional change in energy density of all species between time setups as a function of time.
"""
function FracEnergyDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,only_all=false,TimeUnits::Function=CodeToCodeUnitsTime)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    t_grid = Time.t_grid

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list
    hr_list = Grids.pzr_list

    mass_list = Grids.mass_list

    t_unit_string = TimeUnits()

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel="Eng. Den. Frac. Change",xgridvisible=false,ygridvisible=false) # check units
    else
        ax = Axis(fig,xlabel=xlab,ylabel="Eng. Den. Frac. Change",xgridvisible=false,ygridvisible=false) # check units
    end

    if !isnothing(title)
        ax.title = title
    end

    #j = findfirst(x->x==species,name_list)

    frac_eng = zeros(Float64,length(sol.t))
    frac_eng_total = zeros(Float64,length(sol.t))
    eng = 0.0
    eng_total = zeros(Float64,length(sol.t))

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f1D = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))

            Nᵃ = DiplodocusTransport.FourFlow(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j])
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j]) 

            if i == 1
                frac_eng[i] = 0.0 # initial value
            else
                frac_eng[i] = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num)/eng - 1.0
            end

            eng = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num)

            if species== "All"
                eng_total[i] += eng
            end
        
        end

        if !only_all
            if t_grid == "u"
                scatterlines!(ax,TimeUnits.(sol.t),frac_eng,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
                xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
            elseif t_grid == "l"
                scatterlines!(ax,log10.(TimeUnits.(sol.t)),frac_eng,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
                xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
            end
        end

    end

    if species == "All"
        for i in eachindex(sol.t)
            if i == 1
                frac_eng_total[i] = 0.0 # initial value
            else
                frac_eng_total[i] = eng_total[i]/eng_total[i-1] - 1.0
            end
        end
        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),frac_eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),frac_eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end
    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,position = :rb)

    end # with_theme

    return fig
end