"""
    NumberDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark(),title=nothing)

Returns a plot of the number density of all species as a function of time.
"""
function NumberDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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

    num = zeros(Float64,length(sol.t))
    num_total = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits()

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=L"$n$ $[\mathrm{m}^{-3}]$",xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel=L"$n$ $[\mathrm{m}^{-3}]$",xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f1D = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))

            Nᵃ = DiplodocusTransport.FourFlow(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j])
            #Ua = HydroFourVelocity(Na)
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num[i] = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)

            if species== "All"
                num_total[i] += num[i]
            end
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    if species == "All"
        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),num_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),num_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end
    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    
    axislegend(ax,position = :rb)

    end # with_theme

    return fig

end



"""
    FracNumberDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark(),title=nothing)

Returns a plot of the fractional change in number density of all species between time setups as a function of time.
"""
function FracNumberDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,only_all::Bool=false)

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
        ax = Axis(fig[1,1],xlabel=xlab,ylabel="Num. Den. Frac. Change",xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel="Num. Den. Frac. Change",xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    frac_num = zeros(Float64,length(sol.t))
    num = 0.0
    num_total = zeros(Float64,length(sol.t))
    frac_num_total = zeros(Float64,length(sol.t))

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f1D = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))

            Nᵃ = DiplodocusTransport.FourFlow(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j])
            Uₐ = [-1.0,0.0,0.0,0.0]

            if i == 1
                frac_num[i] = 0.0 # initial value
            else
                frac_num[i] = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)/num - 1.0
            end

            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
        
            if species== "All"
                num_total[i] += num
            end

        end

        if !only_all
            if t_grid == "u"
                scatterlines!(ax,TimeUnits.(sol.t),frac_num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
                xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
            elseif t_grid == "l"
                scatterlines!(ax,log10.(TimeUnits.(sol.t)),frac_num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
                xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
            end
        end

    end

    if species == "All"

        for i in eachindex(sol.t)
            if i == 1
                frac_num_total[i] = 0.0
            else
                frac_num_total[i] = num_total[i]/num_total[i-1] - 1.0
            end
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),frac_num_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),frac_num_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end
    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,position = :rb)

    end # with_theme

    return fig
end

