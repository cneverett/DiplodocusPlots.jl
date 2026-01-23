"""
    function IsThermalPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the sum of squared residuals between the distribution function for each species and an expected thermal distribution based on the current distribution of that species as a function of time.  
"""
function IsThermalPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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

    SumSquaredResiduals = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits()

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
        ylab = L"SSR"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
        ylab = L"\log_{10}(SSR)"
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            name = name_list[j]
            mass = mass_list[j]

            f = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))
            Nᵃ = DiplodocusTransport.FourFlow(f,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass)
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
            Δab = DiplodocusTransport.ProjectionTensor(Uₐ)
            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass)
            Pressure = DiplodocusTransport.ScalarPressure(Tᵃᵇ,Δab)
            Temperature = DiplodocusTransport.ScalarTemperature(Pressure,num)

            if mass != 0.0
                Ther = DiplodocusTransport.MaxwellJuttner_Distribution(PhaseSpace,name,Temperature;n=num)
            else
                Ther = DiplodocusTransport.BlackBody_Distribution(PhaseSpace,name,Temperature;n=num)
            end

            f1D = dropdims(sum(reshape(f,(p_num_list[j],u_num_list[j],h_num_list[j])),dims=(2,3)),dims=(2,3))

            residuals = (f1D .- Ther).^2
            residuals = residuals[isfinite.(residuals)]

            SumSquaredResiduals[i] = sum(residuals)
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),SumSquaredResiduals,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    ylims!(ax,0.0,nothing)

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,"Particles")

    end # with_theme

    return fig

end

"""
    function IsIsotopicPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the sum of squared residuals between the distribution function for each species and an expected isotropic distribution based on the current distribution of that species as a function of time.  
"""
function IsIsotropicPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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

    SumSquaredResiduals = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits()
    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
        ylab = L"SSR"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
        ylab = L"\log_{10}(SSR)"
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))
            f3D = reshape(f,(p_num_list[j],u_num_list[j],h_num_list[j]))
            # sum over angles and divide by number of bins
            f_avg = dropdims(sum(f3D,dims=(2,3)),dims=(2,3)) / (u_num_list[j]*h_num_list[j])
 
            residuals = (f3D .- f_avg).^2
            residuals = residuals[isfinite.(residuals)]

            SumSquaredResiduals[i] = sum(residuals)
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),SumSquaredResiduals,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    ylims!(ax,0.0,nothing)

    axislegend(ax)

    end # with_theme

    return fig

end

"""
    function IsThermalAndIsotropicPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the sum of squared residuals between the distribution function for each species and an expected Maxwell-Juttner distribution and an isotropic distribution based on the current distribution of that species as a function of time.  
"""
function IsThermalAndIsotropicPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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

    SumSquaredResiduals = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits()

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
        ylab = L"SSR"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
        ylab = L"\log_{10}(SSR)"
    end

    if isnothing(fig)
        fig = Figure(size=(288,216))
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    # Is Thermal?
    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))
            Nᵃ = DiplodocusTransport.FourFlow(f,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j])
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
            Δab = DiplodocusTransport.ProjectionTensor(Uₐ)
            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j])
            Pressure = DiplodocusTransport.ScalarPressure(Tᵃᵇ,Δab)
            Temperature = DiplodocusTransport.ScalarTemperature(Pressure,num)

            if mass_list[j] != 0.0
                Ther = DiplodocusTransport.MaxwellJuttner_Distribution(PhaseSpace,name_list[j],Temperature;n=num)
            else
                Ther = DiplodocusTransport.BlackBody_Distribution(PhaseSpace,name_list[j],Temperature;n=num)
            end

            f1D = dropdims(sum(reshape(f,(p_num_list[j],u_num_list[j],h_num_list[j])),dims=(2,3)),dims=(2,3))

            residuals = (f1D .- Ther).^2
            residuals = residuals[isfinite.(residuals)]

            SumSquaredResiduals[i] = sum(residuals)
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),SumSquaredResiduals,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j]*" Thermal?")
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j]*" Thermal?")
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    # Is Isotropic?
    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))
            f3D = reshape(f,(p_num_list[j],u_num_list[j],h_num_list[j]))
            # sum over angles and divide by number of bins
            f_avg = dropdims(sum(f3D,dims=(2,3)),dims=(2,3)) / (u_num_list[j]*h_num_list[j])
 
            residuals = (f3D .- f_avg).^2
            residuals = residuals[isfinite.(residuals)]

            SumSquaredResiduals[i] = sum(residuals)
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),SumSquaredResiduals,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,linestyle=:dash,label=name_list[j]*" Isotropic?")
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,linestyle=:dash,label=name_list[j]*" Isotropic?")
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    #ylims!(ax,0.0,nothing)
    
    axislegend(ax,position=:lb)

    end # with_theme

    return fig

end

"""
    function TemperaturePlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the temperature associated with each particles distribution function as a function of time.  
"""
function ThermalPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    name_list_plot = similar(name_list)
    for (idx, name) in enumerate(name_list)
        if name == "Ele"
            name = "Electron"
        elseif name == "Pho"
            name = "Photon"
        elseif name == "Pos"
            name = "Positron"
        else
            name = name
        end
        name_list_plot[idx] = name
    end
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

    Temperature = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits()

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
        ylab = L"$T$ $[\mathrm{K}]$"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
        ylab = L"$T$ $[\mathrm{K}]$"
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            name = name_list[j]
            mass = mass_list[j]

            f = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j))
            Nᵃ = DiplodocusTransport.FourFlow(f,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass)
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
            Δab = DiplodocusTransport.ProjectionTensor(Uₐ)
            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass)
            Pressure = DiplodocusTransport.ScalarPressure(Tᵃᵇ,Δab)
            Temperature[i] = DiplodocusTransport.ScalarTemperature(Pressure,num)
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),Temperature,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),Temperature,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    ylims!(ax,0.0,nothing)

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,position=:lb)

    end # with_theme

    return fig

end