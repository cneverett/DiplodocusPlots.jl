"""
    NumberDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark(),title=nothing)

Returns a plot of the number density of all species as a function of time.
"""
function NumberDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,x_idx=1,y_idx=1,z_idx=1)

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

            f1D = copy(Location_Species_To_StateVector(sol.f[i],PhaseSpace,species_index=j,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))

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
function FracNumberDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),frac_num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),frac_num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,position = :rb)

    end # with_theme

    return fig
end

"""
    NumberDensity1DPlot(sol,PhaseSpace,direction;species="All",fig=nothing,theme=DiplodocusDark(),title=nothing)

Returns a plot of the number density of all species as a function of time along a specific spatial coordinate direction*.
"""
function NumberDensity1DPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct,direction,type::Static,species::Vector{String};fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,step::Int=1,wide=false,logt::Bool=false,plot_limits=(nothing,nothing))

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Grids = PhaseSpace.Grids
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    dx = Grids.dx
    dy = Grids.dy
    dz = Grids.dz
    mx = Grids.mx
    my = Grids.my
    mz = Grids.mz

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    ylab = L"$n\, [\mathrm{m}^{-3}]$"

    if direction == "x"
        xlab = L"$\left(x\,[\mathrm{m}]\right)$"
        n_dir = x_num
        mdir_plot = mx
    elseif direction == "y"
        xlab = L"$\left(y\,[\mathrm{m}]\right)$"
        n_dir = y_num
        mdir_plot = my
    elseif direction == "z"
        xlab = L"$\left(z\,[\mathrm{m}]\right)$"
        n_dir = z_num
        mdir_plot = mz
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]

    legend_elements = []
    line_labels = []

    t_save = length(sol.t)
    t_plot = ceil(Int64,t_save/step)

    values = (1:t_save)*step .+ 2 # add 2 to skip initial and kernel steps

    num = zeros(Float64,n_dir)
    num_total = zeros(Float64,n_dir)

    for (species_idx, species_name) in enumerate(species) 

        species_index = findfirst(x->x==species_name,name_list)

        p_num = Momentum.px_num_list[species_index]  
        u_num = Momentum.py_num_list[species_index]
        h_num = Momentum.pz_num_list[species_index]
        pr = Grids.pxr_list[species_index]
        ur = Grids.pyr_list[species_index]
        hr = Grids.pzr_list[species_index]
        mass = Grids.mass_list[species_index]

        f1D = zeros(Float32,p_num*u_num*h_num)

        t_min = logt ? sol.t[2]/10 : sol.t[1]
        t_max = sol.t[end]

        for i in 1:t_save

            if (i in values || i == 1 || i == 2) # plot first step for initial conds, second for kernel 

                t = sol.t[i]
                #println("t=$(CodeToSIUnitsTime(t))")
                if Time.t_grid == "l" || logt
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                elseif Time.t_grid == "u"
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                end

                # sum over other 2 spatial directions to get average number density along desired direction
                fill!(num,0.0)

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    Nᵃ = DiplodocusTransport.FourFlow(sol.f[i],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z)
                    #Ua = HydroFourVelocity(Na)
                    Uₐ = [-1.0,0.0,0.0,0.0] # static observer
                    n = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)

                    if direction == "x"
                        num[x] += n * dy[y] * dz[z] / (sum(dy)*sum(dz))
                        if species== "All"
                            num_total[x] += num[x]
                        end
                    elseif direction == "y"
                        num[y] += n * dx[x] * dz[z] / (sum(dx)*sum(dz))
                        if species== "All"
                            num_total[y] += num[y]
                        end
                    elseif direction == "z"
                        num[z] += n * dx[x] * dy[y] / (sum(dx)*sum(dy))
                        if species== "All"
                            num_total[z] += num[z]
                        end
                    end

                end

                scatter!(ax,mdir_plot,num,color = color,markersize=5.0)

            end # if times should be plotted
        
        end # time loop
       
    end # species loop 

    t_unit_string = TimeUnits()

    if Time.t_grid == "u"
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(TimeUnits(sol.t[1]),TimeUnits(sol.t[end])),label=L"$t\,$ $%$t_unit_string$")
    elseif Time.t_grid == "l"
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(round(TimeUnits(sol.t[1]),sigdigits=5)),log10(round(TimeUnits(sol.t[end]),sigdigits=5))),label=L"$\log_{10}\left(t\,%$t_unit_string \right)$")
    end
    
    #axislegend(ax,position = :rb)

    end # with_theme

    return fig

end