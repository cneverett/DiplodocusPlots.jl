"""
    BulkGamma2DPlot(type,sol,PhaseSpace,surface,species,time_idx;...)

Returns a density plot of the bulk gamma of a `species` at a given `time_idx` on a specific spatial coordinate `surface`.
"""
function BulkGamma2DPlot(type::Static,sol::OutputStruct,PhaseSpace::PhaseSpaceStruct,surface::String,species::String,time_idx::Int64;fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,wide=false,plot_limits=(nothing,nothing))

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
    xr = Grids.xr
    yr = Grids.yr
    zr = Grids.zr

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    if surface == "xy"
        xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        ylab = L"$\left(y\,[\text{Code Units}]\right)$"
        nx_dir = x_num
        ny_dir = y_num
        xgrid_type = Space.x_grid  
        ygrid_type = Space.y_grid
        if xgrid_type == "u"
            mxdir_plot = mx
            xr_dir = xr
            xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        elseif xgrid_type == "l"
            mxdir_plot = log10.(mx)
            xr_dir = log10.(xr)
            xlab = L"$\log_{10}\left(x\,[\text{Code Units}]\right)$"
        end
        if ygrid_type == "u"
            mydir_plot = my
            yr_dir = yr
            ylab = L"$\left(y\,[\text{Code Units}]\right)$"
        elseif ygrid_type == "l"
            mydir_plot = log10.(my)
            yr_dir = log10.(yr)
            ylab = L"$\log_{10}\left(y\,[\text{Code Units}]\right)$"
        end
    elseif surface == "yz"
        xlab = L"$\left(y\,[\text{Code Units}]\right)$"
        nx_dir = y_num
        ny_dir = z_num
        xgrid_type = Space.y_grid
        ygrid_type = Space.z_grid
        if xgrid_type == "u"
            mxdir_plot = my
            xr_dir = yr
            xlab = L"$\left(y\,[\text{Code Units}]\right)$"
        elseif xgrid_type == "l"
            mxdir_plot = log10.(my)
            xr_dir = log10.(yr)
            xlab = L"$\log_{10}\left(y\,[\text{Code Units}]\right)$"
        end
        if ygrid_type == "u"
            mydir_plot = mz
            yr_dir = zr
            ylab = L"$\left(z\,[\text{Code Units}]\right)$"
        elseif ygrid_type == "l"
            mydir_plot = log10.(mz)
            yr_dir = log10.(zr)
            ylab = L"$\log_{10}\left(z\,[\text{Code Units}]\right)$"
        end
    elseif surface == "xz"
        nx_dir = x_num
        ny_dir = z_num
        xgrid_type = Space.x_grid
        ygrid_type = Space.z_grid
        if xgrid_type == "u"
            mxdir_plot = mx
            xr_dir = xr
            xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        elseif xgrid_type == "l"
            mxdir_plot = log10.(mx)
            xr_dir = log10.(xr)
            xlab = L"$\log_{10}\left(x\,[\text{Code Units}]\right)$"
        end
        if ygrid_type == "u"
            mydir_plot = mz
            yr_dir = zr
            ylab = L"$\left(z\,[\text{Code Units}]\right)$"
        elseif ygrid_type == "l"
            mydir_plot = log10.(mz)
            yr_dir = log10.(zr)
            ylab = L"$\log_{10}\left(z\,[\text{Code Units}]\right)$"
        end
    else
        error("Invalid surface specified. Must be one of 'xy', 'yz', or 'xz'.")
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    ax.limits = plot_limits

    Gamma = zeros(Float64,nx_dir,ny_dir)

    species_index = findfirst(x->x==species,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]

    f1D = zeros(Float32,p_num*u_num*h_num)

    t = sol.t[time_idx]

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        Nᵃ = DiplodocusTransport.FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z)
        #Ua = HydroFourVelocity(Na)
        Uₐ = [-1.0,0.0,0.0,0.0] # static observer
        n = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)

        HydroFourVelocity = DiplodocusTransport.HydroFourVelocity(Nᵃ)

        println(HydroFourVelocity)

        Γ = -transpose(Uₐ) * HydroFourVelocity

        #println("n = $n at (x,y,z) = ($x,$y,$z)")

        if surface == "xy"
            Gamma[x,y] += Γ * dz[z] / (sum(dz))
        elseif surface == "yz"
            Gamma[y,z] += Γ * dx[x] / (sum(dx))
        elseif surface == "xz"
            Gamma[x,z] += Γ * dy[y] / (sum(dy))
        end

    end

    replace!(Gamma,0.0=>NaN)

    hm = heatmap!(ax,xr_dir,yr_dir,Gamma,colormap=theme.colormap_var)       

    t_unit_string = TimeUnits()

    Colorbar(fig[1,2],hm,label=L"$\Gamma\,$")

    end # with_theme

    return fig

end