function ObserverFlux(PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,ObserverAngles::Vector{Float64},ObserverDistance::Float64)

    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    name_list = PhaseSpace.name_list

    if typeof(Space.space_coordinates) != Cylindrical
        error("ObserverFlux is only implemented for cylindrical coordinates.")
    end

    photon_index = findfirst(x->x=="Pho",name_list)

    β = Space.space_coordinates.β # local frame angle
    sβ, cβ = sincospi(β)

    d = ObserverDistance
    to_r = ObserverAngles

    dz = Grids.dz[1]
    R = Grids.xr[2]

    ur = Grids.pyr_list[photon_index]
    mp = Grids.mpx_list[photon_index]
    dp = Grids.dpx_list[photon_index]
    du = Grids.dpy_list[photon_index]

    Fν = zeros(length(sol.t),length(ObserverAngles),Momentum.px_num_list[photon_index])

    for (θ_idx, θ) in enumerate(ObserverAngles)

        sθ, cθ = sincospi(θ)

        a = sθ*sβ
        b = cθ*cβ

        u_low = cospi(θ+β)
        u_up = cospi(θ-β)

        #println("$a,$b")
        #println("$u_low,$u_up")
        #println("")

        for t in 1:length(sol.t)

            f1D = copy(Location_Species_To_StateVector(sol.f[t],PhaseSpace,species_index=photon_index))

            photon_f = reshape(f1D,(
                Momentum.px_num_list[photon_index],
                Momentum.py_num_list[photon_index],
                Momentum.pz_num_list[photon_index]
            ))

            for px in 1:(Momentum.px_num_list[photon_index]-1)

                val_total = 0.0

                for py in 1:(length(ur)-1)

                     if ur[py] <= u_low < ur[py+1] <= u_up  
                    u0 = u_low
                    u1 = ur[py+1]
                    elseif ur[py] <= u_low < u_up <= ur[py+1]
                        u0 = u_low
                        u1 = u_up
                    elseif u_low <= ur[py] < ur[py+1] <= u_up
                        u0 = ur[py]
                        u1 = ur[py+1]
                    elseif u_low <= ur[py] < u_up <= ur[py+1]
                        u0 = ur[py]
                        u1 = u_up
                    elseif u_low == u_up
                        u0 = u1 = 1.0 # i.e. no integration
                    elseif ur[py] <= ur[py+1] < u_low <= u_up
                        u0 = u1 = 1.0 # i.e. no integration
                    elseif  u_low < u_up <= ur[py] < ur[py+1] 
                        u0 = u1 = 1.0 # i.e. no integration
                    elseif  ur[py] < ur[py+1] <= u_low < u_up
                        u0 = u1 = 1.0 # i.e. no integration
                    elseif  u_low < u_up <= ur[py] < ur[py+1]
                        u0 = u1 = 1.0 # i.e. no integration
                    else
                        println("$u_low, $u_up, $(ur[py]), $(ur[py+1])")
                        error("Didn't account for this")
                    end

                    #println("$u_low, $u_up, $(ur[py]), $(ur[py+1]), $u0, $u1, $a, $b")
                    #println("$(a+b)")
                    #println("$β")
                    if u0 != u1 && u0 != u_low && u1 != u_up
                        val = mp[px] * photon_f[px,py,1] * (u1-u0) / dp[px] / du[py]
                    elseif u0 != u1 && u0 == u_low 
                        val = mp[px] * photon_f[px,py,1] * (u1-u0) / dp[px] / du[py]
                    elseif u0 != u1 && u1 == u_up
                        val = mp[px] * photon_f[px,py,1] * (u1-u0) / dp[px] / du[py]
                    else
                        val = 0.0
                    end

                    if β == 0.0
                        val = mp[px] * photon_f[px,py,1] * ((ur[py] <= cθ <= ur[py+1]) ? 1.0 : 0.0) * 2 * sθ / dp[px] / du[py]
                    elseif β == 1.0
                        val = mp[px] * photon_f[px,py,1] * ((ur[py] <= -cθ <= ur[py+1]) ? 1.0 : 0.0) * 2 / dp[px] / du[py] #* sθ
                    end

                    if val < 0
                        println("Negative value encountered: $val at px=$px, py=$py, θ=$θ")
                        println("$u_low, $u_up, $(ur[py]), $(ur[py+1]), $u0, $u1, $a, $b")
                        println("$(a+b)")
                        println("$β")
                        println("$(photon_f[px,py,1])")
                    end

                    val_total += val

                end  
                
                Fν[t,θ_idx,px] = val_total
                Fν[t,θ_idx,px] *= R*dz/(4*pi*d^2)

            end

        end

    end

    return Fν

end


"""
    ObserverFluxPlot(PhaseSpace,sol,time_idx,ObserverAngles,ObserverDistance;plot_limits=(nothing,nothing),theme=DiplodocusDark(),title=nothing)

"""
function ObserverFluxPlot(PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,time_idx::Int64,ObserverAngles::Vector{Float64},ObserverDistance::Float64;plot_limits=(nothing,nothing),theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

    CairoMakie.activate!(inline=true)

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    β = Space.space_coordinates.β # local frame angle

    photon_index = findfirst(x->x=="Pho",name_list)

    ur = Grids.pyr_list[photon_index]
    pr = Grids.pxr_list[photon_index]
    mp = Grids.mpx_list[photon_index]
    dp = Grids.dpx_list[photon_index]

    Fν = ObserverFlux(PhaseSpace,sol,ObserverAngles,ObserverDistance)

    fig = Figure()

    t = round(TimeUnits(sol.t[time_idx]),sigdigits=3)
    t_unit_string = TimeUnits()

    xlab = L"$\log_{10}\left(p [m_ec]\right)$"
    ylab = L"$\log_{10}\left(pF_{p}\right)$ $[\text{m}^{-3}m_ec]$"
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    
    ax.limits = plot_limits

    if !isnothing(title)
        titlestr = L" β =%$(β)\pi, t=%$(t) %$t_unit_string"
        ax.title = titlestr
    end

    max_total = -Inf32

    for θ in 1:length(ObserverAngles)

        flux_val = log10.(mp .* Fν[time_idx,θ,:])
        max_f = maximum(x for x in flux_val if !isnan(x))
        max_total = max(max_total,max_f)

        scatterlines!(ax,log10.(mp),flux_val,color=theme.palette.color[][mod(2*θ-1,7)+1],markersize=0.0,label= L"θ=%$(ObserverAngles[θ])\pi")

    end

    if plot_limits == (nothing,nothing)
        xlims!(ax,(log10(pr[1]),log10(pr[end])))
        ylims!(ax,(max_total-9.0,max_total+1.0)) 
    end

    axislegend(ax,location=:lt)

    return fig

    end # theme

end

