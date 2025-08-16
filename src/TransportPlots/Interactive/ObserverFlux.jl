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
                        val = mp[px] * photon_f[px,py,1] * (u1-u0)#* (-atan(sqrt(a^2-(b-u1)^2),(u1-b))+atan(sqrt(a^2-(b-u0)^2),(u0-b)))
                    elseif u0 != u1 && u0 == u_low 
                        val = mp[px] * photon_f[px,py,1] * (u1-u0) #* (-atan(sqrt(a^2-(b-u1)^2),(u1-b)) + (sign(u0-b)==1 ? 0 : pi))
                    elseif u0 != u1 && u1 == u_up
                        val = mp[px] * photon_f[px,py,1] * (u1-u0) #* (-(sign(u1-b)==1 ? 0 : pi)+atan(sqrt(a^2-(b-u0)^2),(u0-b)))
                    else
                        val = 0.0
                    end

                    if β == 0.0
                        val = mp[px] * photon_f[px,py,1] * ((ur[py] <= cθ <= ur[py+1]) ? 1.0 : 0.0) * 2 * sθ 
                    elseif β == 1.0
                        val = mp[px] * photon_f[px,py,1] * ((ur[py] <= -cθ <= ur[py+1]) ? 1.0 : 0.0) * 2 #* sθ
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


function ObserverFluxPlot(PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,ObserverAngles::Vector{Float64},ObserverDistance::Float64;plot_limits=(nothing,nothing),theme=DiplodocusDark(),title=nothing)

    with_theme(theme) do

    GLMakie.activate!(inline=false)

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
    
    sg = SliderGrid(fig[2,1],
    (label = "t_idx", range = 1:length(sol.t), startvalue = 1, update_while_dragging = false),
    )

    t_idx = sg.sliders[1].value

    t = @lift(sol.t[$t_idx])

    t_v = t[]

    ax = Axis(fig[1,1],ylabel=L"$\log_{10}\left(pF_{p}\right)$ $[\text{m}^{-3}]$",aspect=DataAspect())
    ax.limits = plot_limits

    if !isnothing(title)
        titlestr = @lift("Observer Flux at distance $ObserverDistance, with B-Field at an angle of β =$(β)π, at t=$(sol.t[$t_idx])")
        ax.title = titlestr
    end

    for θ in 1:length(ObserverAngles)

        flux_val = @lift(log10.(mp .* Fν[$t_idx,θ,:]))

        scatterlines!(ax,log10.(mp),flux_val,color=theme.palette.color[][mod(2*θ-1,7)+1],markersize=0.0,label= "θ=$(ObserverAngles[θ])π")

    end

    axislegend(ax)

    end # theme

    return fig

end