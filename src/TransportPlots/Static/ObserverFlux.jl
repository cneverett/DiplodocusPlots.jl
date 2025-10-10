function ObserverFlux(PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,ObserverAngles::Vector{Float64},ObserverDistance::Float64;R=nothing,Z=nothing)

    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    name_list = PhaseSpace.name_list

    c = getfield(DC,Symbol("c"))
    mEle = getfield(DC,Symbol("mEle"))

    if typeof(Space.space_coordinates) != Cylindrical
        error("ObserverFlux is only implemented for cylindrical coordinates.")
    end

    photon_index = findfirst(x->x=="Pho",name_list)

    β = Space.space_coordinates.β # local frame angle
    sβ, cβ = sincospi(β)

    d = ObserverDistance
    to_r = ObserverAngles

    
    if isnothing(R)
        R = Grids.xr[2]
    end 
    if isnothing(Z)
        Z = Grids.dz[1]
    end 

    ur = Grids.pyr_list[photon_index]
    mp = Grids.mpx_list[photon_index]
    dp = Grids.dpx_list[photon_index]
    du = Grids.dpy_list[photon_index]

    Fp = zeros(length(sol.t),length(ObserverAngles),Momentum.px_num_list[photon_index])

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
                
                Fp[t,θ_idx,px] = val_total
                Fp[t,θ_idx,px] *= R*Z/(4*pi*d^2) * c^2 * mEle * c

            end

        end

    end

    return Fp

end


"""
    ObserverFluxPlot(PhaseSpace,sol,time_idx,ObserverAngles,ObserverDistance;plot_limits=(nothing,nothing),theme=DiplodocusDark(),title=nothing)

"""
function ObserverFluxPlot(PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,time_idx::Int64,ObserverAngles::Vector{Float64},ObserverDistance::Float64;plot_limits=(nothing,nothing),theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,R=nothing,Z=nothing,fig_size=(3.25inch,2.4375inch))

    CairoMakie.activate!(inline=true)

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    c = getfield(DC,Symbol("c"))
    mEle = getfield(DC,Symbol("mEle"))
    ħ = getfield(DC,Symbol("ħ"))
    h = ħ*2pi

    β = Space.space_coordinates.β # local frame angle

    photon_index = findfirst(x->x=="Pho",name_list)

    ur = Grids.pyr_list[photon_index]
    pr = Grids.pxr_list[photon_index]
    mp = Grids.mpx_list[photon_index]
    dp = Grids.dpx_list[photon_index]

    if isnothing(R)
        R = Grids.xr[2]
    end 
    if isnothing(Z)
        Z = Grids.dz[1]
    end 

    Fp = ObserverFlux(PhaseSpace,sol,ObserverAngles,ObserverDistance;R=R,Z=Z)

    fig = Figure(size=fig_size)

    t = round(TimeUnits(sol.t[time_idx]),sigdigits=3)
    t_unit_string = TimeUnits()

    xlab = L"$\log_{10}\left(p [m_ec]\right)$"
    xlab_t = L"$\log_{10}\left(\nu [\text{Hz}]\right)$"
    ylab = L"$\log_{10}\left(pF_{p}\,[\text{J}\text{m}^{-2}\text{s}^{-1}]\right)$"
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax_t = Axis(fig[1,1],xlabel=xlab_t,xaxisposition=:top)
    #hidespines!(ax_t)
    linkyaxes!(ax,ax_t)
    #hidedecorations!(ax_t)

    if !isnothing(title)
        titlestr = L" β =%$(β)\pi, t=%$(t) %$t_unit_string"
        ax.title = titlestr
    end

    max_total = -Inf32

    for θ in 1:length(ObserverAngles)

        flux_val = log10.(mp .* Fp[time_idx,θ,:])
        max_f = maximum(x for x in flux_val if !isnan(x))
        max_total = max(max_total,max_f)

        scatterlines!(ax,log10.(mp),flux_val,color=theme.palette.color[][mod(2*θ-1,7)+1],markersize=0.0,label= L"θ=%$(ObserverAngles[θ])\pi")

    end

    if plot_limits == (nothing,nothing)
        xlims!(ax,(log10(pr[1]),log10(pr[end])))
        ylims!(ax,(max_total-9.0,max_total+1.0)) 
        xlims!(ax_t,(log10(pr[1]*mEle*c^2/h),log10(pr[end]*mEle*c^2/h)))
    else
        ax.limits = plot_limits
        plot_limits_x = plot_limits[1]
        plot_limits_y = plot_limits[2]
        ax_t.xlimits = log10.(10^plot_limits_x .* mEle*c^2/h)
    end

    axislegend(ax,position=:lt)

    return fig

    end # theme

end

function BFieldObserverPlot(sols::Vector{OutputStruct},PhaseSpaces::Vector{PhaseSpaceStruct},time_idx::Int64,ObserverAngle::Float64;ObserverDistance::Float64=1.0,theme=DiplodocusDark(),plot_limits=(nothing,nothing),TimeUnits::Function=CodeToCodeUnitsTime,title=nothing,R=nothing,Z=nothing)

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

        fig = Figure(size=(450,300)) # 6:3 ratio
        t = round(TimeUnits(sols[1].t[time_idx]),sigdigits=3)
        t_unit_string = TimeUnits()
        xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
        xlab_t = L"$\log_{10}\left(\nu\,[\text{Hz}]\right)$"
        ylab = L"$\log_{10}\left(\nu F_{\nu}\,[\text{J}\text{m}^{-2}\text{s}^{-1}]\right)$"
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect(),xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(2),xminorgridvisible=true,yminorgridvisible=true,xminorticksvisible=true,yminorticksvisible=true)
        ax_t = Axis(fig[1,1],xlabel=xlab_t,xaxisposition=:top,aspect=DataAspect())
        hidespines!(ax_t)
        hidexdecorations!(ax_t,ticklabels=false,ticks=false,label=false,minorgrid=false,minorticks=false)
        hideydecorations!(ax_t,ticklabels=false,ticks=false,minorgrid=false,minorticks=false)

        if !isnothing(title)
            titlestr = L"t=%$(t)\, %$t_unit_string, \theta_\text{Obs}=%$ObserverAngle \pi"
            ax.title = titlestr
        end

        max_total = -Inf32

        name_list = PhaseSpaces[1].name_list
        Space = PhaseSpaces[1].Space
        Momentum = PhaseSpaces[1].Momentum
        Grids = PhaseSpaces[1].Grids
        Time = PhaseSpaces[1].Time

        photon_index = findfirst(x->x=="Pho",name_list)

        ur = Grids.pyr_list[photon_index]
        pr = Grids.pxr_list[photon_index]
        mp = Grids.mpx_list[photon_index]
        dp = Grids.dpx_list[photon_index]

        c = getfield(DC,Symbol("c"))
        mEle = getfield(DC,Symbol("mEle"))
        ħ = getfield(DC,Symbol("ħ"))
        h = ħ*2pi

        if isnothing(R)
            R = Grids.xr[2]
        end 
        if isnothing(Z)
            Z = Grids.dz[1]
        end 

        legend_elements = []
        
        for sol in eachindex(sols)

            Fp = ObserverFlux(PhaseSpaces[sol],sols[sol],[ObserverAngle],ObserverDistance,R=R,Z=Z)
            pFp = log10.(mp .* Fp[time_idx,1,:])
            max_f = maximum(x for x in pFp if !isnan(x))
            max_total = max(max_total,max_f)

            println(max_total)

            if sol != length(sols)
                scatterlines!(ax,log10.(mp),pFp,color=theme.palette.color[][mod(sol,7)+1],markersize=0.0)
                push!(legend_elements,LineElement(color = theme.palette.color[][mod(sol,7)+1], linestyle = :solid,linewidth = 2.0))
            else
                scatterlines!(ax,log10.(mp),pFp,color=theme.textcolor[],markersize=0.0,linestyle=(:dash,:dense))
                push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = (:dash,:dense),linewidth = 2.0))

                scatterlines!(ax_t,log10.(mp .* (mEle*c^2/h)),pFp,color=:transparent,markersize=0.0,linestyle=(:dash,:dense))
            end  

        end

        if plot_limits == (nothing,nothing)
            xlims!(ax,(log10(pr[1]),log10(pr[end])))
            ylims!(ax,(max_total-9.0,max_total+1.0)) 
            xlims!(ax_t,(log10(pr[1]*mEle*c^2/h),log10(pr[end]*mEle*c^2/h)))
        else
            ax.limits = plot_limits
            plot_limits_x = plot_limits[1]
            plot_limits_y = plot_limits[2]
            ax_t.limits = (log10.(10 .^ plot_limits_x .* (mEle*c^2/h)),plot_limits_y)
        end

        #line_labels=[L"\theta_\text{B}=0.0\pi",L"\theta_\text{B}=0.1\pi",L"\theta_\text{B}=0.2\pi",L"\theta_\text{B}=0.3\pi",L"\theta_\text{B}=0.4\pi",L"\theta_\text{B}=0.5\pi",L"\text{Iso}"]

        line_labels=[L"\theta_B=0",L"\theta_B=\pi/6",L"\theta_B=\pi/3",L"\theta_B=\pi/2",L"\text{Iso}"]

        axislegend(ax,legend_elements,line_labels,position = :lt)

        return fig

    end # with theme

end

