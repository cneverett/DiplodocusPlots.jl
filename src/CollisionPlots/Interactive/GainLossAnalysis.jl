function InteractiveGainLossPlot(Output::Tuple)

    GLMakie.activate!(inline=false)

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Output[1]
    #SMatrix3 = Output[2]
    #SMatrix4 = Output[3]
    #TMatrix1 = Output[4]
    #TMatrix2 = Output[5]
    SMatrix3 = dropdims(sum(Output[2],dims=(2,3,5,6,8,9)),dims=(2,3,5,6,8,9));
    SMatrix4 = dropdims(sum(Output[3],dims=(2,3,5,6,8,9)),dims=(2,3,5,6,8,9));
    TMatrix1 = dropdims(sum(Output[4],dims=(2,3,5,6)),dims=(2,3,5,6));
    TMatrix2 = dropdims(sum(Output[5],dims=(2,3,5,6)),dims=(2,3,5,6));

    p2_r = DC.bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = DC.deltaVector(p2_r);
    p2_m = DC.meanVector(p2_r);
    p1_r = DC.bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = DC.deltaVector(p1_r);
    p1_m = DC.meanVector(p1_r);


    fig = Figure(size = (960,540))
    ax1 = Axis(fig[1, 1],xlabel=L"$\log_{10} \text{momentum}[m_\text{Ele}c]$",ylabel= L"$\log_{10} Sp [\text{some units}]$")

    sg = SliderGrid(fig[2,1],
    (label = "p1", range = 1:length(p1_m), startvalue = 25, update_while_dragging = false),
    (label = "p2", range = 1:length(p2_m), startvalue = 25, update_while_dragging = false)
    )

    p1_idx = sg.sliders[1].value
    p2_idx = sg.sliders[2].value

    p1_val = @lift(log10.(p1_m[$p1_idx]))
    p2_val = @lift(log10.(p2_m[$p2_idx]))

    vlines!(ax1,p1_val,color=:black,linestyle=:dash)
    vlines!(ax1,p2_val,color=:black,linestyle=:dot)

    # analytic kernels
    p3_val_true = @lift(log10.(replace(IC_kernel.(p2_m[$p2_idx],p1_m[$p1_idx],p2_m[$p2_idx]+sqrt(p1_m[$p1_idx]^2+1e0^2).- sqrt.(p1_m.^2 .+1e0^2)),0.0 => NaN) .* (p1_m ./ sqrt.(p1_m.^2 .+ 1e0^2)).* p1_d))
    p4_val_true = @lift(log10.(replace(IC_kernel.(p2_m[$p2_idx],p1_m[$p1_idx],p2_m),0.0 => NaN) .* p2_d))
    scatterlines!(ax1,log10.(p1_m),p3_val_true,color=:green,strokewidth=1.0,markersize=0.0)
    scatterlines!(ax1,log10.(p2_m),p4_val_true,markersize=0.0,color=:red, strokewidth=1.0);

    # Loss Matrix   
    p1_loss = @lift(log10.(TMatrix1[$p1_idx,$p2_idx]))
    p2_loss = @lift(log10.(TMatrix2[$p2_idx,$p1_idx]))
    hlines!(ax1,p1_loss)
    hlines!(ax1,p2_loss)

    # Gain Matrix
    p3_gain = @lift(log10.(SMatrix3[1:end-1,$p1_idx,$p2_idx])) # ignore overflow
    p4_gain = @lift(log10.(SMatrix4[1:end-1,$p1_idx,$p2_idx])) # ignore overflow
    scatter!(ax1,log10.(p1_m),p3_gain)
    scatter!(ax1,log10.(p2_m),p4_gain)

    ax1.aspect = DataAspect()

    return fig

end