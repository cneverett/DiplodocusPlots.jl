module DiplodocusPlots

    export DiplodocusDark, DiplodocusLight
    export MomentumDistributionPlot

    using CairoMakie
    using GLMakie

    import DiplodocusCollisions as DC
    import DiplodocusTransport as DT

    include("Themes.jl")
    include("CollisionPlots/Interactive/AnalyticalSpectra.jl/InverseCompton.jl")
    include("CollisionPlots/Interactive/GainLossAnalysis.jl")

    # Transport Plots
    # Static Plots
    include("TransportPlots/Static/Distribution.jl")
    #Interactive Plots
    include("TransportPlots/Interactive/ObserverFlux.jl")

end
