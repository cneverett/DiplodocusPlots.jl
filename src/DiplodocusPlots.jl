module DiplodocusPlots

    #using CairoMakie
    using GLMakie

    import DiplodocusCollisions as DC
    import DiplodocusTransport as DT

    include("CollisionPlots/Interactive/AnalyticalSpectra.jl/InverseCompton.jl")
    include("CollisionPlots/Interactive/GainLossAnalysis.jl")

    include("TransportPlots/ObserverFlux.jl")

end
