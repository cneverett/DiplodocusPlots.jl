module DiplodocusPlots

    #using CairoMakie
    using GLMakie

    import DiplodocusCollisions as DC
    import DiplodocusTransport as DT

    include("..\\..\\BoltzmannCollisionIntegral\\CompTests.jl") # remove later
    include("CollisionPlots/Interactive/GainLossAnalysis.jl")

    include("TransportPlots/ObserverFlux.jl")

end
