module DiplodocusPlots

    export DiplodocusDark, DiplodocusLight
    export MomentumDistributionPlot
    export InteractiveBinaryGainLossPlot
    export FracNumberDensityPlot, NumberDensityPlot
    export FracEnergyDensityPlot, EnergyDensityPlot
    export IsThermalPlot, IsIsotropicPlot, IsThermalAndIsotropicPlot

    using CairoMakie
    using GLMakie

    import DiplodocusCollisions as DC
    using DiplodocusTransport

    include("Themes.jl")
    include("CollisionPlots/Interactive/AnalyticalSpectra.jl/InverseCompton.jl")
    include("CollisionPlots/Interactive/GainLossAnalysis.jl")

    # Transport Plots
    # Static Plots
    include("TransportPlots/Static/Distribution.jl")
    include("TransportPlots/Static/NumberDensity.jl")
    include("TransportPlots/Static/EnergyDensity.jl")
    include("TransportPlots/Static/ThermalAndIsotropic.jl")
    #Interactive Plots
    include("TransportPlots/Interactive/ObserverFlux.jl")

end
