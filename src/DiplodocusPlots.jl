module DiplodocusPlots

    export DiplodocusDark, DiplodocusLight
    export Static, Animated, Interactive
    export MomentumDistributionPlot, AngleDistributionPlot, MomentumAndPolarAngleDistributionPlot, MomentumComboAnimation, MomentumAndAzimuthalAngleDistributionPlot, AzimuthalAngleDistributionPlot
    export InteractiveBinaryGainLossPlot
    export FracNumberDensityPlot, NumberDensityPlot
    export FracEnergyDensityPlot, EnergyDensityPlot
    export IsThermalPlot, IsIsotropicPlot, IsThermalAndIsotropicPlot
    export InteractiveEmissionGainLossPlot
    export CodeToSIUnitsTime, CodeToCodeUnitsTime, SIToCodeUnitsTime, SyncToCodeUnitsTime, CodeToSyncUnitsTime
    export ObserverFluxPlot

    #using Makie
    using CairoMakie
    using GLMakie
    #using Interpolations
    #using ScatteredInterpolation
    #using DIVAnd

    import DiplodocusCollisions as DC
    using DiplodocusTransport

    include("Themes.jl")
    include("PlotTypes.jl")
    include("TimeUnits.jl")
    include("CollisionPlots/Interactive/AnalyticalSpectra.jl/InverseCompton.jl")
    include("CollisionPlots/Interactive/GainLossAnalysis.jl")

    # Transport Plots
    #Interactive Plots
    #include("TransportPlots/Interactive/ObserverFlux.jl")
    # Static Plots
    include("TransportPlots/Static/Distribution.jl")
    include("TransportPlots/Static/NumberDensity.jl")
    include("TransportPlots/Static/EnergyDensity.jl")
    include("TransportPlots/Static/ThermalAndIsotropic.jl")
    include("TransportPlots/Static/ObserverFlux.jl")
    include("TransportPlots/Static/TimeScale.jl")
    

end
