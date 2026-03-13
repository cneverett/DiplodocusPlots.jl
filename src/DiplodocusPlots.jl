module DiplodocusPlots

    export DiplodocusDark, DiplodocusLight
    export Static, Animated, Interactive
    export MomentumDistributionPlot0D, MomentumAndPolarAngleDistributionPlot0D
    # export AngleDistributionPlot, MomentumComboAnimation, MomentumAndAzimuthalAngleDistributionPlot, AzimuthalAngleDistributionPlot
    export TimeScalePlot
    export InteractiveBinaryGainLossPlot
    export NumberDensityPlot0D
    export EnergyDensityPlot0D
    export IsThermalAndIsotropicPlot0D
    #export IsThermalPlot, IsIsotropicPlot, 
    export InteractiveEmissionGainLossPlot
    export CodeToSIUnitsTime, CodeToCodeUnitsTime, SIToCodeUnitsTime, SyncToCodeUnitsTime, CodeToSyncUnitsTime
    export ObserverFluxPlot

    using CairoMakie
    using GLMakie

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
    include("TransportPlots/Static/ObserverFlux.jl")
    include("TransportPlots/Static/TimeScale.jl")

    # 0D (space) plots
    include("TransportPlots/Plots_0D/Distribution0D.jl")
    include("TransportPlots/Plots_0D/NumberDensity0D.jl")
    include("TransportPlots/Plots_0D/EnergyDensity0D.jl")
    include("TransportPlots/Plots_0D/ThermalAndIsotropic0D.jl")

    # 1D (space) plots
    include("TransportPlots/Static/Distribution1D.jl")
    include("TransportPlots/Static/NumberDensity1D.jl")

    # 2D (space) plots 
    include("TransportPlots/Static/NumberDensity2D.jl")
    include("TransportPlots/Static/BulkGamma2D.jl")

end
