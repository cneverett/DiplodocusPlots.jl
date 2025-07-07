function DiplodocusDark()
    # inspired from Makie's theme_black but with latex fonts and some modifications
    return Theme(
        backgroundcolor=:transparent,
        textcolor=:white,
        linecolor=:white,
        #palette = Makie.generate_default_palette(:black),
        fonts = Attributes(
            :bold => Makie.texfont(:bold),
            :bolditalic => Makie.texfont(:bolditalic),
            :italic => Makie.texfont(:italic),
            :regular => Makie.texfont(:regular)
        ),
        Figure = (
            backgroundcolor = :transparent,
            size=(1200,1200),
            fontsize = 20.0f0
        ),
        Axis = (
            backgroundcolor=:transparent,
            bottomspinecolor = :white,
            topspinecolor = :white,
            leftspinecolor = :white,
            rightspinecolor = :white,
            xgridcolor = :grey45,
            ygridcolor = :grey45,
            xtickcolor = :white,
            ytickcolor = :white,
            xtickwidth = 2.0,
            ytickwidth = 2.0,
            # modification from theme_black
            titlesize = 20.0f0,
            xlabelsize = 20.0f0,
            ylabelsize = 20.0f0,
            titlecolor= :white,
            ylabelcolor= :white,
            xlabelcolor= :white,
            spinewidth = 2.0,
            xgridwidth = 1.0,
            ygridwidth = 1.0,
            xminorgridcolor=:white,
            yminorgridcolor=:white,
            xminortickcolor=:white,
            yminortickcolor=:white,
            xminorgridwidth = 2.0,
            yminorgridwidth = 2.0,
            xticklabelcolor=:white,
            yticklabelcolor=:white,
        ),
        Legend = (
            framecolor = :white,
            backgroundcolor = :transparent,
            framewidth = 2.0,
        ),
        Colorbar = (
            tickcolor = :white,
            spinecolor = :white,
            topspinecolor = :white,
            bottomspinecolor = :white,
            leftspinecolor = :white,
            rightspinecolor = :white,
            # modifications from theme_black
            spinewidth= 2.0,
            labelsize = 20.0f0
        ),
        ScatterLines = (
            linewidth = 3.0f0,
        ),
        PolarAxis = (
            backgroundcolor = :transparent,
            rticksvisible = true,
            rticksmirrored = true,
            rtickcolor = :white,
            rtickwidth = 2.0,
            gridcolor = :grey45,
            gridwidth = 1.0,
            tickcolor = :white,
            tickwidth = 2.0,
            labelsize = 20.0f0,
            labelcolor = :white,
            spinewidth = 2.0,
            spinecolor = :white,
            clip = false,
        )
    )
end

function DiplodocusLight()
    # slight modifications from the default theme but with Latex fonts and some modifications
    return Theme(
        backgroundcolor=:transparent,
        textcolor=:black,
        linecolor=:black,
        #palette = Makie.generate_default_palette(:white),
        fonts = Attributes(
            :bold => Makie.texfont(:bold),
            :bolditalic => Makie.texfont(:bolditalic),
            :italic => Makie.texfont(:italic),
            :regular => Makie.texfont(:regular)
        ),
        Figure = (
            backgroundcolor = :transparent,
            size=(1200,1200),
            fontsize = 20.0f0
        ),
        Axis = (
            backgroundcolor=:transparent,
            titlesize = 20.0f0,
            xlabelsize = 20.0f0,
            ylabelsize = 20.0f0,
            titlecolor= :black,
            ylabelcolor= :black,
            xlabelcolor= :black,
            spinewidth = 2.0,
            xgridwidth = 1.0,
            ygridwidth = 1.0,
            xgridcolor = :grey45,
            ygridcolor = :grey45,
            xminorgridcolor=:black,
            yminorgridcolor=:black,
            xminorgridwidth = 2.0,
            yminorgridwidth = 2.0,
            xtickcolor=:black,
            ytickcolor=:black,
            xtickwidth = 2.0,
            ytickwidth = 2.0,
            xticklabelcolor=:black,
            yticklabelcolor=:black,
        ),
        Colorbar = (
            backgroundcolor = :transparent,
            spinewidth=1.5,
            labelsize = 20.0f0
        ),
        Legend = (
            backgroundcolor = :transparent,
            framecolor = :black,
            framewidth = 2.0,
        ),
        ScatterLines = (
            linewidth = 3.0f0,
        ),
        PolarAxis = (
            rticksvisible = true,
            rticksmirrored = true,
            rtickcolor = :black,
            rtickwidth = 2.0,
            gridcolor = :grey45,
            gridwidth = 1.0,
            tickcolor = :black,
            tickwidth = 2.0,
            labelsize = 20.0f0,
            labelcolor = :black,
            spinewidth = 2.0,
            spinecolor = :black,
            clip = false,
        )
    )
end