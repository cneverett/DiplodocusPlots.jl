inch = 96
pt = 4/3
cm = inch / 2.54 

#= 
Latex double column width is 236pt which is 3.29 inches or 8.36cm

Deault figure size will be designed for single column with a width of 3 inches or 7.62cm and a aspect ratio of 4:3, therefore height of 2.25 inches or 5.72cm

=#

function DiplodocusDark()
    # inspired from Makie's theme_black but with latex fonts and some modifications
    return Theme(
        backgroundcolor=:transparent,
        textcolor=:white,
        linecolor=:white,
        size=(3inch,2.25inch),
        fontsize = 9pt,
        #palette = Makie.generate_default_palette(:black),
        fonts = Attributes(
            :bold => Makie.texfont(:bold),
            :bolditalic => Makie.texfont(:bolditalic),
            :italic => Makie.texfont(:italic),
            :regular => Makie.texfont(:regular)
        ),
        Figure = (
            backgroundcolor = :transparent,
            fontsize = 9pt
        ),
        Axis = (
            fontsize = 9pt,
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
            titlesize = 11pt,
            #xlabelsize = 10pt, # defults if :fontsize
            #ylabelsize = 10pt, # defults if :fontsize
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
            xticklabelsize = 9pt,
            yticklabelsize = 9pt,
        ),
        Legend = (
            fontsize = 9pt,
            labelsize = 9pt,
            framecolor = :white,
            backgroundcolor = :transparent,
            framewidth = 2.0,
        ),
        Colorbar = (
            fontsize = 9pt,
            labelsize = 9pt,
            ticklabelsize = 9pt,
            tickcolor = :white,
            spinecolor = :white,
            topspinecolor = :white,
            bottomspinecolor = :white,
            leftspinecolor = :white,
            rightspinecolor = :white,
            # modifications from theme_black
            spinewidth= 2.0,
        ),
        ScatterLines = (
            linewidth = 3.0,
        ),
        PolarAxis = (
            fontsize = 9pt,
            labelsize = 9pt,
            rticklabelsize = 9pt,
            backgroundcolor = :transparent,
            rticksvisible = true,
            rticksmirrored = true,
            rtickcolor = :white,
            rtickwidth = 2.0,
            gridcolor = :grey45,
            gridwidth = 1.0,
            tickcolor = :white,
            tickwidth = 2.0,
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
        size=(3inch,2.25inch),
        fontsize = 9pt,
        #palette = Makie.generate_default_palette(:white),
        fonts = Attributes(
            :bold => Makie.texfont(:bold),
            :bolditalic => Makie.texfont(:bolditalic),
            :italic => Makie.texfont(:italic),
            :regular => Makie.texfont(:regular)
        ),
        Figure = (
            backgroundcolor = :transparent,
            fontsize = 9pt
        ),
        Axis = (
            fontsize = 9pt,
            backgroundcolor=:transparent,
            titlesize = 11pt,
            xlabelsize = 9pt, # defults if :fontsize
            ylabelsize = 9pt, # defults if :fontsize
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
            xticklabelsize = 9pt,
            yticklabelsize = 9pt,
        ),
        Colorbar = (
            fontsize = 9pt,
            labelsize = 9pt,
            ticklabelsize = 9pt,
            backgroundcolor = :transparent,
            spinewidth=1.5,
        ),
        Legend = (
            fontsize = 9pt,
            labelsize = 9pt,
            backgroundcolor = :transparent,
            framecolor = :black,
            framewidth = 2.0,
        ),
        ScatterLines = (
            linewidth = 3.0f0,
        ),
        PolarAxis = (
            fontsize = 9pt,
            labelsize = 9pt,
            rticklabelsize = 9pt,
            thetaticklabelsize = 5pt,
            rticksvisible = true,
            rticksmirrored = true,
            rtickcolor = :black,
            rtickwidth = 2.0,
            gridcolor = :grey45,
            gridwidth = 1.0,
            tickcolor = :black,
            tickwidth = 2.0,
            labelcolor = :black,
            spinewidth = 2.0,
            spinecolor = :black,
            clip = false,
        )
    )
end