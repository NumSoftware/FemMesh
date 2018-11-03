using PyCall
function cplot(X, Y, filename=""; xlabel="\$x\$", ylabel="\$y\$", lw=0.7, ms=2, marker=nothing, color="", legend=[], legendloc="best",
               xbins=6, ybins=6, grid=false, figsize=(3,2), legendexpand=false)

    @eval import PyPlot:plt, matplotlib, figure

    labels = Dict( 
                  "sxx"=>raw"$\sigma_{xx}$",
                  "syy"=>raw"$\sigma_{yy}$",
                  "szz"=>raw"$\sigma_{zz}$",
                  "sxy"=>raw"$\sigma_{xy}$",
                  "syz"=>raw"$\sigma_{yz}$",
                  "sxz"=>raw"$\sigma_{xz}$",
                  "exx"=>raw"$\varepsilon_{xx}$",
                  "eyy"=>raw"$\varepsilon_{yy}$",
                  "ezz"=>raw"$\varepsilon_{zz}$",
                  "exy"=>raw"$\varepsilon_{xy}$",
                  "eyz"=>raw"$\varepsilon_{yz}$",
                  "exz"=>raw"$\varepsilon_{xz}$",
                 )
    xlabel = get(labels, xlabel, xlabel)
    ylabel = get(labels, ylabel, ylabel)

    # Configure plot
    plt[:close]("all")

    plt[:rc]("mathtext", fontset="cm")
    plt[:rc]("font", family="STIXGeneral", size=7)
    plt[:rc]("lines", scale_dashes=true)

    if filename!=""
        plt[:ioff]()
        plt[:rc]("xtick", labelsize=7)
        plt[:rc]("ytick", labelsize=7)
        plt[:rc]("lines", lw=lw)
        plt[:rc]("lines", markersize=2)
        plt[:rc]("axes" , linewidth=0.5)
        plt[:rc]("figure", figsize=(3, 2))
        plt[:rc]("legend", fontsize=7)
    end

    # Check the size of Y
    n = length(X)
    @assert n==size(Y,1)
    if size(Y,2)==1
        Y = reshape(n,1)
    end
    m = size(Y,2)

    # Fix markers
    if typeof(marker)==String || marker==nothing
        marker = repeat([marker], m)
    end
    @assert length(marker)==m

    # Fix markers size
    if typeof(ms) <: Real
        ms = repeat([ms], m)
    end
    @assert length(ms)==m

    # Fix line colors
    if typeof(color)==String
        if color==""
            color = nothing
        end
        color = repeat([color], m)
    end
    @assert length(color)==m

    # Fix legend labels
    haslegend = true
    if length(legend)==0
        haslegend = false
        legend = repeat([""],m)
    end
    @assert length(legend)==m


    # Plot curves
    for i=1:m
        plt[:plot](X, Y[:,i], marker=marker[i], ms=ms[i], color=color[i], lw=lw, label=legend[i])
    end

    plt[:grid](grid, ls="dotted", lw=0.3)
    ax = plt[:axes]()
    plt[:locator_params](axis="x", nbins=xbins)
    plt[:locator_params](axis="y", nbins=ybins)

    # Print axes labels
    plt[:xlabel](xlabel)
    plt[:ylabel](ylabel)

    # plot legend
    if haslegend
        mode = legendexpand ? "expand" : nothing

        if legendloc=="top"
            leg = plt[:legend](loc="lower left", bbox_to_anchor=(-0.02, 1.01, 1.04, 0.2), edgecolor="k", ncol=m, mode=mode)
        elseif legendloc=="right"
            leg = plt[:legend](loc="upper left", bbox_to_anchor=(1.01, 1), edgecolor="k")
        elseif legendloc=="bottom"
            leg = plt[:legend](loc="upper left", bbox_to_anchor=(-0.02, -0.02, 1.04, -0.2), edgecolor="k", ncol=m, mode=mode)
        else
            leg = plt[:legend](loc=legendloc, edgecolor="k")
        end

        frame = leg[:get_frame]()
        frame[:set_linewidth](0.5)
    end

    # show or save plot
    if filename==""
        plt[:show]()
    else
        plt[:savefig](filename, bbox_inches="tight", pad_inches=0.01, format="pdf")
    end
end
