# Basic plotting blocks for finite size scaling



# without confidence bands


function finite_size_eigstat(LList, WList, quantity, cmap)
    ############################################
    # finite size scaling
    ############################################

    aspect_ratio = 4/3
    fig_width = 800
    fig_height = fig_width / aspect_ratio

    fig = Figure(resolution = (fig_width, fig_height))
    ax = Axis(fig[1, 1])
    
    ax.xreversed = true
    xlims!(ax, (WList[1], WList[end]))
    # define a consistent color palette
    pal = cgrad(cmap, 2 * length(LList), categorical = true)

    for (i, L) in enumerate(LList)
        col = pal[2 * length(LList) - 2 * i + 1]

   
        # mean line on top
        lines!(ax, WList[1:length(quantity[i])], quantity[i][:];
               color = col, linewidth = 3, label = "L=$(L)")
    end

    axislegend(ax)
    return fig
end



# with confidence bands


function finite_size_eigstat(LList, WList, quantity, lower, upper, cmap)

 
    ############################################
    # finite size scaling
    ############################################

    aspect_ratio = 4/3
    fig_width = 800
    fig_height = fig_width / aspect_ratio

    fig = Figure(resolution = (fig_width, fig_height))
    ax = Axis(fig[1, 1])

    ax.xreversed = true
    xlims!(ax, (WList[1], WList[end]))
    #ylims!(ax, (1e-3, 1e0))
    # define a consistent color palette
    pal = cgrad(cmap, 2 * length(LList), categorical = true)


    for (i, L) in enumerate(LList)
        col = pal[2 * length(LList) - 2 * i + 1]



        # confidence band first (semi-transparent)
        band!(ax, WList[1:length(lower[i][:])], lower[i][:], upper[i][:];
              color = (col, 0.15))#, strokewidth = 0)      

        # mean line on top
        lines!(ax, WList[1:length(quantity[i][:])], quantity[i][:];
               color = col, linewidth = 3, label = "L=$(L)")
    end

    axislegend(ax)
    return fig
end


# with two confidence bands

function finite_size_eigstat(LList, WList, quantity, lower_1, upper_1, lower_2, upper_2 , cmap)
    ############################################
    # finite size scaling
    ############################################

    aspect_ratio = 4/3
    fig_width = 800
    fig_height = fig_width / aspect_ratio

    fig = Figure(resolution = (fig_width, fig_height))
    
    ax = Axis(fig[1, 1], hidedecorations=:all)

   
    #ax.xreversed = true
    xlims!(ax, (WList[1], WList[end]))
    ylims!(ax, )
    # define a consistent color palette
    pal = cgrad(cmap, 2 * length(LList), categorical = true)

    for (i, L) in enumerate(LList)
        col = pal[2 * length(LList) - 2 * i + 1]

                   
        # confidence band second (semi-transparent)
        band!(ax, WList[1:length(lower_2[i])], lower_2[i][:], upper_2[i][:];
              color = (col, 0.1))#, strokewidth = 0)      


        # confidence band first (semi-transparent)
        band!(ax, WList[1:length(lower_1[i])], lower_1[i][:], upper_1[i][:];
              color = (col, 0.1))#, strokewidth = 0)            
        # mean line on top
        lines!(ax, WList[1:length(quantity[i])], quantity[i][:];
               color = col, linewidth = 3, label = "L=$(L)")
    end

    axislegend(ax)
    return fig
end


#########################################################################
# Wrappers for Default ColorSchemes
#########################################################################


function finite_size_eigstat(LList, WList, quantity)
    finite_size_eigstat(LList, WList, quantity,:viridis)
end

function finite_size_eigstat(LList, WList, quantity, lower, upper )
    finite_size_eigstat(LList, WList, quantity, lower, upper , :viridis)
end

function finite_size_eigstat(LList, WList, quantity, lower_1, upper_1, lower_2, upper_2 )
    finite_size_eigstat(LList, WList, quantity, lower_1, upper_1, lower_2, upper_2 , :viridis)
end