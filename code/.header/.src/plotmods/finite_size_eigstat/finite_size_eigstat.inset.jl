using CairoMakie
using ColorSchemes

#------------------------------------------------------------------------------
# First we define the plots as insets
# Then we will use the multiple dispatch to create a function with default axis
# The other options should be automatically followed through;
#------------------------------------------------------------------------------


#----------------------------------------------------------------
#  Inset:
#----------------------------------------------------------------

#  Without Error Bands:


function finite_size_eigstat_in(ax, LList, WList, quantity, cmap)
    ############################################
    # finite size scaling
    ############################################

    WList = reverse(WList)
    
    # define a consistent color palette
    pal = cgrad(cmap, 6 * length(LList), categorical = true)

    for (i, L) in enumerate(LList)
        col = pal[6 * i - 3 ]


    
        # mean line on top
        lines!(ax, WList[1:length(quantity[i])], quantity[i][:];
               color = col, linewidth = 3, label = "L=$(L)")
        scatter!(ax, WList[1:length(quantity[i])], quantity[i][:];
               color = col)
   end

    axislegend(ax)
    return fig
end



# With Error Bands:


function finite_size_eigstat_in(ax, LList, WList, quantity, lower, upper, cmap)
    ############################################
    # finite size scaling
    ############################################

    #cmap=cgrad(cmap, rev=true)

    #xlims!(ax, (WList[1], WList[end]))
    #ylims!(ax, (1e-3, 1e0))
    WList = reverse(WList)
        #ax.xreversed = true
    # define a consistent color palette
    pal = cgrad(cmap, 2 * length(LList), categorical = true)


    for (i, L) in enumerate(LList)
        col = pal[2 * length(LList) - 2 * i + 1]



        # confidence band first (semi-transparent)
        band!(ax, WList[1:length(lower[i][:])], lower[i][:], upper[i][:];
              color = (col, 0.15))#, strokewidth = 0)      

         # mean line on top
        lines!(ax, WList[1:length(quantity[i])], quantity[i][:];
               color = col, linewidth = 3, label = "L=$(L)")
        scatter!(ax, WList[1:length(quantity[i])], quantity[i][:];
               color = col)
    end

    axislegend(ax)
    return fig
end



#-------------------------------------------------------
# Default ColorSchemes
#-------------------------------------------------------

# Insets

        function finite_size_eigstat_in(ax, LList, WList, quantity)
            finite_size_eigstat_in(ax, LList, WList, quantity, :viridis)
        end

        function finite_size_eigstat_in(ax, LList, WList, quantity, lower, upper)
            finite_size_eigstat_in(ax, LList, WList, quantity, lower, upper, :viridis)
        end





#-------------------------------------------------------
# Specific Cases
#-------------------------------------------------------


        
function finite_size_eigstat_lsr(ax, LList, WList, quantity, lower, upper, cmap)
    ############################################
    # finite size scaling
    ############################################

    #cmap=cgrad(cmap, rev=true)

    xlims!(ax, (WList[1], WList[end]))
    ylims!(ax, (0.35, 0.65))

     WList = reverse(WList)
     ax.xreversed = true
    # define a consistent color palette
    pal = cgrad(cmap, 6 * length(LList), categorical = true)
    
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


        
function finite_size_eigstat_hist(ax, LList, WList, quantity, cmap)
    ############################################
    # finite size scaling
    ############################################

    #cmap=cgrad(cmap, rev=true)

    xlims!(ax, (WList[1], WList[end]))
    #ylims!(ax, (0.35, 0.65))

     WList = reverse(WList)
     ax.xreversed = true
    # define a consistent color palette
    pal = cgrad(cmap, 2 * length(LList), categorical = true)
    
    for (i, L) in enumerate(LList)
        col = pal[2 * length(LList) - 2 * i + 1]


        # mean line on top
        stairs!(ax, WList[1:length(quantity[i][:])], quantity[i][:];
               color = col, linewidth = 2, label = "L=$(L)")
    end


    axislegend(ax)
    return fig
end
