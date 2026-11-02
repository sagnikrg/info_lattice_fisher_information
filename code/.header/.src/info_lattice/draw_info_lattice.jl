using CairoMakie
using GeometryBasics#: Point2f
using ColorSchemes
using HDF5




function inset_info_lattice(infolattice, ax_inset)
        L=length(infolattice[1])

        info_lattice_sum=fill(0.0,L)
        for i in 1:L
            info_lattice_sum[i]=sum(infolattice[i])
        end

        file = h5open("C:\\Users\\Ritam\\OneDrive\\Documents\\Codespaces\\Codespace_Th\\.header\\.header\\.src\\info_lattice\\info_lattice_inset_data.hdf5", "r")
        info_lattice_sum_GUE=file["GUE/info_lattice_sum/L$(L)"][:] 
        close(file)

        xlims!(ax_inset, 1, L)
        lines!(ax_inset, 1:L , info_lattice_sum_GUE, color=reds[12], linewidth=1, linestyle=:dash, label="GUE")
        lines!(ax_inset, 1:L , (1:L).^(-2).*L./sum((1:L).^(-2)), color=viridis[42], linewidth=1, linestyle=:dash, label="conformal")
        lines!(ax_inset, 1:L , info_lattice_sum, color=viridis[24], linewidth=2, label="data")

    
end



function draw_info_lattice(infolattice, cmap)
    fig = Figure(size = (800, 600), fontsize = 12)
    ax = Axis(fig[1, 1], aspect = DataAspect(),
            xlabel=L"i", ylabel=L"l", 
            xticks = (1.5:0.5:(length(infolattice[1])+0.5), string.(1:0.5:(length(infolattice[1])))), 
            yticks = ((1:(length(infolattice[1]))).*(sqrt(3)/2), string.(1:(length(infolattice[1]))))
            )

    a = Point2f(1, 0)
    b = Point2f(0.5, sqrt(3)/2)

    xlims!(ax, 0, length(infolattice[1])+2)
    ylims!(ax, 0, length(infolattice))

    radius = 0.55  # Size of the hexagons
    hex_vertices = [Point2f(radius * cos(θ), radius * sin(θ)) for θ in π/6:π/3:2π]

    for (nb, line) in enumerate(infolattice)
        for (na, val) in enumerate(line)
            coord = na * a + nb * b
            poly = Polygon([v + coord for v in hex_vertices])
            p=val
           # Add rounded effect using stroke
           # poly!(ax, poly,  color = (viridis[36]*(p)+magma[1]*(1-p)),strokewidth = 1.5, strokecolor = :black)
           poly!(ax, poly,  color = val, colormap = cmap, colorrange = (1, 0), strokewidth = 1.5, strokecolor = :black)
        end
           
                  
    end


        
        ax_inset = Axis(fig[1, 1],
        width=Relative(0.28),
        height=Relative(0.22),
        halign=0.12,
        valign=0.95,
        title="Total bits of information per level",
        xlabel=L"l",
        ylabel=L"I(l)",
        xticks = ((1:(length(infolattice[1]))), string.(1:(length(infolattice[1])))))
        
    inset_info_lattice(infolattice, ax_inset)

    axislegend(ax_inset, position=:rt, framevisible=false)
    Colorbar(fig[1,2] , colormap=cmap,                      #Colorbar
                 ticks = 0.0:0.1:1.0)            #Colorbar ticks
    
    fig
    
    return fig
end



;
;

function draw_info_lattice(infolattice)
    cmap=Reverse(:navia)
    fig = draw_info_lattice(infolattice, cmap) 
    return fig
end

;
;