


function draw_error_band(ax, WList, lower, upper, col) # defined as a function for uniform transparency acrross plots

     band!(ax, WList, lower, upper;
              color = (col, 0.1))#, strokewidth = 0)
end