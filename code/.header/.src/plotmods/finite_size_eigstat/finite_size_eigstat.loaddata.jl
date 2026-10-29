

function plotdata_load(file_plotdata::HDF5.File, quantity_name::String)
    
    
    file_attrs=HDF5.attributes(file_plotdata)

    LList=read(file_attrs["[Parameters] LList"])
    
    Ltot=length(LList)

    quantity=Array{Array{Float64,1},1}(undef, Ltot)

    for i in 1:Ltot
        quantity[i]=read(file_plotdata, "L$(LList[i])/$quantity_name")
    end

    return quantity
end




function plotdata_load(file_plotdata::HDF5.File, path::String ,quantity_name::String)
    
    func="quantile80"
    
    file_attrs=HDF5.attributes(file_plotdata)

    LList=read(file_attrs["[Parameters] LList"])
    WList=read(file_attrs["[Parameters] WList"])

    Ltot=length(LList)

    quantity=Array{Array{Float64,1},1}(undef, Ltot)

    for i in 1:Ltot
        quantity[i]=read(file_plotdata, "L$(LList[i])/$path/$quantity_name")
    end

    return quantity
end



function plotdata_transpose(quantity::Array{Array{Float64,1},1})


    Ltot=length(quantity)
    Wtot=length(quantity[1])

    quant_transposed=Array{Array{Float64,1},1}(undef, Wtot)

    for j in 1:Wtot
        quant_transposed[j]=zeros(Ltot)
    end

    for i in 1:Ltot
        for j in 1:Wtot
            quant_transposed[j][i]=quantity[i][j]
        end
    end
    return quant_transposed
end