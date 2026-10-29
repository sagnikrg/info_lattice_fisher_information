### functions for automated extraction of environment details from printed metadata




function read_hdf5_attributes(file::HDF5.File)
    attrs = HDF5.attributes(file)
    for k in keys(attrs)
        println("$(k) => $(string(attrs["$k"][]))")
    end
end



function extract_metadata(file::HDF5.File, filename::String = "metadata.txt")
    attrs = HDF5.attributes(file)
    open(filename, "w") do io
        for k in keys(attrs)
            println(io, "$(k) => $(string(attrs["$k"][]))\n")
        end
    end
end

function source_code(file::HDF5.File)
    attribs = HDF5.attributes(file)
    code=attribs["[ENV] Code"][];
    open("source_code.jl", "w") do io
        write(io, code)
    end
end