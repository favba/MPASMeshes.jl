for N in 10:16
    precompile(Tuple{typeof(NCDatasetsExt.read_tanVelRecon), Base.Val{N}, Array{UInt8, 1}, NCDatasets.NCDataset{Nothing, Base.Missing}})
end
precompile(Tuple{typeof(Base.setindex!), Base.Dict{String, Union{Float64, Int32, Int64, String}}, String, String})
precompile(Tuple{typeof(Base.setindex!), Base.Dict{String, Union{Float64, Int32, Int64, String}}, Float64, String})
for N in 10:16
    precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.CFVariable{Int32, 2, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Base.OneTo{Int64}})
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_sizecheck!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
    precompile(Tuple{typeof(NCDatasets._write_data_to_nc), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
    precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.CFVariable{Float64, 2, NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.ReinterpretArray{Float64, 2, NTuple{N, Float64}, Array{NTuple{N, Float64}, 1}, true}, Base.OneTo{Int64}, Base.OneTo{Int64}})
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_sizecheck!), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Float64, 2, NTuple{N, Float64}, Array{NTuple{N, Float64}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock!), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Float64, 2, NTuple{N, Float64}, Array{NTuple{N, Float64}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
    precompile(Tuple{typeof(NCDatasets._write_data_to_nc), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Float64, 2, NTuple{N, Float64}, Array{NTuple{N, Float64}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
end
for N in 6:9
    for N2 in 8:(2N)
        precompile(VoronoiMeshes.save, (String, MPASMesh{false, N, N2, Int32, Float64, Zeros.Zero}))
        precompile(VoronoiMeshes.save, (String, MPASMesh{true, N, N2, Int32, Float64, Float64}))
    end
end

