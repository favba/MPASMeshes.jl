module NCDatasetsExt

using MPASMeshes, VoronoiOperators, NCDatasets, TensorsLite, Zeros, ImmutableVectors
import MPASMeshes: VoronoiMeshes
import VoronoiMeshes: copy_matrix_to_tuple_vector!, save_to_netcdf!

function read_tanVelRecon(::Val{maxEdgesOnEdges}, nEdges::AbstractVector, ncfile) where {maxEdgesOnEdges}
    n = length(nEdges)
    inds_tuple = Vector{NTuple{maxEdgesOnEdges, Int32}}(undef, n)
    inds = ImmutableVectorArray(inds_tuple, nEdges)
    edgesOnEdgeArray = ncfile["edgesOnEdge"][:, :]::Matrix{Int32}
    copy_matrix_to_tuple_vector!(inds_tuple, edgesOnEdgeArray)

    w_tuple = Vector{NTuple{maxEdgesOnEdges, Float64}}(undef, n)
    w = ImmutableVectorArray(w_tuple, nEdges)
    weightsArray = ncfile["weightsOnEdge"][:, :]::Matrix{Float64}
    copy_matrix_to_tuple_vector!(w_tuple, weightsArray)

    method = if haskey(ncfile.attrib,"tangential_velocity_reconstruction_method")
        ncfile.attrib["tangential_velocity_reconstruction_method"]::String
    else
        "Thuburn"
    end

    return MPASMeshes.TangentialVelocityReconstructionGeneric(inds, w, inds.length, method)
end

function read_tanVelRecon(ncfile)
    nEdges = UInt8.(ncfile["nEdgesOnEdge"][:]::Vector{Int32})
    max_n_edges = Int(maximum(nEdges))
    return read_tanVelRecon(Val(max_n_edges), nEdges, ncfile)
end

function _MPASMesh(ncfile::NCDatasets.NCDataset)
    voro_mesh = VoronoiMesh(ncfile, false)
    tanVelRecon = read_tanVelRecon(ncfile)
    attrib = convert(Dict{String, Union{String, Float64, Int32, Int64}}, ncfile.attrib)
    return MPASMesh(voro_mesh.cells, voro_mesh.vertices, voro_mesh.edges, tanVelRecon, attrib)
end

function MPASMeshes.MPASMesh(ncfile::NCDatasets.NCDataset, warn_issues::Bool = true)
    mesh = _MPASMesh(ncfile)
    if warn_issues
        Threads.@spawn VoronoiMeshes.warn_mesh_issues($mesh)
    end
    return mesh
end

function MPASMeshes.MPASMesh(file_name::String, warn_issues::Bool = true)
    f = NCDataset(file_name)
    try
        MPASMeshes.MPASMesh(f, warn_issues)
    finally
        close(f)
    end
end

function write_tanVelRecon_data!(ds::NCDataset, tanVelRecon::MPASMeshes.TangentialVelocityReconstructionGeneric{N, TI, TF}) where {N, TI, TF}
    ds.dim["maxEdges2"] = N

    defVar(
        ds, "nEdgesOnEdge", TI.(tanVelRecon.nEdges),
        ("nEdges",), attrib = [
            "units" => "-",
            "long_name" => "Number of edges involved in reconstruction of tangential velocity for an edge.",
        ]
    )

    defVar(
        ds, "edgesOnEdge", reinterpret(reshape, TI, tanVelRecon.indices.data),
        ("maxEdges2", "nEdges"), attrib = [
            "units" => "-",
            "long_name" => "IDs of edges involved in reconstruction of tangential velocity for an edge.",
        ]
    )

    defVar(
        ds, "weightsOnEdge", reinterpret(reshape, TF, tanVelRecon.weights.data),
        ("maxEdges2", "nEdges"), attrib = [
            "units" => "-",
            "long_name" => "Weights used in reconstruction of tangential velocity for an edge.",
        ]
    )

    return ds
end

function save_to_netcdf!(ds::NCDataset, mesh::MPASMesh{S, N, N2, TI}; force3D::Bool = true, write_computed::Bool = true) where {S, N, N2, TI}

    save_to_netcdf!(ds, VoronoiMesh(mesh.cells, mesh.vertices, mesh.edges), force3D = true, write_computed = true)

    defVar(
        ds, "indexToCellID", TI.(1:mesh.cells.n),
        ("nCells",), attrib = [
            "units" => "-",
            "long_name" => "Mapping from local array index to global cell ID",
        ]
    )

    defVar(
        ds, "indexToVertexID", TI.(1:mesh.vertices.n),
        ("nVertices",), attrib = [
            "units" => "-",
            "long_name" => "Mapping from local array index to global vertex ID",
        ]
    )

    defVar(
        ds, "indexToEdgeID", TI.(1:mesh.edges.n),
        ("nEdges",), attrib = [
            "units" => "-",
            "long_name" => "Mapping from local array index to global edge ID",
        ]
    )

    write_tanVelRecon_data!(ds, mesh.tanVelRec)

    for (key, val) in pairs(mesh.attributes)
        ds.attrib[key] = val
    end

    return ds
end

include("precompile_NCDatasetsExt.jl")

end # module
