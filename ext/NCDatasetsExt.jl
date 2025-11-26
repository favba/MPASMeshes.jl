module NCDatasetsExt

using TensorsLite, Zeros, SmallCollections, VoronoiMeshes, VoronoiOperators, MPASMeshes
import VoronoiMeshes: copy_matrix_to_fixedvector_vector!, save_to_netcdf!, SmallVectorArray
using NCDatasets
using PrecompileTools

function read_tanVelRecon(::Val{maxEdgesOnEdges}, nEdges::AbstractVector, ncfile) where {maxEdgesOnEdges}
    n = length(nEdges)
    inds_tuple = Vector{FixedVector{maxEdgesOnEdges, Int32}}(undef, n)
    inds = SmallVectorArray(inds_tuple, nEdges)
    edgesOnEdgeArray = ncfile["edgesOnEdge"][:, :]::Matrix{Int32}
    copy_matrix_to_fixedvector_vector!(inds_tuple, edgesOnEdgeArray)

    w_tuple = Vector{FixedVector{maxEdgesOnEdges, Float64}}(undef, n)
    w = SmallVectorArray(w_tuple, nEdges)
    weightsArray = ncfile["weightsOnEdge"][:, :]::Matrix{Float64}
    copy_matrix_to_fixedvector_vector!(w_tuple, weightsArray)

    method = if haskey(ncfile.attrib,"tangential_velocity_reconstruction_method")
        ncfile.attrib["tangential_velocity_reconstruction_method"]::String
    else
        "Thuburn"
    end

    return MPASMeshes.TangentialVelocityReconstructionGeneric(inds, w, inds.length, method)
end

function read_tanVelRecon(ncfile)
    nEdges = Int16.(ncfile["nEdgesOnEdge"][:]::Vector{Int32})
    max_n_edges = Int(maximum(nEdges))
    return read_tanVelRecon(Val(max_n_edges), nEdges, ncfile)
end

function _MPASMesh(ncfile::NCDatasets.NCDataset)
    voro_mesh = VoronoiMesh(ncfile, warn_issues=false)
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

function MPASMeshes.write_coeffs_reconstruct_to_grid_netcdf(filename::String, velRecon::CellVelocityReconstruction{N_MAX, TI, TF, TZ}) where {N_MAX, TI, TF, TZ}

    NCDataset(filename, "a") do ds

        ncells = length(velRecon.weights)

        weights = if (TZ === TF)
            reshape(reinterpret(TF, velRecon.weights.data), (3, N_MAX, ncells))
        else
            wdata = velRecon.weights.data
            w = Array{TF}(undef, (3, N_MAX, ncells))
            @inbounds for k in Base.OneTo(ncells)
                t = wdata[k]
                for j in Base.OneTo(N_MAX)
                    v = t[j]
                    w[1, j, k] = v.x
                    w[2, j, k] = v.y
                    w[3, j, k] = zero(TF)
                end
            end
            w
        end

        defVar(
            ds, "coeffs_reconstruct", weights,
            ("R3", "maxEdges", "nCells"), attrib = [
                "units" => "-",
                "long_name" => "Coefficients to reconstruct velocity vectors at cell centers",
            ]
        )

        ds.attrib["cell_velocity_reconstruction_method"] = method_name(velRecon)
    end

    return
end

function MPASMeshes.write_coeffs_reconstruct_to_grid_netcdf(filename::String, velRecon::VertexVelocityReconstruction{TI, TF, TZ}) where {TI, TF, TZ}

    NCDataset(filename, "a") do ds

        nvertices = length(velRecon.weights)

        weights = if (TZ === TF)
            reshape(reinterpret(TF, velRecon.weights), (3, 3, nvertices))
        else
            wdata = velRecon.weights
            w = Array{TF}(undef, (3, 3, nvertices))
            @inbounds for k in Base.OneTo(nvertices)
                t = wdata[k]
                for j in Base.OneTo(3)
                    v = t[j]
                    w[1, j, k] = v.x
                    w[2, j, k] = v.y
                    w[3, j, k] = zero(TF)
                end
            end
            w
        end

        defVar(
            ds, "vertex_coeffs_reconstruct", weights,
            ("R3", "vertexDegree", "nVertices"), attrib = [
                "units" => "-",
                "long_name" => "Coefficients to reconstruct velocity vectors at vertices",
            ]
        )

        ds.attrib["vertex_velocity_reconstruction_method"] = method_name(velRecon)
    end

    return
end


function MPASMeshes.write_coeffs_scalar_reconstruct_to_grid_netcdf(filename::String, edgeToCell::EdgeToCellTransformation{N_MAX, TI, TF}) where {N_MAX, TI, TF}

    NCDataset(filename, "a") do ds

        ncells = length(edgeToCell.weights)

        weights = reshape(reinterpret(TF, edgeToCell.weights.data), (N_MAX, ncells))

        defVar(
            ds, "coeffs_scalar_reconstruct", weights,
            ("maxEdges", "nCells"), attrib = [
                "units" => "-",
                "long_name" => "Coefficients to reconstruct edge scalar field at cell centers",
            ]
        )

        ds.attrib["edge_to_cell_method"] = method_name(edgeToCell)
    end

    return
end

@setup_workload begin
    files = ("../test/spherical_unif_grid_4000km.nc", "../test/spherical_unif_grid_2000km.nc", "../test/spherical_grid_500km.nc")
    for f in files
        bdir = string(Base.@__DIR__, "/")
        fin = string(bdir, f)
        fout1 = string(bdir, "asda1_.nc")
        fout1_graph = string(bdir, "asda1_.graph.info")
        fout2 = string(bdir, "asda2_.nc")
        fout2_graph = string(bdir, "asda2_.graph.info")
        fout2_1 = string(bdir, "asda2_1.nc")
        fout2_1_graph = string(bdir, "asda2_1.graph.info")
        fout3 = string(bdir, "asda3_.nc")
        fout3_graph = string(bdir, "asda3_.graph.info")
        fout4 = string(bdir, "asda4_.nc")
        fout4_graph = string(bdir, "asda4_.graph.info")

        @compile_workload begin
            regenerate_mesh(fin, fout1)
            regenerate_mesh(fin, fout2, reconstruction_method = "peixoto")
            regenerate_mesh(fin, fout2_1, reconstruction_method = "peixoto_old")
            regenerate_mesh(fin, fout3, area_type = "geometric")
            regenerate_mesh(fin, fout4, reconstruction_method = "lsq2")
        end

        Base.Filesystem.rm(fout1)
        Base.Filesystem.rm(fout1_graph)
        Base.Filesystem.rm(fout2)
        Base.Filesystem.rm(fout2_graph)
        Base.Filesystem.rm(fout2_1)
        Base.Filesystem.rm(fout2_1_graph)
        Base.Filesystem.rm(fout3)
        Base.Filesystem.rm(fout3_graph)
        Base.Filesystem.rm(fout4)
        Base.Filesystem.rm(fout4_graph)
    end

end

end # module
