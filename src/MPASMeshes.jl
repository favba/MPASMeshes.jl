module MPASMeshes

using Zeros, TensorsLite, ImmutableVectors, VoronoiMeshes, VoronoiOperators
using Reexport

@reexport using VoronoiMeshes
export MPASMesh

import Random

include("tangential_velocity_reconstruction.jl")

struct MPASMesh{OnSphere, max_nEdges, max_nEdges2, TI, TF, TZ} <: AbstractVoronoiMesh{OnSphere, max_nEdges, TI, TF, TZ}
    cells::Cells{OnSphere, max_nEdges, TI, TF, TZ}
    vertices::Vertices{OnSphere, max_nEdges, TI, TF, TZ}
    edges::Edges{OnSphere, max_nEdges, TI, TF, TZ}
    tanVelRec::TangentialVelocityReconstructionGeneric{max_nEdges2, TI, TF}
    attributes::Dict{String, Union{String, TF, Int32, Int64}}

    function MPASMesh(
            cells::Cells{OnSphere, max_nEdges, TI, TF, TZ},
            vertices::Vertices{OnSphere, max_nEdges, TI, TF, TZ},
            edges::Edges{OnSphere, max_nEdges, TI, TF, TZ},
            tanVelRec::TangentialVelocityReconstruction{max_nEdges2, TI, TF},
            attributes::Dict{String, Union{String, TF, Int32, Int64}}
        ) where {OnSphere, max_nEdges, TI, TF, TZ, max_nEdges2}

        if !(get_diagram(cells) === get_diagram(vertices) === get_diagram(edges))
            throw(DimensionMismatch("`Cell`, `Vertices` and `Edges` structs are not based on the same Voronoi diagram"))
        end

        tvr = TangentialVelocityReconstructionGeneric(tanVelRec)
        attributes["tangential_velocity_reconstruction_method"] = tvr.type

        return new{OnSphere, max_nEdges, max_nEdges2, TI, TF, TZ}(cells, vertices, edges, tvr, attributes)
    end
end

for N in 6:9
    for N2 in 10:(2N)
        precompile(MPASMesh, (Cells{false, N, Int32, Float64, Zero}, Vertices{false, N, Int32, Float64, Zero}, Edges{false, N, Int32, Float64, Zero}, TangentialVelocityReconstructionThuburn{N2, Int32, Float64}, Dict{String, Union{String, Float64, Int32, Int64}}))
        precompile(MPASMesh, (Cells{true, N, Int32, Float64, Float64}, Vertices{true, N, Int32, Float64, Float64}, Edges{true, N, Int32, Float64, Float64}, TangentialVelocityReconstructionThuburn{N2, Int32, Float64}, Dict{String, Union{String, Float64, Int32, Int64}}))
    end
end

function generate_attributes(mesh::AbstractVoronoiMesh{S, N, TI, TF}) where {S, N, TI, TF}
    attrib = Dict{String, Union{String, TF, Int32, Int64}}()

    attrib["mesh_spec"] = "1.0"
    attrib["file_id"] = Random.randstring(8)

    if S
        attrib["on_a_sphere"] = "YES"
        attrib["sphere_radius"] = mesh.sphere_radius
        attrib["is_periodic"] = "NO"
    else
        attrib["on_a_sphere"] = "NO"
        attrib["sphere_radius"] = zero(TF)
        attrib["is_periodic"] = "YES"
        attrib["x_period"] = mesh.x_period
        attrib["y_period"] = mesh.y_period
    end

    return attrib
end

MPASMesh(mesh::AbstractVoronoiMesh, tanVelRec::TangentialVelocityReconstruction) = MPASMesh(mesh.cells, mesh.vertices, mesh.edges, tanVelRec, generate_attributes(mesh))
MPASMesh(mesh::AbstractVoronoiMesh) = MPASMesh(mesh, TangentialVelocityReconstructionThuburn(mesh))

MPASMesh(N::Integer, lx::Real, ly::Real; density::F = x -> 1, max_iter::Integer = 20000, rtol::Real = 1.0e-10) where {F} = MPASMesh(VoronoiMesh(N, lx, ly; density = density, max_iter = max_iter, rtol = rtol))
MPASMesh(points::AbstractVector{<:Vec}, lx::Real, ly::Real; density::F = x -> 1, max_iter::Integer = 20000, rtol::Real = 1.0e-10) where {F} = MPASMesh(VoronoiMesh(points, lx, ly; density = density, max_iter = max_iter, rtol = rtol))

function compute_mpas_fields!(mesh::MPASMesh)
    mesh.cells.area
    mesh.cells.longitude
    mesh.cells.latitude
    mesh.vertices.area
    mesh.vertices.kiteAreas
    mesh.vertices.longitude
    mesh.vertices.latitude
    mesh.edges.length
    mesh.edges.lengthDual
    mesh.edges.angle
    mesh.edges.longitude
    mesh.edges.latitude
    return mesh
end

function VoronoiMeshes.save(filename, obj::MPASMesh; kwds...)
    if VoronoiMeshes.is_netcdf_ext(filename)
        compute_mpas_fields!(obj)
        VoronoiMeshes.save_to_netcdf(filename, obj; kwds...)
    else
        error("Unsupported file extension: $filename")
    end
    return nothing
end

for N in 6:9
    for N2 in 8:(2N)
        precompile(compute_mpas_fields!, (MPASMesh{false, N, N2, Int32, Float64, Zeros.Zero},))
        precompile(compute_mpas_fields!, (MPASMesh{true, N, N2, Int32, Float64, Float64},))
    end
end

function write_coeffs_reconstruct_to_grid(velRecon::CellVelocityReconstruction, filename::AbstractString)
    if VoronoiMeshes.is_netcdf_ext(filename)
        write_coeffs_reconstruct_to_grid_netcdf(filename, velRecon)
    else
        error("Unsupported file extension: $filename")
    end
    return nothing
end

#implemented in NCDatasetsExt.jl
function write_coeffs_reconstruct_to_grid_netcdf end


end
