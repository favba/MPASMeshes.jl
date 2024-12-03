module MPASMeshes

using Zeros, TensorsLite, VoronoiMeshes, VoronoiOperators
using Reexport

@reexport using VoronoiMeshes
export MPASMesh

import Random

const TVRecon{TI, TF} = VoronoiOperators.TangentialVelocityReconstruction{max_Edges2, TI, TF} where {max_Edges2}

struct MPASMesh{OnSphere, max_nEdges, TI, TF, TZ, TVR <: TVRecon{TI, TF}} <: AbstractVoronoiMesh{OnSphere, max_nEdges, TI, TF, TZ}
    cells::Cells{OnSphere, max_nEdges, TI, TF, TZ}
    vertices::Vertices{OnSphere, max_nEdges, TI, TF, TZ}
    edges::Edges{OnSphere, max_nEdges, TI, TF, TZ}
    tanVelRec::TVR
    attributes::Dict{String, Union{String, TF, Int32, Int64}}

    function MPASMesh(
            cells::Cells{OnSphere, max_nEdges, TI, TF, TZ},
            vertices::Vertices{OnSphere, max_nEdges, TI, TF, TZ},
            edges::Edges{OnSphere, max_nEdges, TI, TF, TZ},
            tanVelRec::TVR,
            attributes::Dict{String, Union{String, TF, Int32, Int64}}
        ) where {OnSphere, max_nEdges, TI, TF, TZ, TVR <: TVRecon{TI, TF}}

        if !(get_diagram(cells) === get_diagram(vertices) === get_diagram(edges))
            throw(DimensionMismatch("`Cell`, `Vertices` and `Edges` structs are not based on the same Voronoi diagram"))
        end

        return new{OnSphere, max_nEdges, TI, TF, TZ, TVR}(cells, vertices, edges, tanVelRec, attributes)
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

MPASMesh(mesh::AbstractVoronoiMesh, tanVelRec::TVRecon) = MPASMesh(mesh.cells, mesh.vertices, mesh.edges, tanVelRec, generate_attributes(mesh))
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
    for N2 in 10:(2N)
        precompile(compute_mpas_fields!, (MPASMesh{false, N, Int32, Float64, Zeros.Zero, TangentialVelocityReconstructionThuburn{N2, Int32, Float64}},))
        precompile(compute_mpas_fields!, (MPASMesh{true, N, Int32, Float64, Float64, TangentialVelocityReconstructionThuburn{N2, Int32, Float64}},))
        precompile(Tuple{typeof(VoronoiMeshes.save), String, MPASMesh{false, N, Int32, Float64, Zeros.Zero, VoronoiOperators.TangentialVelocityReconstructionThuburn{N2, Int32, Float64}}})
        precompile(Tuple{typeof(VoronoiMeshes.save), String, MPASMesh{true, N, Int32, Float64, Float64}, VoronoiOperators.TangentialVelocityReconstructionThuburn{N2, Int32, Float64}})
    end
end

end
