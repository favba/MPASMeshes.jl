module MPASMeshes

using Reexport, Zeros, TensorsLite, SmallCollections
@reexport using VoronoiMeshes
using  VoronoiOperators

export MPASMesh
export regenerate_mesh

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
        attributes["tangential_velocity_reconstruction_method"] = tvr.method_name

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

MPASMesh(N::Integer, lx::Real, ly::Real; kwd...) = MPASMesh(VoronoiMesh(N, lx, ly; kwd...))
MPASMesh(points::AbstractVector{<:AbstractVec}, lx::Real, ly::Real; kwd...) = MPASMesh(VoronoiMesh(points, lx, ly; kwd...))

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

function VoronoiMeshes.save(filename::String, obj::MPASMesh; kwds...)
    name, ext = Base.Filesystem.splitext(filename)
    if ext == ".nc"
        compute_mpas_fields!(obj)
        VoronoiMeshes.save_to_netcdf(filename, obj; kwds...)
        write(name*".graph.info", String(take!(graph_partition(obj))))
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

function write_coeffs_reconstruct_to_grid(velRecon::Union{<:CellVelocityReconstruction,<:VertexVelocityReconstruction}, filename::String)
    _, ext = Base.Filesystem.splitext(filename)
    if ext == ".nc"
        write_coeffs_reconstruct_to_grid_netcdf(filename, velRecon)
    else
        error("Unsupported file extension: $filename")
    end
    return nothing
end

function write_coeffs_scalar_reconstruct_to_grid(edgeToCell::EdgeToCellTransformation, filename::String)
    _, ext = Base.Filesystem.splitext(filename)
    if ext == ".nc"
        write_coeffs_scalar_reconstruct_to_grid_netcdf(filename, edgeToCell)
    else
        error("Unsupported file extension: $filename")
    end
    return nothing
end

#implemented in NCDatasetsExt.jl
function write_coeffs_reconstruct_to_grid_netcdf end
function write_coeffs_scalar_reconstruct_to_grid_netcdf end

"""
    regenerate_mesh(input_mesh_name::String, out_file_name::String, [method::String = "trisk"]) -> nothing

Regenerate the mesh given by `input_mesh_name` and write it to `out_file_name`.
The new mesh will have the same Voronoi Diagram and cells / vertices ordering of the original mesh.
The edge information will be completely recreated, and any indexing problem will be fixed.

Opitionally, the `method` string specifies which method to use to compute the tangential velocity reconstruction weights, cell velocity
reconstruction weights (if needed), and vertex velocity reconstruction weights (if needed).
Currently, valid options are "trisk", "peixoto", "peixoto_perot_perot", "peixoto_perot_lsq2", "peixoto_lsq2_lsq2", "thuburn_perot",
"thuburn_perot_perot", "thuburn_perot_lsq2", "thuburn_lsq2_lsq2", and "lsq2_lsq2_lsq2".
"""
function regenerate_mesh(inputfile::String, outputname::String, method="trisk")
    v_mesh = VoronoiMesh(fix_diagram!(VoronoiDiagram(inputfile)))
    if method == "trisk" #equivalent to "thuburn_ringler_gassman"
        mpas_mesh = MPASMesh(v_mesh)
        save(outputname, mpas_mesh)
    elseif method == "peixoto" # equivalent to "peixoto_perot"
        mpas_mesh = MPASMesh(v_mesh, TangentialVelocityReconstructionPeixoto(v_mesh))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
    elseif method == "peixoto_old"
        mpas_mesh = MPASMesh(v_mesh, VoronoiOperators.TangentialVelocityReconstructionPeixotoOld(v_mesh))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(VoronoiOperators.CellVelocityReconstructionPerotOld(mpas_mesh), outputname)
    elseif method == "peixoto_perot_perot"
        mpas_mesh = MPASMesh(v_mesh, TangentialVelocityReconstructionPeixoto(v_mesh))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
        write_coeffs_reconstruct_to_grid(VertexVelocityReconstructionPerot(mpas_mesh), outputname)
    elseif method == "peixoto_perot_lsq2"
        mpas_mesh = MPASMesh(v_mesh, TangentialVelocityReconstructionPeixoto(v_mesh))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
        write_coeffs_reconstruct_to_grid(VertexVelocityReconstructionLSq2(mpas_mesh), outputname)
    elseif method == "peixoto_lsq2_lsq2"
        mpas_mesh = MPASMesh(v_mesh, TangentialVelocityReconstructionPeixoto(v_mesh))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionLSq2(mpas_mesh), outputname)
        write_coeffs_reconstruct_to_grid(VertexVelocityReconstructionLSq2(mpas_mesh), outputname)
    elseif method == "thuburn_perot"
        mpas_mesh = MPASMesh(v_mesh)
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
    elseif method == "thuburn_perot_perot"
        mpas_mesh = MPASMesh(v_mesh)
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
        write_coeffs_reconstruct_to_grid(VertexVelocityReconstructionPerot(mpas_mesh), outputname)
    elseif method == "thuburn_perot_lsq2"
        mpas_mesh = MPASMesh(v_mesh)
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
        write_coeffs_reconstruct_to_grid(VertexVelocityReconstructionLSq2(mpas_mesh), outputname)
    elseif method == "thuburn_lsq2_lsq2"
        mpas_mesh = MPASMesh(v_mesh)
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionLSq2(mpas_mesh), outputname)
        write_coeffs_reconstruct_to_grid(VertexVelocityReconstructionLSq2(mpas_mesh), outputname)
    elseif method == "lsq2_lsq2_lsq2"
        cR = CellVelocityReconstructionLSq2(v_mesh)
        mpas_mesh = MPASMesh(v_mesh, TangentialVelocityReconstructionVelRecon(v_mesh, cR))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(cR, outputname)
        write_coeffs_reconstruct_to_grid(VertexVelocityReconstructionLSq2(mpas_mesh), outputname)
    elseif method == "trisk_tangent_ke_area"
        mpas_mesh = MPASMesh(v_mesh)
        save(outputname, mpas_mesh)
        write_coeffs_scalar_reconstruct_to_grid(EdgeToCellArea(mpas_mesh), outputname)
    elseif method == "trisk_tangent_ke_lsq2"
        mpas_mesh = MPASMesh(v_mesh)
        save(outputname, mpas_mesh)
        write_coeffs_scalar_reconstruct_to_grid(EdgeToCellLSq2(mpas_mesh), outputname)
    elseif method == "trisk_tangent_ke_lsq3"
        mpas_mesh = MPASMesh(v_mesh)
        save(outputname, mpas_mesh)
        write_coeffs_scalar_reconstruct_to_grid(EdgeToCellLSq3(mpas_mesh), outputname)
    elseif method == "peixoto_tangent_ke_area"
        mpas_mesh = MPASMesh(v_mesh, TangentialVelocityReconstructionPeixoto(v_mesh))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
        write_coeffs_scalar_reconstruct_to_grid(EdgeToCellArea(mpas_mesh), outputname)
    elseif method == "peixoto_tangent_ke_lsq2"
        mpas_mesh = MPASMesh(v_mesh, TangentialVelocityReconstructionPeixoto(v_mesh))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
        write_coeffs_scalar_reconstruct_to_grid(EdgeToCellLSq2(mpas_mesh), outputname)
    elseif method == "peixoto_tangent_ke_lsq3"
        mpas_mesh = MPASMesh(v_mesh, TangentialVelocityReconstructionPeixoto(v_mesh))
        save(outputname, mpas_mesh)
        write_coeffs_reconstruct_to_grid(CellVelocityReconstructionPerot(mpas_mesh), outputname)
        write_coeffs_scalar_reconstruct_to_grid(EdgeToCellLSq3(mpas_mesh), outputname)
    else
        error("Method '$method' not implemented")
    end
    return nothing
end

end
