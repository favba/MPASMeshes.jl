module ComoniconExt

using NCDatasets, Comonicon, VoronoiMeshes
import MPASMeshes


"""
    regenerate_mesh(grid::String, output::String; [tangent_reconstruction::String = "thuburn", area_type::String = "mimetic") -> nothing

# Introduction
 
Regenerate the mesh given by `grid` and write it to `output`.
The new mesh will have the same Voronoi Diagram and cells / vertices ordering of the original mesh.
Computable information, such as areas, and all edge information will be completely recreated, and any indexing problem with the original mesh will be fixed.
Optionally, the `tangent_reconstruction` string specifies which method to use to compute the tangential velocity reconstruction weights, and the `area_type` string specifies the area computation method.
The "thuburn" `tangent_reconstruction` scheme is the standard scheme used in MPAS (see Thuburn (2009) and Ringler et. al. (2010)).
The "peixoto" and "peixoto_old" schemes are based on the work of Peixoto (2016) and differ only for spherical meshes.
The former first performs a full projection of the whole cell onto the plane tangential to the cell's generator point and then computes the reconstruction weights of the projected planar cell.
The latter first computes the reconstruction weights using 3D vectors that connect the relevant points and projects the end result onto the tangential direction of the edge.
The "lsq2" scheme projects the cell on the tangent plane and computes the weights by using a quadratic (penalized when needed) least squares regression.

# Arguments

- `grid`: NetCDF file containing Voronoi grid information.
- `output`: NetCDF file name for output mesh.

# Options

- `-t, --tangent_reconstruction`: Edge tangential velocity reconstruction method. Available options: "peixoto", "peixoto_old", "lsq2", and "thuburn". Default is "thuburn".
- `-a, --area_type`: Areas computation method. Available options: "geometric" for the true geometrical area; "mimetic" for area values suitable for the TRiSK mimetic scheme. Default is "mimetic".
"""
Comonicon.@main function regenerate_mesh(grid::String, output::String;
    tangent_reconstruction::String="thuburn",
    area_type::String="mimetic")

    MPASMeshes.regenerate_mesh(grid, output; reconstruction_method=tangent_reconstruction, area_type=area_type)

    return 0
end

precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:tangent_reconstruction, :area_type),Tuple{String,String}},typeof(regenerate_mesh),String,String})
precompile(regenerate_mesh, (String, String, @NamedTuple{tangent_reconstruction::String, area_type::String}))

precompile(Tuple{typeof(command_main),Vector{String}})
precompile(Tuple{typeof(command_main)})

MPASMeshes.regenerate_mesh(args::Vector{String}) = command_main(args)

precompile(Tuple{typeof(MPASMeshes.regenerate_mesh),Vector{String}})

end
