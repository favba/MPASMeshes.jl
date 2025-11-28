module ComoniconExt

using NCDatasets, Comonicon, VoronoiMeshes
import MPASMeshes


"""
    regenerate_mesh(grid::String, output::String; [tangent_reconstruction::String = "trisk", area_type::String = "mimetic") -> nothing

# Intro 
 
Regenerate the mesh given by `grid` and write it to `output`.
The new mesh will have the same Voronoi Diagram and cells / vertices ordering of the original mesh.
The edge information will be completely recreated, and any indexing problem will be fixed.
Opitionally, the `tangent_reconstruction` string specifies which method to use to compute the tangential velocity reconstruction weights.
and the `area_type` string specifies the area computation methods.

# Arguments

- `grid`: NetCDF file containing Voronoi grid information.
- `output`: NetCDF file name for output mesh.

# Options

- `-t, --tangent_reconstruction`: Edge tangential velocity reconstruction method. Available options: "peixoto", "peixoto_old", "lsq2", and "thuburn". Default is "thuburn".
- `-a, --area_type`: Areas computation method. Available options: "geometric" for the true geometrical area; "mimetic" for area values suitable for the trisk mimetic scheme. Default is "mimetic".
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
