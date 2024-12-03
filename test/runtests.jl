using NCDatasets, DelaunayTriangulation, MPASMeshes
using Test

@testset "MPASMeshes.jl" begin
    for mp in (MPASMesh(33, 1.0, 1.0, rtol = 1.0e-7), MPASMesh("spherical_grid_500km.nc"))
        save("lakdf.nc", mp)
        mp_read = MPASMesh("lakdf.nc")
        Base.Filesystem.rm("lakdf.nc")

        c1 = mp.cells
        c2 = mp_read.cells
        @test c1.position == c2.position
        @test c1.area == c2.area
        @test c1.longitude == c2.longitude
        @test c1.latitude == c2.latitude

        v1 = mp.vertices
        v2 = mp_read.vertices
        @test v1.position == v2.position
        @test v1.area == v2.area
        @test v1.kiteAreas == v2.kiteAreas
        @test v1.longitude == v2.longitude
        @test v1.latitude == v2.latitude

        e1 = mp.edges
        e2 = mp_read.edges
        @test e1.position == e2.position
        @test e1.length == e2.length
        @test e1.lengthDual == e2.lengthDual
        @test e1.angle == e2.angle
        @test e1.longitude == e2.longitude
        @test e1.latitude == e2.latitude
    end
end
