struct TangentialVelocityReconstructionGeneric{N_MAX, TI, TF} <: VoronoiOperators.TangentialVelocityReconstruction{N_MAX, TI, TF}
    indices::ImVecArray{N_MAX, TI, 1}
    weights::ImVecArray{N_MAX, TF, 1}
    nEdges::Vector{UInt8}
    type::String
end

TangentialVelocityReconstructionGeneric(tvr::TangentialVelocityReconstructionGeneric) = tvr

function TangentialVelocityReconstructionGeneric(tvr::TangentialVelocityReconstructionThuburn)
    return TangentialVelocityReconstructionGeneric(tvr.indices, tvr.weights, tvr.weights.length, "Thuburn")
end

function TangentialVelocityReconstructionGeneric(tvr::TangentialVelocityReconstructionPeixoto)
    return TangentialVelocityReconstructionGeneric(tvr.indices, tvr.weights, tvr.weights.length, "Peixoto")
end
