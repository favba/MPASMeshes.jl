struct TangentialVelocityReconstructionGeneric{N_MAX, TI, TF} <: VoronoiOperators.TangentialVelocityReconstruction{N_MAX, TI, TF}
    indices::SmVecArray{N_MAX, TI, 1}
    weights::SmVecArray{N_MAX, TF, 1}
    nEdges::Vector{Int16}
    method_name::String
end

TangentialVelocityReconstructionGeneric(tvr::TangentialVelocityReconstructionGeneric) = tvr

function TangentialVelocityReconstructionGeneric(tvr::TangentialVelocityReconstruction)
    return TangentialVelocityReconstructionGeneric(tvr.indices, tvr.weights, tvr.weights.length, method_name(tvr))
end
