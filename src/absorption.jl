using Interpolations
using StaticArrays
export AbstractAbsorptionModel

abstract type AbstractAbsorptionModel end

absorption_length(model::AbstractAbsorptionModel, wavelength::Real) = _not_implemented(model)


struct InterpolatedAbsorptionModel{T, ITP} <: AbstractAbsorptionModel
    interpolation::ITP
end

function InterpolatedAbsorptionModel(wavelengths::SVector{N, T}, absorption::SVector{N, T}) where {N, T}
    itp = linear_interpolation(wavelengths, absorption, extrapolation_bc = T(0))
    return InterpolatedAbsorptionModel{T, typeof(itp)}(itp)
end

function absorption_length(model::InterpolatedAbsorptionModel, wavelength::Real)
    model.interpolation(wavelength)
end


"""
    WavelengthIndependentAbsorptionModel{T} <: AbstractAbsorptionModel

Struct for wavelength-independent absorption model.
"""
struct WavelengthIndependentAbsorptionModel{T} <: AbstractAbsorptionModel
    absorption_length::T
end


function absorption_length(model::WavelengthIndependentAbsorptionModel, wavelength::Real)
    return model.absorption_length
end