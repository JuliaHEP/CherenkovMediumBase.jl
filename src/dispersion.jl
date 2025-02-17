export DispersionModel
export QuanFryDispersion

abstract type AbstractDispersionModel end

"""
    phase_refractive_index(disperion_model::AbstractDispersionModel)

Calculate the phase refractive index for a given dispersion model.
"""
phase_refractive_index(disperion_model::AbstractDispersionModel) = _not_implemented(disperion_model)

"""
    dispersion(disperion_model::AbstractDispersionModel)

Calculate the dispersion for a given dispersion model.
"""
dispersion(disperion_model::AbstractDispersionModel) = _not_implemented(disperion_model)


"""
    calc_quan_fry_params(salinity::Real, temperature::Real, pressure::Real)

Helper function to get the parameters for the Quan & Fry formula as function of
salinity, temperature and pressure.
"""
function calc_quan_fry_params(
    salinity::T,
    temperature::T,
    pressure::T) where {T <: Real}

    n0 = 1.31405
    n1 = 1.45e-5
    n2 = 1.779e-4
    n3 = 1.05e-6
    n4 = 1.6e-8
    n5 = 2.02e-6
    n6 = 15.868
    n7 = 0.01155
    n8 = 0.00423
    n9 = 4382
    n10 = 1.1455e6

    a01 = (
        n0
        +
        (n2 - n3 * temperature + n4 * temperature^2) * salinity
        -
        n5 * temperature^2
        +
        n1 * pressure
    )
    a2 = n6 + n7 * salinity - n8 * temperature
    a3 = -n9
    a4 = n10

    return T(a01), T(a2), T(a3), T(a4)
end


"""
    _refractive_index_fry(wavelength, quan_fry_params)

The phase refractive index of sea water according to a model
from Quan & Fry.

wavelength is given in nm, salinity in permille, temperature in °C and pressure in atm

The original model is taken from:
X. Quan, E.S. Fry, Appl. Opt., 34, 18 (1995) 3477-3480.

An additional term describing pressure dependence was included according to:
Wolfgang H.W.A. Schuster, "Measurement of the Optical Properties of the Deep
Mediterranean - the ANTARES Detector Medium.",
PhD thesis (2002), St. Catherine's College, Oxford
downloaded Jan 2011 from: http://www.physics.ox.ac.uk/Users/schuster/thesis0098mmjhuyynh/thesis.ps

Adapted from clsim (https://github.com/claudiok/clsim)
"""
function _refractive_index_fry(
    wavelength::Real,
    quan_fry_params::Tuple{U,U,U,U}
) where {U<:Real}
    a01, a2, a3, a4 = quan_fry_params
    x = one(wavelength) / wavelength
    return oftype(wavelength, a01 + x*a2 + x^2*a3 + x^3 * a4)
end

"""
    dispersion_fry(wavelength, a2, a3, a4)

Calculate the dispersion for the Quan & Fry dispersion model.
"""
function dispersion_fry(wavelength::T, a2, a3, a4) where {T<:Number}
    x = one(T) / wavelength

    return T(a2 + T(2) * x * a3 + T(3) * x^2 * a4) * T(-1) / wavelength^2
end

"""
    QuanFryDispersion{T <: Real}

Struct to hold parameters for the Quan & Fry dispersion model.
"""
struct QuanFryDispersion{T <: Real} <: AbstractDispersionModel
    a01::T
    a2::T
    a3::T
    a4::T
end

"""
    QuanFryDispersion(a01, a2, a3, a4)

Constructor for QuanFryDispersion with given parameters.
"""
function QuanFryDispersion(a01, a2, a3, a4)
    return QuanFryDispersion(promote(a01, a2, a3, a4)...)
end

"""
    QuanFryDispersion(salinity, temperature, pressure)

Constructor for QuanFryDispersion using salinity, temperature, and pressure.
"""
function QuanFryDispersion(salinity, temperature, pressure)
    return QuanFryDispersion(calc_quan_fry_params(salinity, temperature, pressure)...)
end

"""
    phase_refractive_index(disperion_model::QuanFryDispersion, wavelength)

Calculate the phase refractive index for the Quan & Fry dispersion model.
"""
function phase_refractive_index(disperion_model::QuanFryDispersion, wavelength)
    return _refractive_index_fry(wavelength, (disperion_model.a01, disperion_model.a2, disperion_model.a3, disperion_model.a4))
end

"""
    dispersion(dispersion_model::QuanFryDispersion, wavelength)

Calculate the dispersion for the Quan & Fry dispersion model.
"""
function dispersion(dispersion_model::QuanFryDispersion, wavelength)
    return dispersion_fry(wavelength, dispersion_model.a2, dispersion_model.a3, dispersion_model.a4)
end


"""
    group_velocity(dispersion_model, wavelength)
Return the group_velocity in m/ns at `wavelength`.

`wavelength` is expected to be in units nm.
"""
function group_velocity(dispersion_model::AbstractDispersionModel, wavelength)
    global c_vac_m_ns
    T = typeof(wavelength)

    # Explicitely convert everything to the type of wavelength
    # This is useful to avoid double precision if this function is called in a CUDA kernel

    ref_ix::T = phase_refractive_index(dispersion_model, wavelength)
    λ_0::T = ref_ix * wavelength
    T(c_vac_m_ns) / (ref_ix - λ_0 * dispersion(dispersion_model, wavelength))
end

"""
    cherenkov_angle(medium, wavelength)
Calculate the cherenkov angle (in rad) for `wavelength`.

`wavelength` is expected to be in units nm.
"""
function cherenkov_angle(dispersion_model::AbstractDispersionModel, wavelength)
    return acos(one(typeof(wavelength)) / phase_refractive_index(dispersion_model, wavelength))
end

"""
    group_refractive_index(medium, wavelength)
Return the group refractive index at `wavelength`.

`wavelength` is expected to be in units nm.
"""
function group_refractive_index(dispersion_model::AbstractDispersionModel, wavelength)
    T = typeof(wavelength)

    ref_ix = phase_refractive_index(dispersion_model, wavelength)
    return ref_ix / (T(1.0) + dispersion(dispersion_model, wavelength) * wavelength / ref_ix)
end