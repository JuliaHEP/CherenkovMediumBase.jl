module CherenkovMediumBase

export MediumProperties
export pressure, temperature, material_density, mean_scattering_angle
export sample_scattering_function, absorption_length, group_refractive_index, phase_refractive_index
export dispersion, cherenkov_angle, group_velocity, radiation_length
export scattering_function
export scattering_length


const c_vac_m_ns = 0.299792458

abstract type MediumProperties end

_not_implemented(type) = error("Not implemented for type $(typeof(type))")

"""
    get_dispersion_model(medium::MediumProperties)
Return the dispersion model for a given medium.
"""
get_dispersion_model(medium::MediumProperties) = _not_implemented(medium)


"""
    get_scattering_model(medium::MediumProperties)
Return the scattering model for a given medium.
"""
get_scattering_model(medium::MediumProperties) = _not_implemented(medium)

"""
    pressure(medium::MediumProperties)

This function returns the pressure for a given medium.
"""
pressure(medium::MediumProperties) = _not_implemented(medium)

"""
    temperature(medium::MediumProperties)

This function returns the temperature for a given medium.
"""
temperature(medium::MediumProperties) = _not_implemented(medium)

"""
    material_density(medium::MediumProperties)

This function returns the material density for a given medium.
"""
material_density(medium::MediumProperties) = _not_implemented(medium)

"""
    radiation_length(medium::MediumProperties)

This function returns the radiation length for a given medium.
"""
radiation_length(medium::MediumProperties) = _not_implemented(medium)

"""
    sample_scattering_function(medium::MediumProperties)

Return a scattering angle sampled from the scattering function of the medium.
"""
function sample_scattering_function(medium::MediumProperties)
    model = scattering_model(medium)


end

"""
    scattering_length(medium::MediumProperties, wavelength)

Return scattering length at `wavelength` in units m.
`wavelength` is expected to be in units nm. Returned length is in units m.
"""
function scattering_length(medium::MediumProperties, wavelength)
    model = get_scattering_model(medium)
    return scattering_length(model, wavelength)
end

"""
    absorption_length(medium::MediumProperties, wavelength)

Return absorption length at `wavelength` in units m.
`wavelength` is expected to be in units nm.
"""
absorption_length(medium::MediumProperties, wavelength) = _not_implemented(medium)

"""
    phase_refractive_index(medium::MediumProperties, wavelength)

Return the phase refractive index at `wavelength`.
`wavelength` is expected to be in units nm.
"""
function phase_refractive_index(medium::MediumProperties, wavelength)
    model =get_dispersion_model(medium)
    return phase_refractive_index(model, wavelength)
end

"""
    dispersion(medium::MediumProperties, wavelength)

Return the dispersion dn/dλ at `wavelength` in units 1/nm.
`wavelength` is expected to be in units nm.
"""
function dispersion(medium::MediumProperties, wavelength)
    model = get_dispersion_model(medium)
    return dispersion(model, wavelength)
end

"""
    cherenkov_angle(medium, wavelength)
Calculate the cherenkov angle (in rad) for `wavelength`.

`wavelength` is expected to be in units nm.
"""
function cherenkov_angle(medium::MediumProperties, wavelength)
    return acos(one(typeof(wavelength)) / phase_refractive_index(medium, wavelength))
end

"""
    group_velocity(medium, wavelength)
Return the group_velocity in m/ns at `wavelength`.

`wavelength` is expected to be in units nm.
"""
function group_velocity(medium::MediumProperties, wavelength)
    global c_vac_m_ns
    T = typeof(wavelength)

    # Explicitely convert everything to the type of wavelength
    # This is useful to avoid double precision if this function is called in a CUDA kernel

    ref_ix::T = phase_refractive_index(medium, wavelength)
    λ_0::T = ref_ix * wavelength
    T(c_vac_m_ns) / (ref_ix - λ_0 * dispersion(medium, wavelength))
end


"""
    group_refractive_index(medium, wavelength)
Return the group refractive index at `wavelength`.

`wavelength` is expected to be in units nm.
"""
function group_refractive_index(medium, wavelength)
    T = typeof(wavelength)

    ref_ix = phase_refractive_index(medium, wavelength)
    return ref_ix / (T(1.0) + dispersion(medium, wavelength) * wavelength / ref_ix)
end

include("utils.jl")
include("dispersion.jl")
include("scattering.jl")

end # module
