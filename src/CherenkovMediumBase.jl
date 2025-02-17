module CherenkovMediumBase

export MediumProperties

export get_dispersion_model, get_scattering_model, get_absorption_model
export pressure, temperature, material_density, radiation_length
export sample_scattering_function, absorption_length, scattering_length
export group_refractive_index, phase_refractive_index
export dispersion, cherenkov_angle, group_velocity


const c_vac_m_ns = 0.299792458

abstract type MediumProperties end

_not_implemented(type) = error("Not implemented for type $(typeof(type))")

# Interface functions

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
    get_absorption_model(medium::MediumProperties)
Return the absorption model for a given medium.
"""
get_absorption_model(medium::MediumProperties) = _not_implemented(medium)

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

# Utility functions

"""
    sample_scattering_function(medium::MediumProperties)

Return a cos(scattering angle) sampled from the scattering function of the medium.
"""
function sample_scattering_function(medium::MediumProperties)
    model = get_scattering_model(medium)
    return sample_scattering_function(model)
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
function absorption_length(medium::MediumProperties, wavelength)
    model = get_absorption_model(medium)
    return absorption_length(model, wavelength)
end

"""
    phase_refractive_index(medium::MediumProperties, wavelength)

Return the phase refractive index at `wavelength`.
`wavelength` is expected to be in units nm.
"""
function phase_refractive_index(medium::MediumProperties, wavelength)
    model = get_dispersion_model(medium)
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
    disp = get_dispersion_model(medium)
    return cherenkov_angle(disp, wavelength)
end

"""
    group_velocity(medium, wavelength)
Return the group_velocity in m/ns at `wavelength`.

`wavelength` is expected to be in units nm.
"""
function group_velocity(medium::MediumProperties, wavelength)
    disp = get_dispersion_model(medium)
    return group_velocity(disp, wavelength)
end


"""
    group_refractive_index(medium, wavelength)
Return the group refractive index at `wavelength`.

`wavelength` is expected to be in units nm.
"""
function group_refractive_index(medium, wavelength)
    disp = get_dispersion_model(medium)
    return group_refractive_index(disp, wavelength)
end

include("utils.jl")
include("dispersion.jl")
include("scattering.jl")
include("absorption.jl")

end # module
