# CherenkovMediumBase

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliahep.github.io/CherenkovMediumBase.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliahep.github.io/CherenkovMediumBase.jl/dev/)
[![Build Status](https://github.com/JuliaHEP/CherenkovMediumBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaHEP/CherenkovMediumBase.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/juliahep/CherenkovMediumBase.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliahep/CherenkovMediumBase.jl)

CherenkovMediumBase provides abstract types and functions for calculating various optical properties of different media. It provides functionalities to compute properties such as refractive index, scattering length, absorption length, and more, based on the medium's characteristics.

## Main Features

- Calculation of phase and group refractive indices
- Scattering and absorption length computations
- Cherenkov angle and group velocity calculations
- Support for different scattering functions including Henyey-Greenstein and Simplified-Liu models

## Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add("CherenkovMediumBase")
```

## Implementing a Concrete MediumProperties Subtype

To implement a concrete subtype of `MediumProperties`, you need to define a new type and implement the required methods. Here is an example:

```julia

using CherenkovMediumBase
struct Water <: MediumProperties
    salinity::Float64
    temperature::Float64
    pressure::Float64
    vol_conc_small_part::Float64
    vol_conc_large_part::Float64
    mean_scattering_angle::Float64
    dispersion_model::QuanFryDispersion
    scattering_model::KopelevichScatteringModel
    function Water(salinity, temperature, pressure, vol_conc_small_part, vol_conc_large_part, mean_scattering_angle)
        return new(
            salinity,
            temperature,
            pressure,
            vol_conc_small_part,
            vol_conc_large_part,
            mean_scattering_angle,            
            QuanFryDispersion(salinity, temperature, pressure),
            KopelevichScatteringModel(HenyeyGreenStein(mean_scattering_angle), vol_conc_small_part, vol_conc_large_part)
            )
        end
end

CherenkovMediumBase.pressure(medium::Water) = medium.pressure
CherenkovMediumBase.temperature(medium::Water) = medium.temperature
CherenkovMediumBase.material_density(medium::Water) = 1.
CherenkovMediumBase.absorption_length(medium::Water, wavelength) = 20.
CherenkovMediumBase.sample_scattering_function(medium::Water) =  rand(medium.scattering_model.scattering_function)


function CherenkovMediumBase.phase_refractive_index(medium::Water, wavelength)
    return phase_refractive_index(medium.dispersion_model, wavelength)
end

function CherenkovMediumBase.dispersion(medium::Water, wavelength)
    return dispersion(medium.dispersion_model, wavelength)
end

function CherenkovMediumBase.scattering_length(medium::Water, wavelength)
    return scattering_length(medium.scattering_model, wavelength)
end

# Implement other required methods similarly
```

This example defines a `Water` type which uses the Quan-Fry dispersion model and the Kopalevich scattering model

## Usage

Here is a simple example of how to use the package:

```julia
using CherenkovMediumBase

# Define a medium (example)
medium = Water(34.82, 4, 100, 0.005, 0.005, 0.95)

# Calculate phase refractive index at a given wavelength
wavelength = 500.0 # nm
n_phase = phase_refractive_index(medium, wavelength)

println("Phase refractive index at $wavelength nm: $n_phase")
```

For more detailed usage and examples, please refer to the [documentation](https://juliahep.github.io/CherenkovMediumBase.jl/stable/).



