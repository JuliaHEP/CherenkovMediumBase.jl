using CherenkovMediumBase
using CherenkovMediumBase: es_scattering, es_scattering_integral, es_scattering_cumulative
using Test
using Random
using StaticArrays

Random.seed!(1234)


function test_interface(medium::MediumProperties)
    wavelength = 400.

    @test pressure(medium) isa Number
    @test temperature(medium) isa Number
    @test material_density(medium) isa Number
    @test scattering_length(medium, wavelength) isa Number
    @test absorption_length(medium, wavelength) isa Number
    @test group_refractive_index(medium, wavelength) isa Number
    @test phase_refractive_index(medium, wavelength) isa Number
    @test dispersion(medium, wavelength) isa Number
    @test cherenkov_angle(medium, wavelength) isa Number
    @test group_velocity(medium, wavelength) isa Number
    @test sample_scattering_function(medium) isa Number

end

# Define a dispersion model



# Define a concrete MediumProperties subtype for testing
struct MockMediumProperties <: MediumProperties
    salinity::Float64
    temperature::Float64
    pressure::Float64
    dispersion_model::QuanFryDispersion
    scattering_model::KopelevichScatteringModel
    absoption_model::InterpolatedAbsorptionModel
    function MockMediumProperties(salinity, temperature, pressure)
        return new(
            salinity,
            temperature,
            pressure,
            QuanFryDispersion(salinity, temperature, pressure),
            KopelevichScatteringModel(HenyeyGreenStein(0.8), 0.1, 0.2),
            InterpolatedAbsorptionModel(SA[1., 2], SA[3., 4])
            )
    end
end

struct FaultyInterfaceMockMedium <: MediumProperties end

CherenkovMediumBase.get_absorption_model(medium::MockMediumProperties) = medium.absoption_model
CherenkovMediumBase.get_scattering_model(medium::MockMediumProperties) = medium.scattering_model
CherenkovMediumBase.get_dispersion_model(medium::MockMediumProperties) = medium.dispersion_model

CherenkovMediumBase.pressure(medium::MockMediumProperties) = medium.pressure
CherenkovMediumBase.temperature(medium::MockMediumProperties) = medium.temperature
CherenkovMediumBase.material_density(medium::MockMediumProperties) = 3.


@testset "CherenkovMediumBase.jl" begin
    medium = MockMediumProperties(35.0, 15.0, 1.)
    test_interface(medium)

    @testset "Faulty Interface" begin
        medium = FaultyInterfaceMockMedium()
        @test_throws ErrorException get_absorption_model(medium)
        @test_throws ErrorException get_scattering_model(medium)
        @test_throws ErrorException get_dispersion_model(medium)
        @test_throws ErrorException pressure(medium)
        @test_throws ErrorException temperature(medium)
        @test_throws ErrorException material_density(medium)
        @test_throws ErrorException group_velocity(medium, 100)
        @test_throws ErrorException sample_scattering_function(medium)
        @test_throws ErrorException cherenkov_angle(medium, 100)
        @test_throws ErrorException dispersion(medium, 100)
        @test_throws ErrorException group_refractive_index(medium, 100)
        @test_throws ErrorException phase_refractive_index(medium, 100)
        @test_throws ErrorException absorption_length(medium, 100)
        @test_throws ErrorException scattering_length(medium, 100)
    end




end

@testset "Dispersion Models" begin
    @testset "QuanFry Dispersion" begin
        disp_model = QuanFryDispersion(35.0, 15.0, 1.0)

        @testset "Phase Refractive Index" begin
            wavelength = 500.0 # nm
            n_phase = phase_refractive_index(disp_model, wavelength)
            @test isapprox(n_phase, 1.34, atol=0.01)
        end

        @testset "Dispersion" begin
            wavelength = 500.0 # nm
            disp = dispersion(disp_model, wavelength)
            @test isapprox(disp, -0.0001, atol=0.0001)
        end
        
        # Test group_refractive_index function
        @testset "Group Refractive Index" begin
            medium = MockMediumProperties(35.0, 15.0, 1.0)
            wavelength = 500.0 # nm
            n_group = group_refractive_index(medium, wavelength)
            @test isapprox(n_group, 1.368, atol=0.01)
        end

        # Test group_velocity function
        @testset "Group Velocity" begin
            medium = MockMediumProperties(35.0, 15.0, 1.0)
            wavelength = 500.0 # nm
            v_group = group_velocity(medium, wavelength)
            @test isapprox(v_group, 0.22, atol=0.01)
        end

        # Test cherenkov_angle function
        @testset "Cherenkov Angle" begin
            medium = MockMediumProperties(35.0, 15.0, 1.0)
            wavelength = 500.0 # nm
            angle = cherenkov_angle(medium, wavelength)
            @test isapprox(angle, 0.73, atol=0.01)
        end
    end
end

@testset "Scattering" begin

    @testset "Scattering Functions" begin
        @testset "HenyeyGreenStein" begin
            scattering_function = HenyeyGreenStein(0.95)
            @test rand(scattering_function) ≈ 0.99014 atol=1e-5
        end

        @testset "EinsteinSmoluchowsky" begin
            scattering_function = EinsteinSmoluchowsky(0.835)
            @test rand(scattering_function) ≈ -0.41788 atol=1e-5
            @test @inferred rand(scattering_function) isa Float64
        end    

        @testset "Type Stability" begin
            
            let b = 0.835, cos_theta = 0.5

                @test @inferred es_scattering(cos_theta, b) isa Float64
                @test @inferred es_scattering_integral(cos_theta, b) isa Float64
                @test @inferred es_scattering_cumulative(cos_theta, b) isa Float64
            end

            let b = 0.835f0,  cos_theta = 0.5f0

                @test @inferred es_scattering(cos_theta, b) isa Float32
                @test @inferred es_scattering_integral(cos_theta, b) isa Float32
                @test @inferred es_scattering_cumulative(cos_theta, b) isa Float32
            end
        end
    end

    @testset "KopelevichScatteringModel" begin
        scattering_model = KopelevichScatteringModel(HenyeyGreenStein(0.95), 7.5E-3, 7.5E-3)

        @testset "isbits" begin
            @test isbits(scattering_model)
        end

        # Test scattering_length function
        @testset "scattering_length" begin
            @test scattering_length(scattering_model, 400.0) ≈ 37.69 atol=1e-2
        end
    end

    @testset "MixedHGES" begin
        scattering_function = MixedHGES(0.95, 0.835, 0.5)

        @testset "isbits" begin
            @test isbits(scattering_function)
        end
    end

    @testset "WavelengthIndependentScatteringModel" begin
        sca_len = 25.0
        scattering_function = MixedHGES(0.95, 0.835, 0.5)
        scattering_model = WavelengthIndependentScatteringModel(scattering_function, sca_len)

        @test scattering_length(scattering_model, 450.0) == sca_len 
        @test scattering_length(scattering_model, 550.0) == sca_len 
        @test sample_scattering_function(scattering_model) isa Number
    end

    
end

@testset "Absorption Models" begin
    @testset "Interpolated Absorption Model" begin
        wavelengths = SA[400.0, 500.0, 600.0]
        absorption = SA[10.0, 20.0, 30.0]
        abs_model = InterpolatedAbsorptionModel(wavelengths, absorption)

        @test absorption_length(abs_model, 450.0) ≈ 15.0 atol=1e-2
        @test absorption_length(abs_model, 550.0) ≈ 25.0 atol=1e-2
    end

    @testset "Wavelength Independent Absorption Model" begin
        abs_length = 25.0
        abs_model = WavelengthIndependentAbsorptionModel(abs_length)

        @test absorption_length(abs_model, 450.0) ≈ abs_length atol=1e-2
        @test absorption_length(abs_model, 550.0) ≈ abs_length atol=1e-2
    end
end
