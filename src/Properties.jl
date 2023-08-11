module Properties

using CUFSM, CrossSection, Parameters, Serialization


@with_kw struct JoistMaterial

  fy::Float64
  fu::Float64
  E::Float64
  ν::Float64

end

@with_kw struct Section

    chord::CUFSM.SectionPropertiesObject
    diagonals::Vector{CUFSM.SectionPropertiesObject}

end


function calculate_chord_section_properties(chord)

    x = [chord.centerline_cross_section_coordinates[i][1] for i in eachindex(chord.centerline_cross_section_coordinates)]
    y = [chord.centerline_cross_section_coordinates[i][2] for i in eachindex(chord.centerline_cross_section_coordinates)]

    y = -y  #flip so that the cross section looks like top chord
    y = y .- chord.t/2 

    section_properties = CUFSM.cutwp_prop2([x y], [1:length(x)-1 2:length(x) ones(Float64, length(x)-1) * chord.t])

	return section_properties 

end






function calculate_diagonal_local_buckling_load(diagonal_section_geometry, diagonal_dimensions, joist_material_properties)

    #local buckling 
    P = 1.0
    Mxx = 0.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    lengths = range(2.0, 4.0, 5)  #hard coded
    neigs = 1

    diagonal_section_local_buckling = Vector{CUFSM.Model}(undef, 0)
    Pcrℓ = Vector{Float64}(undef, 0)

    for i in eachindex(diagonal_section_geometry)

        x_center = [diagonal_section_geometry[i].center[j][1] for j in eachindex(diagonal_section_geometry[i].center)]
        y_center = [diagonal_section_geometry[i].center[j][2] for j in eachindex(diagonal_section_geometry[i].center)]
        diagonal_section_local_buckling = push!(diagonal_section_local_buckling, CUFSM.Tools.open_section_analysis(x_center, y_center, diagonal_dimensions.t[i], lengths, joist_material_properties.E, joist_material_properties.ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs))

        Pcrℓ = push!(Pcrℓ, minimum([diagonal_section_local_buckling[i].curve[j,1][2] for j=1:length(lengths)]))

    end

    # serialize("/Users/crismoen/.julia/dev/MarkoJIT/assets/diagonal_Pcrl", Pcrℓ)
    # serialize("/Users/crismoen/.julia/dev/MarkoJIT/assets/diagonal_section_local_buckling", diagonal_section_local_buckling)

    return Pcrℓ, diagonal_section_local_buckling

end



function calculate_diagonal_global_buckling_load(diagonal_section_geometry, t, E, ν, L)

    #local buckling 
    P = 1.0
    Mxx = 0.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    lengths = [L]
    neigs = 1

    x_center = [diagonal_section_geometry.center[j][1] for j in eachindex(diagonal_section_geometry.center)]
    y_center = [diagonal_section_geometry.center[j][2] for j in eachindex(diagonal_section_geometry.center)]
    diagonal_global_buckling = CUFSM.Tools.open_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    Pcre = diagonal_global_buckling.curve[1,1][2]


    return Pcre, diagonal_global_buckling

end


function calculate_chord_cross_section_buckling_load(chord_dimensions, joist_material_properties, lengths)

    #local buckling 
    P = 1.0
    Mxx = 0.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    # lengths = collect(0.25*3.0:3.0/20:1.25*3.0)

    x = [chord_dimensions. centerline_cross_section_coordinates[j][1] for j in eachindex(chord_dimensions. centerline_cross_section_coordinates)]
    y = [chord_dimensions. centerline_cross_section_coordinates[j][2] for j in eachindex(chord_dimensions. centerline_cross_section_coordinates)]

    y = -y  #flip so that the cross section looks like top chord
    y = y .- chord_dimensions.t/2 

    chord_buckling = CUFSM.Tools.open_section_analysis(x, y, chord_dimensions.t, lengths, joist_material_properties.E, joist_material_properties.ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    Pcr_chord = minimum([chord_buckling.curve[i][2] for i in eachindex(chord_buckling.curve)])

    return Pcr_chord, chord_buckling

end


end #module 