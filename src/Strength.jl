module Strength

using Parameters, AISIS100, CUFSM, Unitful, Serialization

using ..Properties

@with_kw struct DiagonalTension

    Ag::Float64
    Anet::Float64
    Tny::Float64
    eTny::Float64
    Tnu::Float64
    eTnu::Float64
    eTn::Float64

end

@with_kw struct DiagonalCompression

    Pcre::Float64
    Pcrℓ::Float64
    Py::Float64
    Pne::Float64
    ePne::Float64
    Pnℓ::Float64
    ePnℓ::Float64

end



@with_kw struct ChordCompression

    Ag::Float64 
    Py::Float64 
    Pne::Float64 
    ePne::Float64 
    Pcrℓ::Float64 
    local_buckling_properties::CUFSM.Model
    Pnℓ::Float64  
    ePnℓ::Float64 
    Pcrd::Float64 
    distortional_buckling_properties::CUFSM.Model
    Pnd::Float64 
    ePnd::Float64 
    ePn::Float64 

end

@with_kw struct ChordTension

    Ag::Float64
    Tn_y::Float64
    eTn_y::Float64
    Anet::Float64 
    Tn_u::Float64 
    eTn_u::Float64
    eTn::Float64 

end




@with_kw struct Components

    diagonal_tensile_strength::Vector{DiagonalTension}
    diagonal_compressive_strength::Vector{DiagonalCompression}
    chord_compressive_strength::ChordCompression
    chord_tensile_strength::ChordTension
    top_chord_connection_strength::Vector{Union{ReinforcedConnection, UnreinforcedConnection}}
    bottom_chord_connection_strength::Vector{Union{ReinforcedConnection, UnreinforcedConnection}}

    eRn_diagonals::Vector{Float64}
    eRn_top_chord_connections::Vector{Float64}
    eRn_bottom_chord_connections::Vector{Float64}
    eRn_top_chord::Float64
    eRn_bottom_chord::Float64

end


function calculate_diagonal_tensile_strength(bolt_hole_diameter, t, diagonal_section_properties, joist_material_properties, design_code)

    Ag = diagonal_section_properties.A

    Tny, eTny = AISIS100.v16.d21(Ag=Ag, Fy=joist_material_properties.fy, design_code=design_code)

    Anet = Ag - 2 * bolt_hole_diameter * t

    Tnu, eTnu = AISIS100.v16.d31(An=Anet, Fu=joist_material_properties.fu, design_code=design_code)

    eTn = minimum([eTny, eTnu])

    diagonal_tensile_strength = DiagonalTension(Ag=Ag, Anet=Anet, Tny=Tny, eTny=eTny, Tnu=Tnu, eTnu=eTnu, eTn=eTn)

    return diagonal_tensile_strength

end

function calculate_diagonal_global_buckling(model, diagonal_section_assignments, diagonal_bracing, diagonal_dimensions, diagonal_section_geometry, joist_material_properties)

    diagonal_index = findall(section->occursin("diagonal", section), model.inputs.element.cross_section)

    Pcre_diagonal = Vector{Float64}(undef, length(diagonal_index))
    diagonal_global_buckling = Vector{CUFSM.Model}(undef, length(diagonal_index))

    for i in eachindex(diagonal_index)

        section_index = findfirst(section->section==diagonal_section_assignments[i], diagonal_dimensions.name)

        if diagonal_bracing[i] == "unbraced"
            L = model.properties.L[diagonal_index[i]]
        elseif diagonal_bracing[i] == "braced"
            L = model.properties.L[diagonal_index[i]]/2
        end

        Pcre_diagonal[i], diagonal_global_buckling[i] = Properties.calculate_diagonal_global_buckling_load(diagonal_section_geometry[section_index], diagonal_dimensions.t[section_index], joist_material_properties.E, joist_material_properties.ν, L)

    end

    return Pcre_diagonal, diagonal_global_buckling

end

function calculate_all_diagonal_tensile_strengths(diagonal_section_assignments, diagonal_dimensions, diagonal_section_properties, joist_material_properties, bolt_hole_diameter, design_code)

    diagonal_tensile_strength = Vector{Strength.DiagonalTension}(undef, length(diagonal_section_assignments))

    for i in eachindex(diagonal_section_assignments)

        section_index = findfirst(section->section==diagonal_section_assignments[i], diagonal_dimensions.name)

        diagonal_tensile_strength[i] = Strength.calculate_diagonal_tensile_strength(bolt_hole_diameter, diagonal_dimensions.t[section_index], diagonal_section_properties[section_index], joist_material_properties, design_code)

    end

    return diagonal_tensile_strength

end


function calculate_all_diagonal_compressive_strengths(diagonal_section_assignments, diagonal_dimensions, Pcre, Pcrℓ, diagonal_section_properties, fy, design_code)

    diagonal_compressive_strength = Vector{DiagonalCompression}(undef, length(diagonal_section_assignments))
    for i in eachindex(diagonal_section_assignments)

        section_index = findfirst(section->section==diagonal_section_assignments[i], diagonal_dimensions.name)

        diagonal_compressive_strength[i] = calculate_diagonal_compressive_strength(Pcre[i], Pcrℓ[section_index], diagonal_section_properties[section_index].A, fy, design_code)

    end

    return diagonal_compressive_strength 

end

function calculate_diagonal_compressive_strength(Pcre, Pcrℓ, A, fy, design_code)

    Py = fy * A

	Pne, ePne = AISIS100.v16.e2(Fcre=Pcre/A, Fy=fy, Ag=A, design_code=design_code)
    Pnℓ, ePnℓ = AISIS100.v16.e321(Pne=Pne, Pcrℓ=Pcrℓ, design_code=design_code)

    diagonal_compressive_strength = Strength.DiagonalCompression(Pcre=Pcre, Pcrℓ=Pcrℓ, Py=Py, Pne=Pne, ePne=ePne, Pnℓ=Pnℓ, ePnℓ=ePnℓ)
   
    return diagonal_compressive_strength

		
end

	



function calculate_chord_compressive_strength(chord_section_properties, joist_material_properties, design_code, chord_dimensions)

    Ag = chord_section_properties.A 
    Py_chord = Ag * joist_material_properties.fy

    #Assume the deck full braces the top chord in compression which means that $P_{ne} = P_y$ in AISI S100-16 Section E2
    Pne_chord, ePne_chord = AISIS100.v16.e2(Fcre=joist_material_properties.fy*100000, Fy=joist_material_properties.fy, Ag=Ag, design_code=design_code)


    lengths = collect(0.25*3.0:3.0/20:1.25*3.0)
    chord_local_buckling = deserialize("/Users/crismoen/.julia/dev/MarkoJIT/assets/chord_local_buckling_properties")
    Pcrℓ_chord = deserialize("/Users/crismoen/.julia/dev/MarkoJIT/assets/chord_Pcrl")
    # Pcrℓ_chord, chord_local_buckling =  Properties.calculate_chord_cross_section_buckling_load(chord_dimensions, joist_material_properties, lengths)

    Pnℓ_chord, ePnℓ_chord = AISIS100.v16.e321(Pne=Pne_chord, Pcrℓ=Pcrℓ_chord, design_code=design_code)


    lengths = collect(26.0:1.0:34.0)
    # Pcrd_chord, chord_distortional_buckling =  Properties.calculate_chord_cross_section_buckling_load(chord_dimensions, joist_material_properties, lengths)
        
    chord_distortional_buckling = deserialize("/Users/crismoen/.julia/dev/MarkoJIT/assets/chord_distortional_buckling_properties")
    Pcrd_chord = deserialize("/Users/crismoen/.julia/dev/MarkoJIT/assets/chord_Pcrd")

    Pnd_chord, ePnd_chord = AISIS100.v16.e41(Py=Py_chord, Pcrd=Pcrd_chord, design_code=design_code)

    ePn = minimum([ePnℓ_chord, ePnd_chord])

    chord_compressive_strength = ChordCompression(Ag=Ag, Py=Py_chord, Pne=Pne_chord, ePne=ePne_chord, Pcrℓ=Pcrℓ_chord, local_buckling_properties=chord_local_buckling, Pnℓ=Pnℓ_chord, ePnℓ=ePnℓ_chord, Pcrd=Pcrd_chord, distortional_buckling_properties=chord_distortional_buckling, Pnd=Pnd_chord, ePnd=ePnd_chord, ePn=ePn)

    return chord_compressive_strength

end

function calculate_chord_tensile_strength(chord_section_properties, chord_dimensions, joist_material_properties, design_code)

    Ag = chord_section_properties.A
    Tn_y, eTn_y = AISIS100.v16.d21(Ag=Ag, Fy=joist_material_properties.fy, design_code=design_code)

    Anet = Ag - chord_dimensions.flange_punchout_width * chord_dimensions.t

    Tn_u, eTn_u = AISIS100.v16.d31(An=Anet, Fu=joist_material_properties.fu, design_code=design_code)

    eTn = minimum([eTn_y, eTn_u])

    chord_tensile_strength = ChordTension(Ag=Ag, Tn_y=Tn_y, eTn_y=eTn_y, Anet=Anet, Tn_u=Tn_u, eTn_u=eTn_u, eTn=eTn)

    return chord_tensile_strength

end

end  #module 