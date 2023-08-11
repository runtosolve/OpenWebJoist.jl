module OpenWebJoist

using Parameters, InstantFrame, CrossSection, Serialization 


export Geometry
include("Geometry.jl")
using .Geometry

export Properties
include("Properties.jl")
using .Properties

export Model
include("Model.jl")
using .Model

export Strength
include("Strength.jl")
using .Strength

@with_kw struct Inputs

    design_code::String
    joist_dimensions::OpenWebJoist.Geometry.JoistDimensions
    chord_dimensions::OpenWebJoist.Geometry.ChordDimensions
    diagonal_dimensions::OpenWebJoist.Geometry.DiagonalDimensions
    diagonal_sections::Vector{String}
    bearing_seat_dimensions::OpenWebJoist.Geometry.BearingSeatDimensions
    chord_splice_dimensions::NamedTuple{(:t, :B, :H, :R, :location), Tuple{Float64, Float64, Float64, Float64, NamedTuple{(:top, :bottom), Tuple{Vector{Float64}, Vector{Float64}}}}}
    girder_dimensions::OpenWebJoist.Geometry.Girder
    girder_material_properties::NamedTuple{(:fy, :fu), Tuple{Float64, Float64}}
    bearing_seat_weld_properties::NamedTuple{(:length, :Fxx, :t, :num_lines), Tuple{Float64, Float64, Float64, Int64}}
    joist_material_properties::OpenWebJoist.Properties.JoistMaterial
    top_chord_connections::Vector{String}
    bottom_chord_connections::Vector{String}

end

@with_kw struct ComponentDemands

    diagonals::Vector{Float64}
    top_chord::Vector{Float64}
    bottom_chord::Vector{Float64}
    bearing_seat_weld::Vector{Float64}
    bearing_seat_compression_field::Vector{Float64}
    bearing_seat_chord_connection::Vector{Float64}

end


@with_kw struct DemandToCapacity

    diagonals::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Vector{Float64}, Int64, Float64}}
    top_chord_connections::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Vector{Float64}, Int64, Float64}}
    bottom_chord_connections::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Vector{Float64}, Int64, Float64}}
    top_chord::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Vector{Float64}, Int64, Float64}}
    bottom_chord::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Vector{Float64}, Int64, Float64}}
    top_chord_splice::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Vector{Float64}, Int64, Float64}} 
    bottom_chord_splice::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Vector{Float64}, Int64, Float64}}
    bearing_seat_weld_1::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Float64, Int64, Float64}}
    bearing_seat_weld_2::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Float64, Int64, Float64}}
    bearing_seat_compression_field_1::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Float64, Int64, Float64}}
    bearing_seat_compression_field_2::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Float64, Int64, Float64}}
    bearing_seat_chord_connection_1::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Float64, Int64, Float64}}
    bearing_seat_chord_connection_2::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Float64, Int64, Float64}}
    joist_deflection::NamedTuple{(:DC, :max_location, :max_DC), Tuple{Float64, Int64, Float64}}

    all::Vector{Float64}
    labels::Vector{String}
    max::Float64
    controlling::String
    failure_location::Int64

end

@with_kw struct JoistSpan

    inputs::Inputs
    properties::Properties.Section
    model::InstantFrame.Model
    strength::Strength.Components
    deflection::Float64
    demand::ComponentDemands
    demand_to_capacity::DemandToCapacity
    limit_state_loads::Vector{Tuple{String, Float64}}
    design_load::Float64

end

function calculate_chord_splice_DC(element_type, model, chord_splice_location, eRn_chord_splice, chord_demand, joist)

    chord_elements = findall(element->element==element_type, model.inputs.element.cross_section)

    chord_element_index = Vector{Int}(undef, 0)

    joist_ends = Geometry.define_joist_ends(joist.span_length, joist.node_spacing)
    for i in eachindex(chord_splice_location)

        if element_type == "top chord"
            chord_element_index_all = findall(location->location<chord_splice_location[i], cumsum(model.properties.L[chord_elements]))
        elseif element_type == "bottom chord"
            chord_element_index_all = findall(location->location<chord_splice_location[i], cumsum(model.properties.L[chord_elements]) .+ (joist_ends[1] - joist.node_spacing/2 - joist.bottom_chord_tail_length))
            # print(chord_element_index_all[end]+1)
        end

        if !isempty(chord_element_index_all)
            chord_element_index = push!(chord_element_index, chord_element_index_all[end]+1)  #need +1 here to land in the correct chord element 
        end
    end



    if isempty(chord_element_index)
        chord_splice_DC = calculate_demand_to_capacity([0.0], eRn_chord_splice) #no splice 
    else
        chord_splice_DC = calculate_demand_to_capacity(abs.(chord_demand[chord_element_index])./1000, eRn_chord_splice)
    end

    return chord_splice_DC

end

function calculate_demand_to_capacity(demand, capacity)

    DC = demand ./ capacity
    max_DC_location = argmax(DC)
    max_DC = maximum(DC)

    DC_details = (DC=DC, max_location=max_DC_location, max_DC = max_DC)

    return DC_details

end


function evaluate_joist_span(design_code, joist_dimensions, chord_dimensions, diagonal_dimensions, bolt_properties, shield_plate_dimensions, diagonal_sections, diagonal_bracing, bearing_seat_dimensions, chord_splice_dimensions,  girder_dimensions, girder_material_properties, bearing_seat_weld_properties, joist_material_properties, top_chord_connections, bottom_chord_connections)

    inputs = Inputs(design_code, joist_dimensions, chord_dimensions, diagonal_dimensions, bolt_properties, shield_plate_dimensions, diagonal_sections, diagonal_bracing, bearing_seat_dimensions, chord_splice_dimensions,  girder_dimensions, girder_material_properties, bearing_seat_weld_properties, joist_material_properties, top_chord_connections, bottom_chord_connections)



    ####calculate section properties 

    #chord section properties
    chord_section_properties = Properties.calculate_chord_section_properties(chord_dimensions)

    #diagonal section properties
    diagonal_section_geometry = [Geometry.define_diagonal_cross_section_geometry(diagonal_dimensions.B[i], diagonal_dimensions.H[i], diagonal_dimensions.R[i], diagonal_dimensions.t[i]) for i in eachindex(diagonal_dimensions.t)]

    diagonal_section_properties = [CrossSection.Properties.open_thin_walled(diagonal_section_geometry[i].center, fill(diagonal_dimensions.t[i], length(diagonal_section_geometry[i].center)-1)) for i in eachindex(diagonal_dimensions.t)]


    properties = Properties.Section(chord=chord_section_properties, diagonals=diagonal_section_properties)


    #define joist geometry 
    coordinates = Model.define_joist_model_coordinates(joist_dimensions, chord_section_properties, bearing_seat_dimensions, girder_dimensions, chord_dimensions)

    element_connectivity, elements_by_component, element_labels = Model.define_joist_model_element_connectivity(coordinates)

    elements_by_section, element_sections = Model.define_element_types(elements_by_component, diagonal_sections)

    model = Model.run_joist_analysis_model(coordinates, element_connectivity, elements_by_component, element_labels, diagonal_dimensions, diagonal_section_properties, chord_section_properties, diagonal_sections)


    ########
    #### calculate diagonal strengths 

    #tension 

    diagonal_tensile_strength = Strength.calculate_all_diagonal_tensile_strengths(diagonal_sections, diagonal_dimensions, diagonal_section_properties, joist_material_properties, bolt_properties.hole_diameter, design_code)


    #compression 

    #diagonal cross-section local buckling properties 
    # Pcrℓ_diagonal, diagonal_section_local_buckling = Properties.calculate_diagonal_local_buckling_load(diagonal_section_geometry, diagonal_dimensions, joist_material_properties)

    Pcrℓ_diagonal = deserialize("/Users/crismoen/.julia/dev/OpenWebJoist/assets/diagonal_Pcrl")
    diagonal_section_local_buckling = deserialize("/Users/crismoen/.julia/dev/OpenWebJoist/assets/diagonal_section_local_buckling")


    #diagonal global buckling 
    Pcre_diagonal, diagonal_global_buckling = Strength.calculate_diagonal_global_buckling(model, diagonal_sections, diagonal_bracing, diagonal_dimensions, diagonal_section_geometry, joist_material_properties)


    #diagonal compressive strength 
    diagonal_compressive_strength = Strength.calculate_all_diagonal_compressive_strengths(diagonal_sections, diagonal_dimensions, Pcre_diagonal, Pcrℓ_diagonal, diagonal_section_properties, joist_material_properties.fy, design_code)

    ####chord strength
    chord_compressive_strength = Strength.calculate_chord_compressive_strength(chord_section_properties, joist_material_properties, design_code, chord_dimensions)

    #chord tensile strength 
    chord_tensile_strength = Strength.calculate_chord_tensile_strength(chord_section_properties, chord_dimensions, joist_material_properties, design_code)

    #chord splice 
    chord_splice_strength = Strength.calculate_chord_splice_strength(bolt_properties.diameter, bolt_properties.Fnv, joist_material_properties.fu, design_code, chord_splice_dimensions.t, chord_dimensions.t)

    #bearing seat weld strength 
    bearing_seat_weld_strength = Strength.calculate_bearing_seat_weld_strength(bearing_seat_weld_properties, bearing_seat_dimensions, joist_material_properties, girder_dimensions, girder_material_properties, design_code)

    #bearing seat compression field strength
    bearing_seat_compression_field_strength = Strength.calculate_bearing_seat_compression_field_strength(joist_material_properties, bearing_seat_dimensions, design_code)

    #bearing seat top chord connection 
    num_bolts = 4
    bearing_seat_chord_connection_strength = Strength.calculate_bolted_connection_strength(bolt_properties.diameter, bolt_properties.Fnv, joist_material_properties.fu, design_code, bearing_seat_dimensions.t, chord_dimensions.t, num_bolts)

    #calculate top chord connection strengths 
    top_chord_connection_strengths = Strength.calculate_connection_strength_series(top_chord_connections, diagonal_sections, diagonal_dimensions, bolt_properties, joist_material_properties, design_code, bearing_seat_dimensions, shield_plate_dimensions, chord_dimensions)

    #calculate bottom chord connection strengths 
    bottom_chord_connection_strengths = Strength.calculate_connection_strength_series(bottom_chord_connections, diagonal_sections, diagonal_dimensions, bolt_properties, joist_material_properties, design_code, bearing_seat_dimensions, shield_plate_dimensions, chord_dimensions)



    #calculate diagonal demands
    diagonal_elements = findall(element->element==1, occursin.("diagonal", model.inputs.element.cross_section))
    diagonal_demands = [model.solution.forces[diagonal_elements[i]][7] for i in eachindex(diagonal_elements)]

    #calculate diagonal strength 
    eRn_diagonals = Vector{Float64}(undef, length(diagonal_elements))
    for i in eachindex(diagonal_elements)

        if diagonal_demands[i] > 0.0  #tension 

            eRn_diagonals[i] = diagonal_tensile_strength[i].eTn

        else diagonal_demands[i] < 0.0

            eRn_diagonals[i] = diagonal_compressive_strength[i].ePnℓ

        end

    end

    diagonal_DC = calculate_demand_to_capacity(abs.(diagonal_demands)/1000, eRn_diagonals)

    eRn_top_chord_connections = [top_chord_connection_strengths[i].eRn for i in eachindex(top_chord_connection_strengths)]
    top_chord_connection_DC = calculate_demand_to_capacity(abs.(diagonal_demands)/1000, eRn_top_chord_connections)

    eRn_bottom_chord_connections = [bottom_chord_connection_strengths[i].eRn for i in eachindex(bottom_chord_connection_strengths)]
    bottom_chord_connection_DC = calculate_demand_to_capacity(abs.(diagonal_demands)/1000, eRn_bottom_chord_connections)


    eRn_top_chord = chord_compressive_strength.ePn
    top_chord_elements = findall(element->element==1, occursin.("top chord", model.inputs.element.cross_section))
    top_chord_compression_demand = [model.solution.forces[top_chord_elements[i]][7] for i in eachindex(top_chord_elements)]
    top_chord_compression_DC = calculate_demand_to_capacity(abs.(top_chord_compression_demand)./1000, fill(eRn_top_chord, length(top_chord_compression_demand)))

    eRn_bottom_chord = chord_tensile_strength.eTn
    bottom_chord_elements = findall(element->element==1, occursin.("bottom chord", model.inputs.element.cross_section))
    bottom_chord_tension_demand = [model.solution.forces[bottom_chord_elements[i]][7] for i in eachindex(bottom_chord_elements)]
    bottom_chord_tension_DC = calculate_demand_to_capacity(abs.(bottom_chord_tension_demand)./1000, fill(eRn_bottom_chord, length(bottom_chord_tension_demand)))

    eRn_chord_splice = chord_splice_strength.eRn

    # top_chord_elements = findall(element->element=="top chord", model.inputs.element.cross_section)
    # chord_element_index = findall(location->location<chord_splice_dimensions.location, cumsum(model.properties.L[top_chord_elements]))
    
    # if isempty(chord_element_index)
    #     top_chord_splice_DC = calculate_demand_to_capacity(0.0, eRn_chord_splice)  #no splice 
    # else
    #     top_chord_splice_DC = calculate_demand_to_capacity(abs.(top_chord_compression_demand[chord_element_index[end]])./1000, eRn_chord_splice)
    # end

    bottom_chord_splice_DC = calculate_chord_splice_DC("bottom chord", model, chord_splice_dimensions.location.bottom, eRn_chord_splice, bottom_chord_tension_demand, joist_dimensions)
    top_chord_splice_DC = calculate_chord_splice_DC("top chord", model, chord_splice_dimensions.location.top, eRn_chord_splice, top_chord_compression_demand, joist_dimensions)


    eRn_bearing_seat_weld = bearing_seat_weld_strength.eRn
    bearing_seat_weld_DC_1 = calculate_demand_to_capacity(abs.(model.solution.reactions[1][1])./1000, eRn_bearing_seat_weld)
    bearing_seat_weld_DC_2 = calculate_demand_to_capacity(abs.(model.solution.reactions[2][1])./1000, eRn_bearing_seat_weld)

    top_of_girder_elements = findall(element->element==1, occursin.("top of girder", model.inputs.element.cross_section))
    eRn_bearing_compression_field = bearing_seat_compression_field_strength.eRn
    bearing_seat_compression_field_DC_1 = calculate_demand_to_capacity(abs.(model.solution.forces[top_of_girder_elements[1]][1])./1000, eRn_bearing_compression_field)
    bearing_seat_compression_field_DC_2 = calculate_demand_to_capacity(abs.(model.solution.forces[top_of_girder_elements[2]][1])./1000, eRn_bearing_compression_field)

    eRn_bearing_seat_chord_connection = bearing_seat_chord_connection_strength.eRn  
    bearing_seat_chord_connection_DC_1 = calculate_demand_to_capacity(abs.(model.solution.forces[top_chord_elements[2]][1])./1000, eRn_bearing_seat_chord_connection)
    bearing_seat_chord_connection_DC_2 = calculate_demand_to_capacity(abs.(model.solution.forces[top_chord_elements[end-1]][1])./1000, eRn_bearing_seat_chord_connection)



    bearing_seat_weld_demands = [abs.(model.solution.reactions[1][1])./1000, abs.(model.solution.reactions[2][1])./1000]
    bearing_seat_compression_field_demands = [abs.(model.solution.forces[top_of_girder_elements[1]][1])./1000, abs.(model.solution.forces[top_of_girder_elements[2]][1])./1000]
    bearing_seat_chord_connection_demands = [abs.(model.solution.forces[top_chord_elements[2]][1])./1000, abs.(model.solution.forces[top_chord_elements[end-1]][1])./1000]
    demand = ComponentDemands(diagonal_demands, top_chord_compression_demand, bottom_chord_tension_demand, bearing_seat_weld_demands, bearing_seat_compression_field_demands, bearing_seat_chord_connection_demands)


    strength = Strength.Components(diagonal_tensile_strength, diagonal_compressive_strength, chord_compressive_strength, chord_tensile_strength, chord_splice_strength, bearing_seat_weld_strength, bearing_seat_compression_field_strength, bearing_seat_chord_connection_strength, top_chord_connection_strengths, bottom_chord_connection_strengths, eRn_diagonals, eRn_top_chord_connections, eRn_bottom_chord_connections, eRn_top_chord, eRn_bottom_chord, eRn_chord_splice, eRn_bearing_seat_weld, eRn_bearing_compression_field, eRn_bearing_seat_chord_connection)


    Δ_joist = maximum(abs.([model.solution.displacements[i][2] for i in eachindex(model.solution.displacements)]))

    argmax(abs.([model.solution.displacements[i][2] for i in eachindex(model.solution.displacements)]))

    joist_deflection_DC = calculate_demand_to_capacity(Δ_joist, joist_dimensions.span_length/180)

    DC_all = [diagonal_DC.max_DC, top_chord_connection_DC.max_DC, bottom_chord_connection_DC.max_DC, top_chord_compression_DC.max_DC, bottom_chord_tension_DC.max_DC, bottom_chord_splice_DC.max_DC, top_chord_splice_DC.max_DC, bearing_seat_weld_DC_1.max_DC, bearing_seat_weld_DC_2.max_DC, bearing_seat_compression_field_DC_1.max_DC, bearing_seat_compression_field_DC_2.max_DC, bearing_seat_chord_connection_DC_1.max_DC, bearing_seat_chord_connection_DC_2.max_DC, joist_deflection_DC.max_DC]

    DC_labels = ["diagonal", "top chord to diagonal connection", "bottom chord to diagonal connection", "top chord compression", "bottom chord tension", "bottom chord splice", "top chord splice", "bearing seat weld end 1", "bearing seat weld end 2", "bearing seat compression field end 1", "bearing seat compression field end 2", "bearing seat chord connection end 1", "bearing seat chord connection end 2", "joist deflection"]

    max_DC = maximum(DC_all)

    joist_strength = 1.0/max_DC  #lbs/ft 
    controlling_DC = DC_labels[argmax(DC_all)]

    if controlling_DC == "diagonal"

        failure_location = diagonal_DC.max_location

    elseif controlling_DC == "top chord to diagonal connection"

        failure_location = top_chord_connection_DC.max_location

    elseif controlling_DC == "bottom chord to diagonal connection"

        failure_location = bottom_chord_connection_DC.max_location

    elseif controlling_DC == "top chord compression"

        failure_location = top_chord_compression_DC.max_location  

    elseif controlling_DC == "bottom chord tension"

        failure_location = bottom_chord_tension_DC.max_location 

    elseif controlling_DC == "bottom chord splice"

        failure_location = bottom_chord_splice_DC.max_location 

    elseif controlling_DC == "top chord splice"

        failure_location = top_chord_splice_DC.max_location
        
    else

        failure_location = 0

    end

    limit_state_loads= [(DC_labels[i], 1 ./DC_all[i]) for i in eachindex(DC_labels)]

    demand_to_capacity = DemandToCapacity(diagonal_DC, top_chord_connection_DC, bottom_chord_connection_DC, top_chord_compression_DC, bottom_chord_tension_DC, top_chord_splice_DC, bottom_chord_splice_DC, bearing_seat_weld_DC_1, bearing_seat_weld_DC_2, bearing_seat_compression_field_DC_1, bearing_seat_compression_field_DC_2, bearing_seat_chord_connection_DC_1, bearing_seat_chord_connection_DC_2, joist_deflection_DC, DC_all, DC_labels, max_DC, controlling_DC, failure_location)

    joist_span = JoistSpan(inputs, properties, model, strength, Δ_joist, demand, demand_to_capacity, limit_state_loads, joist_strength)

    return joist_span

end




end # module
