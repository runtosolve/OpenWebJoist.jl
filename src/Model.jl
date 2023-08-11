module Model

using InstantFrame

using ..Properties, ..Geometry 

function run_joist_analysis_model(coordinates, element_connectivity, elements_by_component, element_labels, diagonal_dimensions, diagonal_section_properties, chord_section_properties, diagonal_section_assignments)

    #####define nodes 
    joist_nodes = define_joist_model_nodes(coordinates)

    # element_connectivity, elements_by_component = define_joist_model_element_connectivity(coordinates)

    #define element types 
    elements_by_section, element_sections = Model.define_element_types(elements_by_component, diagonal_section_assignments)

    material = InstantFrame.Material(names=["steel"], E=[29500000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  ##ρ = lbs * s^2 / in^4

    cross_section = InstantFrame.CrossSection(names=[diagonal_dimensions.name; ["top chord", "bottom chord", "top of girder", "rigid"]], A=[[diagonal_section_properties[i].A for i in eachindex(diagonal_section_properties)]; [chord_section_properties.A, chord_section_properties.A, 10.0E1, 10.0E1]], Iy=[[diagonal_section_properties[i].Iyy for i in eachindex(diagonal_section_properties)]; [chord_section_properties.Iyy, chord_section_properties.Iyy, 10.0E1, 10.0E1]], Iz=[[diagonal_section_properties[i].Ixx for i in eachindex(diagonal_section_properties)]; [chord_section_properties.Ixx, chord_section_properties.Ixx, 10.0E1, 10.0E1]], J=[[diagonal_section_properties[i].A for i in eachindex(diagonal_section_properties)]; [chord_section_properties.J, chord_section_properties.J, 10.0E1, 10.0E1]])

    connection = InstantFrame.Connection(names=["rigid", "pinned"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, 0.0], rz=[Inf, 0.0]))

    node = InstantFrame.Node(numbers=1:size(joist_nodes, 1), coordinates=joist_nodes)

    element_connections = [fill(("pinned", "pinned"), length(elements_by_section.diagonals)); fill(("rigid","rigid"), size(element_connectivity, 1) - length(elements_by_section.diagonals))]

    #element_connections = [fill(("rigid", "rigid"), length(elements_by_section.diagonals)); fill(("rigid","rigid"), size(element_connectivity, 1) - length(elements_by_section.diagonals))]


    # element = InstantFrame.Element(numbers=1:size(element_connectivity, 1), nodes=element_connectivity, orientation=zeros(Float64, size(element_connectivity, 1)), connections=element_connections, cross_section=element_sections, material=fill("steel", size(element_connectivity, 1)), types=fill("frame", size(element_connectivity, 1)))

    element = InstantFrame.Element(numbers=1:size(element_connectivity, 1), nodes=element_connectivity, orientation=zeros(Float64, size(element_connectivity, 1)), connections=element_connections, cross_section=element_sections, material=fill("steel", size(element_connectivity, 1)), types=element_labels)


    support = InstantFrame.Support(nodes=[node.numbers[end-1], node.numbers[end]], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[0.0,0.0]))

    # support = InstantFrame.Support(nodes=[54, 67], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[0.0,0.0]))


    uniform_load = InstantFrame.UniformLoad(labels=["gravity loads"], elements=findall(type->type=="top chord", element_sections), magnitudes=(qX=fill(0.0, length(elements_by_component.top_chord)), qY=fill(-1.0/12, length(elements_by_component.top_chord)), qZ=fill(0.0, length(elements_by_component.top_chord)), mX=fill(0.0, length(elements_by_component.top_chord)), mY=fill(0.0, length(elements_by_component.top_chord)), mZ=fill(0.0, length(elements_by_component.top_chord))))

    point_load = InstantFrame.PointLoad(nothing)

    analysis_type = "first order"
    model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type)

    return model

end

function define_joist_model_coordinates(joist, chord_section_properties, bearing_seat, girder, chord)

    #define model coordinates
    
    # #determine joist end node spacing
    # joist_ends = define_joist_ends(joist.span_length, joist.node_spacing)

    #joist diagonals 
    diagonal_coordinates = Geometry.calculate_joist_diagonal_coordinates(joist.span_length, joist.depth, joist.bottom_chord_tail_length, joist.node_spacing, bearing_seat, chord)

    #top chord
    top_chord_coordinates = Geometry.define_top_chord_coordinates(chord_section_properties, diagonal_coordinates, joist)

    #bottom chord
    bottom_chord_coordinates = Geometry.define_bottom_chord_coordinates(joist, chord_section_properties, diagonal_coordinates)

    #top of girder
    top_of_girder_coordinates = Geometry.define_top_of_girder_coordinates(bearing_seat, girder, chord, joist)

    #collect all coordinates 
    coordinates = Geometry.Coordinates(diagonals = diagonal_coordinates, top_chord = top_chord_coordinates, bottom_chord = bottom_chord_coordinates, top_of_girder=top_of_girder_coordinates)

    return coordinates

end

    function define_joist_model_nodes(coordinates)

        joist_coordinates = [(coordinates.diagonals.top[i][1], coordinates.diagonals.top[i][2], 0.0) for i in eachindex(coordinates.diagonals.top)]
    
        joist_coordinates = vcat(joist_coordinates, [(coordinates.diagonals.bottom[i][1], coordinates.diagonals.bottom[i][2], 0.0) for i in eachindex(coordinates.diagonals.bottom)])
    
        joist_coordinates = vcat(joist_coordinates, [(coordinates.top_chord[i][1], coordinates.top_chord[i][2], 0.0) for i in eachindex(coordinates.top_chord)])
    
        joist_coordinates = vcat(joist_coordinates, [(coordinates.bottom_chord[i][1], coordinates.bottom_chord[i][2], 0.0) for i in eachindex(coordinates.bottom_chord)])
    
        joist_coordinates = vcat(joist_coordinates, [(coordinates.top_of_girder[i][1], coordinates.top_of_girder[i][2], 0.0) for i in eachindex(coordinates.top_of_girder)])
    
        joist_nodes = joist_coordinates

        return joist_nodes
    
    end

    function define_joist_model_element_connectivity(coordinates)

        #number of nodes 
        num_top_chord_nodes = size(coordinates.top_chord, 1)
        num_bottom_chord_nodes = size(coordinates.bottom_chord, 1)
        num_top_diagonal_nodes = size(coordinates.diagonals.top, 1)
        num_bottom_diagonal_nodes = size(coordinates.diagonals.bottom, 1)
        num_diagonal_nodes = num_top_diagonal_nodes + num_bottom_diagonal_nodes
        num_top_of_girder_nodes = size(coordinates.top_of_girder, 1)
    
        node_num_offsets = cumsum([num_top_diagonal_nodes, num_bottom_diagonal_nodes, num_top_chord_nodes, num_bottom_chord_nodes, num_top_of_girder_nodes])
    
        #number of elements 
        num_top_chord_elements = num_top_chord_nodes - 1
        num_bottom_chord_elements = num_bottom_chord_nodes - 1
        num_diagonal_elements = size(coordinates.diagonals.bottom, 1)
        # num_top_of_girder_elements = 2
        # elem_num_offsets = cumsum([num_diagonal_elements, num_top_chord_elements, num_bottom_chord_elements, num_top_of_girder_elements]) 
    
        offset = node_num_offsets[1]
        diagonal_elements = [(i, i+offset) for i=1:num_diagonal_elements]
        offset = node_num_offsets[2]
        top_chord_elements = [(i, i+1) .+ offset for i=1:num_top_chord_elements]
        offset = node_num_offsets[3]
        bottom_chord_elements = [(i, i+1) .+ offset for i=1:num_bottom_chord_elements]
    
        offset = node_num_offsets[4]+1
        top_of_girder_elements = [(offset, top_chord_elements[1][2]), (offset+1, top_chord_elements[end][1])]
    
        # top_of_girder_elements = [top_of_girder_elements; (diagonal_elements[1][1], offset); (diagonal_elements[end][1], offset+1)]
    
    
        #rigid links 
    
        #first chord to diagonal link 
        rigid_elements_bearing_seat_1 = [(top_chord_elements[1][2], diagonal_elements[1][1])]

       
    
        #add all typical top chord connections 
        num_typical_top_chord_connections = num_top_chord_nodes - 2 - 2
        diagonal_element_range = [[i, i+1] for i=2:2:num_diagonal_elements-1]
        rigid_elements_top_chord = reduce(vcat, [[(top_chord_elements[i+2][1], diagonal_elements[diagonal_element_range[i][1]][1]), (top_chord_elements[i+2][1], diagonal_elements[diagonal_element_range[i][2]][1])] for i=1:num_typical_top_chord_connections])
    
        #add last top chord to diagonal link
        rigid_elements_bearing_seat_2 = [(top_chord_elements[end][1], diagonal_elements[end][1])]
    
        rigid_elements_top_chord = [rigid_elements_bearing_seat_1; rigid_elements_top_chord; rigid_elements_bearing_seat_2]
    
        #bottom chord
    
        #first two rigid links at bearing seat piece in bottom chord
        rigid_elements_bottom_chord_seat_1 = [(bottom_chord_elements[1][2], diagonal_elements[1][2]), (bottom_chord_elements[1][2], diagonal_elements[2][2])]
    
        # num_typical_bottom_chord_connections = num_top_chord_nodes - 2 - 2
        # diagonal_element_range = [[i, i+1] for i=3:2:num_diagonal_elements-1]
        diagonal_element_range = [[i, i+1] for i=3:2:num_diagonal_elements-1]
        rigid_elements_bottom_chord = [[(bottom_chord_elements[i+2][1], diagonal_elements[diagonal_element_range[i][1]][2]), (bottom_chord_elements[i+2][1], diagonal_elements[diagonal_element_range[i][2]][2])] for i=1:num_typical_top_chord_connections]
    
        rigid_elements_bottom_chord_seat_2 = [(bottom_chord_elements[end][1], diagonal_elements[end-1][2]), (bottom_chord_elements[end][1], diagonal_elements[end][2])]
    
        # rigid_elements_bottom_chord = reduce(vcat, [rigid_elements_bottom_chord_seat_1; rigid_elements_bottom_chord; rigid_elements_bottom_chord_seat_2])
        rigid_elements_bottom_chord = reduce(vcat, [rigid_elements_bottom_chord_seat_1; rigid_elements_bottom_chord])
    
        rigid_elements = [rigid_elements_top_chord; rigid_elements_bottom_chord]
    
        # num_rigid_elements = size(rigid_elements, 1)
    
    
        elements_by_component = (diagonals = diagonal_elements, top_chord=top_chord_elements, bottom_chord = bottom_chord_elements, top_of_girder=top_of_girder_elements, rigid_elements = rigid_elements)

        element_labels = [fill("diagonals", length(diagonal_elements)); fill("top chord", length(top_chord_elements)); fill("bottom chord", length(bottom_chord_elements)); fill("top of girder", length(top_of_girder_elements)); fill("rigid connections", length(rigid_elements))]

        element_connectivity = [diagonal_elements; top_chord_elements; bottom_chord_elements; top_of_girder_elements; rigid_elements]
    
        return element_connectivity, elements_by_component, element_labels
    
    end

    function define_element_types(elements_by_component, diagonal_section_assignments)


        diagonal_sections = diagonal_section_assignments
        top_chord_sections = fill("top chord", length(elements_by_component.top_chord))
        bottom_chord_sections = fill("bottom chord", length(elements_by_component.bottom_chord))
        top_of_girder_sections = fill("top of girder", length(elements_by_component.top_of_girder))
        rigid_element_sections = fill("rigid", length(elements_by_component.rigid_elements))
    
        elements_by_section = (diagonals = diagonal_sections, top_chord = top_chord_sections, bottom_chord=bottom_chord_sections, top_of_girder = top_of_girder_sections, rigid_elements = rigid_element_sections)
    
        element_sections = [diagonal_sections; top_chord_sections; bottom_chord_sections; top_of_girder_sections; rigid_element_sections]
    
        return elements_by_section, element_sections
    
    end



end  #module 