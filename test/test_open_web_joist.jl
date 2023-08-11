
using OpenWebJoist

span_length = 50.0 * 12.0;
 
joist = OpenWebJoist.Geometry.JoistDimensions(span_length = span_length, depth = 35.44, connection_node_locations = (top_chord=[], bottom_chord=[]), bottom_chord_offset=[XX, XX], top_of_girder_offset=[XX, XX]);

# chord 
chord = OpenWebJoist.Geometry.ChordDimensions(name = ["top chord", "bottom chord"], t=[0.061, 0.061], centerline_cross_section_coordinates = [XXXX, XXXX])

# diagonal cross sections
diagonal = OpenWebJoist.Geometry.DiagonalDimensions(name=["16 gauge diagonal", "14 gauge diagonal", "12 gauge diagonal"], t=[0.061, 0.078, 0.105], centerline_cross_section_coordinates = [XXXX, XXXX, XXXX]);

# chord connections
connection = (name =[XX, XX], strength = [XXXXX])

#assign chord sections
chord_sections = (top = "top chord", bottom = "bottom chord")

#assign diagonal sections 
diagonal_sections = [["12 gauge diagonal", "12 gauge diagonal"]; fill("12 gauge diagonal", num_diagonals - 4); ["12 gauge diagonal", "12 gauge diagonal"]];

# define top chord connections 
top_chord_connections = ["bearing seat"; fill("reinforced connection", num_diagonals-2); "bearing seat"];

# define bottom chord connections 
bottom_chord_connections = ["bearing seat"; fill("reinforced connection", num_diagonals-2); "bearing seat"];

# joist material properties 
joist_material_properties = OpenWebJoist.Properties.JoistMaterial(fy=50.0, fu=65.0, E=29500.0, ν=0.30);


# Calculate span limit state loads.
joist_span = OpenWebJoist.evaluate_joist_span(design_code, joist_dimensions, chord_dimensions, diagonal_dimensions, bolt_properties, shield_plate_dimensions, diagonal_sections, diagonal_bracing, bearing_seat_dimensions, chord_splice_dimensions,  girder_dimensions, girder_material_properties, bearing_seat_weld_properties, joist_material_properties, top_chord_connections, bottom_chord_connections);
