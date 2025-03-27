from opentrons import protocol_api
import csv
import io
import re

metadata = {
    "protocolName": "WTA Assembly Validation",
    "description": "DNA WTA Computer Assembly Validation Assay",
    "author": "Yaman Al Janaideh <yaman.aljanaideh@mail.mcgill.ca>"
}
requirements = {"robotType": "OT-2", "apiLevel": "2.20"}

def load_bundled_data():
    import os
    """
    Load bundled CSV data for local simulation. It searches for a file matching the pattern 'SDR_computer*.csv'
    in the current directory and returns the file contents and its name.

    Returns:
        tuple: A dictionary (`bundled_data`) where the key is the file name and the value is the file content, 
               and the name of the file (`file_local_name`).

    Raises:
        ValueError: If no file matching the pattern 'SDR_computer*.csv' is found in the directory.
    """
    # Define the pattern to search for SDR_computer CSV files
    pattern = re.compile(r'SDR_computer.*\.csv')
    file_name = None

    # Step 1: Search for a file matching the pattern in the current directory
    for file in os.listdir('.'):
        if pattern.match(file):
            file_name = file  # Found a matching file, store its name
            break
    
    # Step 2: Raise an error if no matching file is found
    if not file_name:
        raise ValueError("No valid SDR_computer file found in the current directory.")
    
    # Step 3: Load the content of the found file
    bundled_data = {}
    with open(file_name, 'r') as file:
        file_local_name = os.path.basename(file_name)  # Get the file's base name
        bundled_data[file_local_name] = file.read()    # Store the file's content in the dictionary

    # Step 4: Return the bundled data and the file name
    return bundled_data, file_local_name

def calculate_liquid_height(liquid_volume, tube_type):
    """
    Calculate the height of a liquid in a tube based on its volume and the type of tube.

    Args:
        liquid_volume (float): The volume of the liquid in microliters (µL).
        tube_type (str): The type of tube, either 'buffer' (50 mL Falcon) or 'liquid' (1.5 mL Eppendorf).

    Returns:
        float: The calculated liquid height in millimeters (mm).

    Raises:
        ValueError: If the liquid volume exceeds the capacity of the tube.
    """
    # Define the cone and cylinder dimensions based on the tube type
    if tube_type == "buffer":
        cone_height = 3.0  # cm, height of the conical bottom of the Falcon tube
        cone_radius = 0.75  # cm, radius of the conical bottom
        cylinder_height = 12  # cm, height of the cylindrical section of the Falcon tube
        cylinder_radius = 1.5  # cm, radius of the cylindrical section
        offset = 10  # Correction offset for Falcon tubes
    elif tube_type == "liquid":
        cone_height = 0.5  # cm, height of the conical bottom of the Eppendorf tube
        cone_radius = 0.5  # cm, radius of the conical bottom
        cylinder_height = 1.5  # cm, height of the cylindrical section of the Eppendorf tube
        cylinder_radius = 0.6  # cm, radius of the cylindrical section
        offset = 2.5  # Correction offset for Eppendorf tubes
    else:
        raise ValueError("Unknown tube type. Valid types are 'buffer' or 'liquid'.")

    # Convert the liquid volume from µL to mL (1 cm³ = 1 mL)
    liquid_volume = liquid_volume / 1000

    # Calculate the volume of the conical section (V_cone = 1/3 * pi * r^2 * h)
    V_cone = (1 / 3) * 3.14 * (cone_radius ** 2) * cone_height

    # If the liquid volume fits entirely within the conical section
    if liquid_volume <= V_cone:
        # Calculate the height of the liquid in the cone (h_cone = 3V / (πr²))
        height_in_cone = (3 * liquid_volume) / (3.14 * (cone_radius ** 2))
        return (height_in_cone * 10) - offset  # Convert from cm to mm and apply offset

    # Otherwise, calculate the height in the cylindrical section
    else:
        # Volume left to fill in the cylindrical section
        volume_in_cylinder = liquid_volume - V_cone
        # Calculate the height in the cylindrical section (h_cylinder = V / (πr²))
        height_in_cylinder = volume_in_cylinder / (3.14 * (cylinder_radius ** 2))

        # If the calculated height exceeds the cylinder height, raise an error
        if height_in_cylinder > cylinder_height:
            raise ValueError("Liquid volume exceeds tube capacity.")

        # Return the total height (cone height + cylindrical height), converted to mm and with the offset applied
        return ((cone_height + height_in_cylinder) * 10) - offset  # Convert from cm to mm and apply offset
    
def load_computer(number_of_memories, conditions, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack, liquid_eppendorf_rack_2=None, file_local_name=None):
    """
    Loads and parses the liquid data from the CSV file for the protocol, organizing it based on the assay conditions and provided number of memories into components, 
    buffers, and possibly annihilators or inputs, depending on the assay.

    Args:
        number_of_memories (int): The number of memories (full networks, each with a unique input for its corresponding memory) to plate (1-4).
        conditions (dict): A dictionary containing assay conditions {condition_number (int): [list of liquids in condition (str)]}.
        bundled_data (dict): A dictionary containing file data, {file_name (str): file_data (bytes or text)}.
        liquid_eppendorf_rack_1 (Labware): Labware object representing the first rack of liquid tubes.
        buffer_tube_rack (Labware): Labware object representing the rack holding the buffer tubes.
        liquid_eppendorf_rack_2 (Labware, optional): Labware object representing the second rack of liquid tubes (optional). Defaults to None.
        file_local_name (str, optional): Local name of the CSV file if bundled_data contains decoded text. Defaults to None.

    Returns:
        tuple: Four dictionaries containing the following:
                - inputs (dict):  {memopry_number: {srouce: (str), total_volume: (float), reaction_volume: (float)}}
                - components (dict): {memory_number: {liquid_type: {source: (str), total_volume: (float), reaction_volume: (float)}}}
                - buffer (dict): {source: (str), total_volume: (float), memory_number: {condition_number: total buffer volume (float)}}
                - annihilators (list): [{source: (str), total_volume: (float), reaction_volume: (float)}]

                    where:
                        from CSV file:
                            - memory_number (int or str): The memory number corresponding to the liquid components (0 for annihlation gates and 'b' for buffer).
                            - liquid_type: The type of liquid component (e.g., input, buffer, etc.).
                            - source: The location of the corresponding tube in reference to the appropriate rack.
                            - total_volume: The total volume of the liquid component in the tube.
                        Standardized:
                            - total_reaction_volume: The total volume of the liquid component required for the reaction.
                            - reaction_volume: The volume of the liquid component required for the reaction.

    Raises:
        ValueError: If no valid SDR_computer file is found in the provided data.
    """
    # Step 1: Find and read the SDR_computer CSV file
    pattern = re.compile(r'SDR_computer.*\.csv') # Pattern to identify the correct file in bundled_data
    file_name = None
    
    # Look for the file in bundled_data
    for key in bundled_data.keys():
        if pattern.match(key):
            file_name = key
            break
    
    # If no valid file is found, raise an error
    if not file_name:
        raise ValueError("No valid SDR_computer file found in the provided data.")
    
    # Extract CSV data from bundled_data; check if it needs to be decoded
    liquids_file = bundled_data[file_name].decode('utf-8') if isinstance(bundled_data[file_name], bytes) else bundled_data[file_local_name]
    
    # Step 2: Initialize structures to store parsed data
    components = {}        # Store liquid components of each memory (e.g., strands and gates)
    buffer = {}            # Store buffer data
    inputs = {}            # Store input data if test type is 1 (Multiplication)
    annihilators = []      # Store annihilators if test type is 3 (Annihilation)
    
    total_reaction_volume = 1000.0 # Total reaction volume per well
    standard_liquid_volume = 10.0  # Standard volume for liquid components
    
    # Step 3: Parse CSV data and organize it into components, inputs, buffers, etc.
    liquids_reader = csv.DictReader(io.StringIO(liquids_file))
    row_count = 0 # Keep track of the row number to handle conditions for test 1
    
    for row in liquids_reader:
        row_count += 1
        
        # Get memory number and determine its type (int or str)
        memory_number = int(row["memory_number"]) if isinstance(row["memory_number"], int) else str(row["memory_number"])
        
        # Assign the correct rack for the liquid
        if isinstance(memory_number, str) and memory_number == "b": 
            liquid_rack = buffer_tube_rack # Buffer memory gets the buffer rack
        elif row_count <= 23:   # 2 memories per rack (all anhilation on rack 1)
            liquid_rack = liquid_eppendorf_rack_1 
        elif row_count > 23: 
            liquid_rack = liquid_eppendorf_rack_2
        
        # Extract relevant data from the row    
        liquid_type = row['liquid_type']
        source = row['liquid_source'] 
        liquid_total_volume = float(row['liquid_total_volume'])
        
        if liquid_type == "input": 
            inputs[memory_number] = {
                "source": liquid_rack[source],
                "total_volume": liquid_total_volume,
                "reaction_volume": standard_liquid_volume
            }
            continue
        
        if liquid_type == "annihilation":
            annihilators.append({
                "source": liquid_rack[source],
                "total_volume": liquid_total_volume,
                "reaction_volume": standard_liquid_volume
            })
            continue
        
        if liquid_type == "buffer":
            buffer["source"] = liquid_rack[source],    
            buffer["total_volume"] = liquid_total_volume
            continue
        
        if memory_number not in components: 
            components[memory_number] = {}
            
        components[memory_number][liquid_type] = {
            "source": liquid_rack[source],
            "total_volume": liquid_total_volume,
            "reaction_volume": standard_liquid_volume
        }
    
    # Step 4: Calculate buffer volumes for each condition    
    for memory in inputs.keys():
        for condition, liquids_in_condition in conditions.items():
            # calculate buffer volume for each condition; for condition 3 --> will add all components of all memories except input for 1 memory and annihilation treated as a separate component (-1 -1 = -2)
            num_annihilators = len(annihilators)
            num_liquids = ((len(liquids_in_condition)-2)*number_of_memories + 1 + num_annihilators) if condition == 3 else len(liquids_in_condition)
            buffer_well_volume = total_reaction_volume - (num_liquids * standard_liquid_volume)
            
            if memory not in buffer: 
                buffer[memory] = {}
            buffer[memory][condition] = buffer_well_volume

    # Step 5: Return the dicts containing the parsed data
    return inputs, components, annihilators, buffer
      
def map_plate(conditions, plate, num_replicates, inputs, buffer):
    """
    Maps each memory and its corresponding conditions to specific wells on a 48-well or 96-well plate, 
    assigning destinations for liquids and buffer volumes.

    Args:
        conditions (dict): A dictionary where keys are condition numbers and values are lists of liquid types.
        plate (Labware): The 48-well plate used in the protocol.
        num_replicates (int): The number of replicates for each condition (1 or 3).
        inputs (dict): A dictionary of input components for each unique memory.
        buffer (dict): A dictionary of buffer volumes for each well to be plated.

    Returns:
        dict: A dictionary (plate_map) that maps each memory and condition to specific well destinations {memory: {condition: {destination: well, liquids: [list], buffer_volume: float}}}.

    Notes:
        - The function loops over all memories and conditions, assigning wells based on memory and condition.
        - It adjusts the layout depending on the number of replicates.
        - Conditions are plated along the columns and the replicates along the rows. Unique memories are plated clockewise in the order they appear in the CSV file.
    """
    plate_map = {} # Initialize the mapping of wells to memories and conditions
    plate_rows = ['A', 'B', 'C', 'D', 'E', 'F'] # Define the rows of the 48-well plate (6 rows)
    plate_cols = list(range(1, 9)) # Define the columns of the 48-well plate (8 columns)

    row_index = 0 # Initialize row index
    col_index = 0 # Initialize column index
    
    memories = inputs.keys()
    
    # Step 1: Loop through each unique memory being plated
    for memory in memories:
        
        # Initialize the entry in the plate map for the current memory
        plate_map[memory] = {}
        
        # Step 2: Loop through each condition and its corresponding liquids
        for condition, liquids_in_condition in conditions.items():
            # Step 3: Assign wells based on the number of replicates (1 or 3)
            if num_replicates == 1:
                # Single replicate: Assign one well for the condition
                plate_map[memory][condition] = {
                    'destination': plate[f"{plate_rows[row_index]}{plate_cols[col_index]}"],
                    'liquids': liquids_in_condition,
                    'buffer_volume': buffer[memory][condition]
                }    
            elif num_replicates == 3:
                # Triple replicates: Assign three consecutive wells in the same column for the condition
                plate_map[memory][condition] = {
                    'destination': [
                        plate[f"{plate_rows[row_index]}{plate_cols[col_index]}"], 
                        plate[f"{plate_rows[row_index + 1]}{plate_cols[col_index]}"],
                        plate[f"{plate_rows[row_index + 2]}{plate_cols[col_index]}"]
                    ],
                    'liquids': liquids_in_condition,
                    'buffer_volume': buffer[memory][condition]
                }
            
            # Step 4: Move to the next column after assigning the wells
            col_index += 1   
            
            # Step 5: Reset column after 6 wells (2 unqiue memories * 3 conditions each) to plate the next memories (3 and/or 4) in the bottom half of the plate
            if col_index > 5:
                col_index = 0
                row_index += num_replicates  # Move to the next row(s) (depending on replicates)
    
    # Step 6: Return the completed plate map                
    return plate_map

def plate_buffer(protocol, p1000_single_left, buffer, plate_map):
    """
    Distributes buffer into the designated wells in the plate based on the buffer volumes for each memory and condition.

    Args:
        protocol (ProtocolContext): The Opentrons protocol context 
        p1000_single_left (InstrumentContext): The left-mount P1000 single-channel pipette used for liquid handling.
        buffer (dict): A dictionary of buffer volumes for each well to be plated {source: (str), total_volume: (float), memory_number: {condition_number: total buffer volume (float)}}.
        plate_map (dict): A dictionary mapping each memory and condition to specific well destinations {memory: {condition: {destination: well, liquids: [list], buffer_volume: float}}}.

    Returns:
        None: The function performs buffer transfers directly on the protocol.

    Notes:
        - Buffer is distributed to wells in descending order of buffer volume for optimal pipette usage.
        - The function aspirates and dispenses in chunks based on the maximum pipette volume.
        - Buffer levels in the source tube are updated as the liquid is transferred.
    """
    protocol.comment("Starting buffer plating.")
    
    # Step 1: Define the source location for the buffer (tube in a rack)
    buffer_source = list(buffer["source"])[0] 
    protocol.comment(f"Buffer source: {buffer_source}")
    # Step 2: Calculate the current liquid height in the buffer tube
    buffer_top = buffer_source.bottom(calculate_liquid_height(buffer["total_volume"], "buffer"))
    
    # Step 3: Prepare the destination wells and buffer volumes from the plate_map
    destinations = []
    volumes = []
    for memory, conditions in plate_map.items():
        for condition, volume in conditions.items():
            wells = plate_map[memory][condition]['destination'] # Get the destination wells for this condition
            buffer_volume = plate_map[memory][condition]['buffer_volume'] # Buffer volume required per well
            
            # Handle multiple wells if replicates exist
            if isinstance(wells, list):
                for well in wells:
                    destinations.append(well)
                    volumes.append(buffer_volume) 
            else:
                destinations.append(wells)
                volumes.append(buffer_volume)
    
    # Step 4: Sort destinations by buffer volume in descending order for optimal pipette usage
    destinations_volumes = sorted(zip(destinations, volumes), key=lambda x: x[1], reverse=True)
    destinations, volumes = zip(*destinations_volumes)

    # Step 5: Pick up tip (plate all wells with this tip as there is no need to change it since the buffer is the same)
    p1000_single_left.pick_up_tip()
    # Mix buffer in the source tube to ensure uniformity
    p1000_single_left.mix(3, 500, buffer_top)
    
    # Step 6: Begin distributing the buffer across the plate, 1 destination well at a time
    for destination in range(0, len(destinations)):
        
        # Aspirate buffer with the pipette and dispense into well
        protocol.comment(f"Aspirating {volumes[destination]} µL of buffer from buffer tube.")
        p1000_single_left.aspirate(volumes[destination], buffer_top)
        p1000_single_left.touch_tip(buffer_source, radius=.05, speed=80, v_offset=-1)
        p1000_single_left.move_to(buffer_source.top(10), speed=50)
        
        # Update buffer volume and recalculate the liquid height in the buffer tube
        new_buffer_volume = buffer["total_volume"] - volumes[destination]
        buffer["total_volume"] = new_buffer_volume
        buffer_top = buffer_source.bottom(calculate_liquid_height(new_buffer_volume, "buffer"))
        
        # Dispense buffer into wells 
        protocol.comment(f"Dispensing {volumes[destination]} µL of buffer into well {destinations[destination]}.")
        well = destinations[destination]
        volume = volumes[destination]
        p1000_single_left.move_to(well.top(1))
        p1000_single_left.move_to(well.bottom(1))
        p1000_single_left.dispense(volume, well.bottom(1))    
        p1000_single_left.touch_tip(well, radius=.05, speed=80, v_offset=-1)
        protocol.delay(1)
        p1000_single_left.blow_out()
    
    # Step 7: Drop tip when the process is complete
    p1000_single_left.drop_tip()
    return    
        
def plate_liquids(protocol, liquid_types, p300_single_right, inputs, annihilators, components, plate_map):
    """
    Distributes liquid components into the designated wells according to the specified protocol.

    The function handles the plating of all components (e.g., reporters, inputs, outputs, annihilators, etc.) 
    into wells defined in the plate_map. It considers the unique requirements for each liquid type, memory, and condition.
    
    Plating logic: by a predefined order of liquid types (liquid_types) to avoid early triggering of gates. 
                The pairs are also in a predefined order to avoid contamination between conditions of a memory 
        
        High level:
                    loop over liquid types in liquid_types:
                        loop over each memory:
                            for memory x:
                                special consideratoins:
                                    1. Reporter: one-to-one mapping for conditions 1 and 2 then one-to-many mapping for condition 3
                                    2. Input: one-to-one mapping for condition 3
                                    3. Output: one-to-one mapping for condition 2

                                - find destinations for the liquid type
                                - new tip --> plate in order of conditions (replicates before moving to next condition)
                                
                    The max possible volumes needed to distribute a memory's liquid type will never be more than 300 uL given 
                    the current setup of assays, so there is no need to change tips until moving to the next memory. Note that
                    this is also means that for any liquid to be plated, it will only need to be aspirated once.
                    
    Args:
        protocol (ProtocolContext): The Opentrons protocol context.
        liquid_types (list): A list of liquid types to be plated ('reporter', 'restoration', 'input', etc.) in a pre-defined order.
        p300_single_right (InstrumentContext): The right-hand P300 single-channel pipette used for liquid handling.
        inputs (dict): A dictionary containing input components for each memory {memopry_number: {srouce: (str), total_volume: (float), reaction_volume: (float)}}.
        annihilators (list): A list of annihilator components shared across all memories  [{source: (str), total_volume: (float), reaction_volume: (float)}].
        components (dict): {memory_number: {liquid_type: {source: (str), total_volume: (float), reaction_volume: (float)}}}
        plate_map (dict): A dictionary mapping each memory and condition to specific well destinations {memory: {condition: {destination: well, liquids: [list], buffer_volume: float}}}.

    Returns:
        None: The function performs the liquid transfers directly on the protocol.
    """
    # Constants
    air_gap = 10  # Air gap in µL to prevent dripping during liquid transfers

    def plate(liquid_to_plate, destinations, volumes):
        """
        Performs the actual pipetting of a liquid to multiple destinations.

        Args:
            liquid_to_plate (dict): Contains 'source', 'total_volume', and 'reaction_volume' for the liquid.
            destinations (list): List of wells to dispense into.
            volumes (list): Corresponding volumes to dispense into each well.
        """
        source = liquid_to_plate['source']
        source_volume = liquid_to_plate['total_volume']
        source_top = source.bottom(calculate_liquid_height(source_volume, "liquid"))
        liquid_to_plate['total_volume'] -= sum(volumes)  # Update total volume after dispensing

        # Pick up a new tip for pipetting
        p300_single_right.pick_up_tip()
        p300_single_right.mix(3, 50, source_top)  # Mix the liquid source before aspiration

        # Aspirate the total volume to be dispensed, plus an air gap to avoid drips
        protocol.comment(f"Aspirating {sum(volumes)} µL from source tube.")
        p300_single_right.aspirate(sum(volumes), source_top)
        p300_single_right.touch_tip(source, radius=0.05, speed=80, v_offset=-1)
        p300_single_right.air_gap(air_gap, height=1)  # Add an air gap to prevent dripping
        
        # Dispense the liquid into each destination well
        for volume, well in zip(volumes, destinations):
            protocol.comment(f"Dispensing {volume} µL into well {well}.")
            p300_single_right.move_to(well.top(1))
            p300_single_right.move_to(well.bottom(1), speed=50)
            p300_single_right.dispense(air_gap + volume, well.bottom(1))  # Dispense the liquid with the air gap
            p300_single_right.touch_tip(well, radius=0.05, speed=80, v_offset=-1)
            protocol.delay(1)  # Short delay for liquid to settle
            if well != destinations[-1]:
                p300_single_right.air_gap(air_gap, height=1)  # Add an air gap between wells to avoid contamination
        p300_single_right.blow_out(destinations[-1].top(1))  # Blow out remaining liquid at the last well
        p300_single_right.drop_tip()  # Discard the tip after use
        return
    
    def find_destinations(liquid, liquid_sources, one_to_one=None, one_to_many=False):
        """
        Determines the list of wells and volumes for a given liquid to be dispensed.
        
        For condition 1:
            1 unique reporter
        For condition 2:
            1 unique reporter and 1 unique output
        For condition 3:
            - all components except input (1 unique input) should go to all memories' condition 3 wells
            - all annihilation gates should go to all memories' condition 3 wells

        Args:
            liquid (str): The type of liquid (e.g., 'input', 'output', 'reporter', etc.).
            liquid_sources (dict): The source data for the liquid.
            one_to_one (int, optional): If specified, only consider this memory.
            one_to_many (bool, optional): If True, handles special cases like condition 3 for reporters.

        Returns:
            tuple: (destinations, volumes) - A list of wells and corresponding volumes for the liquid to be dispensed.
        """
        destinations = []
        volumes = []
        
        if one_to_one:
            # One-to-one mapping: A single memory's liquid goes to its corresponding wells
            memory = one_to_one
            for condition, mapping in plate_map[memory].items():
                if (liquid in mapping['liquids'] and liquid in ["reporter", "output"] and condition != 3) or (liquid in mapping['liquids'] and liquid == 'input'):
                    wells = mapping['destination']
                    reaction_volume = liquid_sources[memory]['reaction_volume'] if liquid == "input" else liquid_sources[memory][liquid]['reaction_volume']
                    if isinstance(wells, list):
                        for well in wells:
                          destinations.append(well)
                          volumes.append(reaction_volume)
                    else:
                        destinations.append(wells)
                        volumes.append(reaction_volume)
            return destinations, volumes
        
        if one_to_many:
            # One-to-many mapping: Liquid goes to all condition 3 wells across all memories
            for memory in plate_map.keys():
                if liquid in plate_map[memory][3]['liquids']:
                    mapping = plate_map[memory][3]
                    wells = mapping['destination']
                    reaction_volume = liquid_sources['reaction_volume'] if liquid == "annihilation" else liquid_sources[memory][liquid]['reaction_volume']
                    if isinstance(wells, list):
                        for well in wells:
                          destinations.append(well)
                          volumes.append(reaction_volume)
                    else:
                        destinations.append(wells)
                        volumes.append(reaction_volume)
            return destinations, volumes

    # Plate liquids in specified order
    for liquid in liquid_types:
        
        if liquid == "input":
            # Input goes to its corresponding wells in condition 3 (one-to-one mapping)
            for memory in inputs:
                destinations, volumes = find_destinations(liquid, inputs, one_to_one=memory)
                source = inputs[memory]["source"]
                protocol.comment(f"Plating {liquid} for memory {memory} from srouce: {source} to {destinations}.")
                plate(inputs[memory], destinations, volumes)
            continue
        
        if liquid == "output":
            # Output goes to its corresponding wells in condition 2 (one-to-one mapping)
            for memory in components.keys():
                destinations, volumes = find_destinations(liquid, components, one_to_one=memory)
                source = components[memory][liquid]["source"]
                protocol.comment(f"Plating {liquid} for memory {memory} from srouce: {source} to {destinations}.")
                plate(components[memory][liquid], destinations, volumes)
            continue
                
        if liquid == "annihilation":
            # Annihilation goes to all memories' condition 3 wells (one-to-many mapping)
            for annihilator in annihilators:
                destinations, volumes = find_destinations(liquid, annihilator, one_to_many=True)
                source = annihilator["source"]
                protocol.comment(f"Plating {liquid} from srouce: {source} to {destinations}.")
                plate(annihilator, destinations, volumes)
            continue

        if liquid == "reporter":
            # Reporters are plated into condition 1, 2 (one-to-one) and condition 3 (one-to-many)
            for memory in components.keys():
                # Plate reporter for conditions 1 and 2 (one-to-one)
                destinations, volumes = find_destinations(liquid, components, one_to_one=memory)
                source = components[memory][liquid]["source"]
                protocol.comment(f"Plating {liquid} for memory {memory} from srouce: {source} to {destinations}.")
                plate(components[memory][liquid], destinations, volumes)

            # Plate reporters for condition 3 (one-to-many across all memories)
            for memory in components.keys():
                destinations, volumes = find_destinations(liquid, components, one_to_many=True)
                source = components[memory][liquid]["source"]
                protocol.comment(f"Plating {liquid} for memory {memory} from srouce: {source} to {destinations}.")
                plate(components[memory][liquid], destinations, volumes)
            continue

        # Plate all other liquid types (one-to-many mapping across all conditions)
        for memory in components.keys():
            destinations, volumes = find_destinations(liquid, components, one_to_many=True)
            source = components[memory][liquid]["source"]
            protocol.comment(f"Plating {liquid} for memory {memory} from srouce: {source} to {destinations}.")
            plate(components[memory][liquid], destinations, volumes)

    return

def add_parameters(parameters: protocol_api.Parameters):
    """
    Define configurable parameters for the protocol that allow the user to control various aspects of the test run.
    These include the test type (assay), number of memories (1-4), number of replicates (1 or 3), and whether to run in simulation mode.

    Args:
        parameters (protocol_api.Parameters): The Parameters object provided by the Opentrons API for defining the protocol parameters.

    Parameters:
        - `Test`: Integer representing which test to run during the protocol.
            - Default: 0 (NONE)
            - Choices:
                - 0: NONE
                - 5: WTA (input/trigger: X)

        - `number_of_memories`: Integer representing the number of memories to plate in the protocol.
            - Default: 0 (NONE)
            - Choices: 0 (NONE), 1, 2, 3, or 4 memories.

        - `number_of_replicates`: Integer representing the number of replicates to run for each condition (single or triplicate).
            - Default: 0 (NONE)
            - Choices: 0 (NONE), 1 (Single), 3 (Triplicate)

        - `simulation`: Boolean to indicate whether to run the protocol in simulation mode.
            - Default: False (run on OT-2 robot)
    """
    parameters.add_int(
        variable_name = "Test",
        display_name = "Test",
        description = "Which test to run",
        default = 0,
        choices = [
            {"display_name": "NONE", "value": 0}, # default
            {"display_name": "WTA", "value": 5} #input/trigger: X. multiplication -> summation -> annihilation -> restoration -> reporter
        ]
    )
        
    parameters.add_int(
        variable_name = "number_of_memories",
        display_name = "Number of memories",
        description = "Number of memories/fluorophores to plate for",
        default = 0,
        choices = [
            {"display_name": "NONE", "value": 0},
            {"display_name": "1", "value": 1},
            {"display_name": "2", "value": 2},
            {"display_name": "3", "value": 3},
            {"display_name": "4", "value": 4}
        ]
    )
    
    parameters.add_int(
        variable_name = "number_of_replicates",
        display_name = "Number of replicates",
        description = "Number of replicates for each condition (single or triplicate)",
        default = 0,
        choices = [
            {"display_name": "NONE", "value": 0},
            {"display_name": "Single", "value": 1},
            {"display_name": "Triplicate", "value": 3}
        ]
    )
    
    parameters.add_bool(
        variable_name = "simulation",
        display_name = "simulation",
        description = "Turn on if running locally and not on robot",
        default = False
    )
    
# run protocol
def run(protocol: protocol_api.ProtocolContext):
    """
    Main protocol function for validating the fully assembled DNA Winner-Take-All Neural Network given different inputs. 
    It loads the required labware, sets up the liquid handling instruments, and performs the steps to plate buffer, 
    input, annihilators, and memory components for the SDR test.

    Args:
        protocol (ProtocolContext): The Opentrons protocol context, which provides access to the robot, pipettes, and labware.

    Steps:
        1. Load labware: Sets up racks (1.5 mL and 15 mL tubes, 300uL tips), pipettes (1x 300 single channel gen 2, 1x 1000 single channel gen 2), and plate (48 well, 1.6 mL flat) for the protocol.
        2. Define conditions and liquid types: Specifies which liquids are required for each condition.
        3. Load data: Loads the liquid components, inputs, and annihilators from the provided CSV.
        4. Map plate layout: Assigns wells for each memory and condition.
        5. Plate buffer: Distributes buffer to the appropriate wells.
        6. Plate liquids: Distributes the liquid components into the appropriate wells based on memory and liquid type.

    Notes:
        - The protocol supports both local simulation and execution on the OT-2 robot.
        - The number of memories and replicates are configurable through the protocol parameters.
    """    
    protocol.comment("Starting WTA assembly validation protocol.")
    # Step 1: Retrieve the test type parameter
    test = protocol.params.Test
    number_of_memories = protocol.params.number_of_memories
    if test == 0:
        protocol.comment("No test selected. Exiting the protocol.")
        return
    
    # Step 2: Load labware onto the robot
    liquid_eppendorf_rack_1 = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '1')
    if number_of_memories > 2: liquid_eppendorf_rack_2 = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '4')
    buffer_tube_rack = protocol.load_labware('opentrons_6_tuberack_falcon_50ml_conical', '5')
    tiprack_1000 = protocol.load_labware('opentrons_96_tiprack_1000ul', '7')
    tiprack_300 = protocol.load_labware('opentrons_96_tiprack_300ul', '8')
    p1000_single_left = protocol.load_instrument('p1000_single_gen2', 'left', tip_racks = [tiprack_1000])
    p300_single_right = protocol.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack_300])
    plate_48 = protocol.load_labware('corning_48_wellplate_1.6ml_flat', '3')
    
 
    # Step 3: Define conditions and liquid types 
    conditions = {
                1: ['reporter'], # each memory has a unique condition 1 with only its reporter
                2: ['reporter', 'output'], # each memory has a unique condition 2 with its reporter and output
                3: ['reporter', 'restoration', 'summation', 'annihilation', 'weight', 'input_fuel', 'output_fuel', 'input'] # all components in condition 3 would be of all memories being test except input (each memory's condition 3 has its unique input)
    }
    liquid_types = ["reporter", "restoration", "summation", 'annihilation', "weight", "input_fuel", "output_fuel", 'output', 'input']     
    #  max volume a pipette will handle: when plating any for condition 3 and doing 3 reps --> 1 any x 4 memories x 3 reps = 12 * 10 = 120 uL
            
    # Step 4: Load the liquid components depending on simulation mode or actual execution
    protocol.comment("Loading data.")
    if protocol.params.simulation: 
        # If in simulation mode, load data from local CSV files
        bundled_data, file_local_name = load_bundled_data()
        if number_of_memories>2: inputs, components, annihilators, buffer = load_computer(number_of_memories, conditions, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack, liquid_eppendorf_rack_2=liquid_eppendorf_rack_2, file_local_name=file_local_name)
        else: inputs, components, annihilators, buffer = load_computer(test, number_of_memories, conditions, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack, file_local_name=file_local_name)
    else: 
        # Load bundled data directly when running on the robot
        bundled_data = protocol.bundled_data
        if number_of_memories>2: inputs, components, annihilators, buffer = load_computer(number_of_memories, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack, buffer_tube_rack, liquid_eppendorf_rack_2=liquid_eppendorf_rack_2)
        else: inputs, components, annihilators, buffer = load_computer(test, number_of_memories, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack)
    
    # Step 5: Map plate layout based on test conditions, memories, and replicates
    num_replicates = protocol.params.number_of_replicates
    protocol.comment("Mapping plate layout.")
    plate_map = map_plate(conditions, plate_48, num_replicates, inputs, buffer)
    
    # Step 6: Plate buffer into the designated wells
    plate_buffer(protocol, p1000_single_left, buffer, plate_map)
    
    # Step 7: Plate the liquids
    plate_liquids(protocol, liquid_types, p300_single_right, inputs, annihilators, components, plate_map)

    return