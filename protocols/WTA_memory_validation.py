from opentrons import protocol_api
import csv
import io
import re

metadata = {
    "protocolName": "WTA Memory Validation",
    "description": "DNA WTA Individual Memories Validation Assays",
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
        tube_type (str): The type of tube, either 'buffer' (15 mL Falcon) or 'liquid' (1.5 mL Eppendorf).

    Returns:
        float: The calculated liquid height in millimeters (mm).

    Raises:
        ValueError: If the liquid volume exceeds the capacity of the tube.
    """
    # Define the cone and cylinder dimensions based on the tube type
    if tube_type == "buffer":
        cone_height = 2.0       # cm, height of the conical bottom of the Falcon tube
        cone_radius = 0.65      # cm, radius of the conical bottom
        cylinder_height = 10    # cm, height of the cylindrical section of the Falcon tube
        cylinder_radius = 0.75  # cm, radius of the cylindrical section
        offset = 10             # Correction offset for Falcon tubes
    elif tube_type == "liquid":
        cone_height = 0.5       # cm, height of the conical bottom of the Eppendorf tube
        cone_radius = 0.5       # cm, radius of the conical bottom
        cylinder_height = 1.5   # cm, height of the cylindrical section of the Eppendorf tube
        cylinder_radius = 0.6   # cm, radius of the cylindrical section
        offset = 2.5            # Correction offset for Eppendorf tubes
    else:
        raise ValueError("Unknown tube type. Valid types are 'buffer' or 'liquid'.")

    # Convert the liquid volume from µL to mL (1 cm^3 = 1 mL)
    liquid_volume = liquid_volume / 1000

    # Calculate the volume of the conical section (V_cone = 1/3 * pi * r^2 * h)
    V_cone = (1 / 3) * 3.14 * (cone_radius ** 2) * cone_height

    # If the liquid volume fits entirely within the conical section
    if liquid_volume <= V_cone:
        # Calculate the height of the liquid in the cone (h_cone = 3V / (pi*r^2))
        height_in_cone = (3 * liquid_volume) / (3.14 * (cone_radius ** 2))
        return (height_in_cone * 10) - offset  # Convert from cm to mm and apply offset

    # Otherwise, calculate the height in the cylindrical section
    else:
        # Volume left to fill in the cylindrical section
        volume_in_cylinder = liquid_volume - V_cone
        # Calculate the height in the cylindrical section (h_cylinder = V / (pi*r^2))
        height_in_cylinder = volume_in_cylinder / (3.14 * (cylinder_radius ** 2))

        # If the calculated height exceeds the cylinder height, raise an error
        if height_in_cylinder > cylinder_height:
            raise ValueError("Liquid volume exceeds tube capacity.")

        # Return the total height (cone height + cylindrical height), converted to mm and with the offset applied
        return ((cone_height + height_in_cylinder) * 10) - offset  # Convert from cm to mm and apply offset     
    
def load_computer(test, conditions, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack, liquid_eppendorf_rack_2=None, file_local_name=None):
    """
    Loads and parses the liquid data data from the CSV file for the protocol, organizing it into components, 
    buffers, and possibly annihilators or inputs, depending on the assay.

    Args:
        test (int): The test type / assay (1: Multiplication, 2: Summation, 3: Annihilation) that determines the organization of the liquid data.
        conditions (dict): A dictionary containing assay conditions {condition_number (int): [list of liquids in condition (str)]}.
        bundled_data (dict): A dictionary containing file data, {file_name (str): file_data (bytes or text)}.
        liquid_eppendorf_rack_1 (Labware): Labware object representing the first rack of liquid tubes.
        buffer_tube_rack (Labware): Labware object representing the rack holding the buffer tubes.
        liquid_eppendorf_rack_2 (Labware, optional): Labware object representing the second rack of liquid tubes (optional). Defaults to None.
        file_local_name (str, optional): Local name of the CSV file if bundled_data contains decoded text. Defaults to None.

    Returns:
        tuple: (components, buffer)
   
                - components (dict): {memory_number: {liquid_type: {source: (str), total_volume: (float), reaction_volume: (float)}}}
                - buffer (dict): {source: (str), total_volume: (float), memory_number: {condition_number: total buffer volume (float)}}

                    where:
                        from CSV file:
                            - memory_number (int or str): The memory number corresponding to the liquid components (0 for annihlation gates in test 3 and 'b' for buffer).
                            - liquid_type: The type of liquid component (e.g., input, buffer, etc.).
                            - source: The location of the corresponding tube in reference to the appropriate rack.
                            - total_volume: The total volume of the liquid component in the tube.
                        Standardized:
                            - total_reaction_volume: The total volume of the liquid component required for the reaction.
                            - reaction_volume: The volume of the liquid component required for the reaction.
    Raises:
        ValueError: If no valid SDR_computer CSV file is found in the provided data.
    """

    # Step 1: Find and read the SDR_computer CSV file
    pattern = re.compile(r'SDR_computer.*\.csv')  # Pattern to identify the correct file in bundled_data
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
    components = {}                           # Store liquid components of each memory (e.g., strands and gates)
    buffer = {}                               # Store buffer data
    
    total_reaction_volume = 100     # Total reaction volume per well
    standard_liquid_volume = 10.0   # Standard volume for liquid components

    # Step 3: Parse CSV data and organize it into components, inputs, buffers, etc.
    liquids_reader = csv.DictReader(io.StringIO(liquids_file))
    row_count = 0  # Keep track of the row number to handle conditions for test 1
    
    for row in liquids_reader:
        row_count += 1
        
        # Get memory number and determine its type (int or str)
        memory_number = int(row["memory_number"]) if isinstance(row["memory_number"], int) else str(row["memory_number"])
        
        # Assign the correct rack for the liquid based on memory number and test type
        if isinstance(memory_number, str) and memory_number == "b":
            liquid_rack = buffer_tube_rack         # Buffer memory gets the buffer rack
        elif test == 1 and row_count <= 19:
            liquid_rack = liquid_eppendorf_rack_1  # For test 1, the first 17 rows are in rack 1
        elif test == 1 and row_count > 19:
            liquid_rack = liquid_eppendorf_rack_2  # For test 1, beyond 17 rows, use rack 2 if available
        else:
            liquid_rack = liquid_eppendorf_rack_1  # Default to rack 1 for other cases

        # Extract relevant data from the row
        liquid_type = row['liquid_type']
        source = row['liquid_source']
        liquid_total_volume = float(row['liquid_total_volume'])

        # Handle 'buffer' liquids
        if liquid_type == "buffer":
            buffer["source"] = liquid_rack[source]
            buffer["total_volume"] = liquid_total_volume
            continue
        
        # For all other liquid types, store them as components
        if memory_number not in components:
            components[memory_number] = {}
        
        components[memory_number][liquid_type] = {
            "source": liquid_rack[source],
            "total_volume": liquid_total_volume,
            "reaction_volume": standard_liquid_volume
        }

        
    # Step 4: Calculate buffer volumes for each condition
    for memory in components.keys():
        for condition, liquids_in_condition in conditions.items():
            num_liquids = len(liquids_in_condition)  # Get the number of liquids in each condition
            buffer_well_volume = total_reaction_volume - (num_liquids * standard_liquid_volume)  # Calculate buffer volume
            
            if memory not in buffer:
                buffer[memory] = {}
            buffer[memory][condition] = buffer_well_volume

    # Step 5: Return the  dicts containing the parsed data based on assay
    return components, buffer


def map_plate(test, conditions, plate, num_replicates, memories, buffer):
    """
    Maps the destination wells on a 96-well plate for each memory and condition, taking into account the number 
    of replicates for each condition and the buffer volume required for each well.

    Args:
        test (int): The test type / assay, which affects the layout of wells on the plate.
        conditions (dict): A dictionary where keys are condition numbers and values are lists of liquid types.
        plate (Labware): The 96-well plate labware object used for the experiment.
        num_replicates (int): The number of replicates for each condition (1 or 3).
        memories (iterable): An iterable of memories being plated (components.keys()).
        buffer (dict): A dictionary containing buffer volume data for each well to be plated.

    Returns:
        dict: A dictionary (plate_map) that maps each memory and condition to specific well destinations, 
              the liquids for that condition, and the required buffer volume for each well.

    Notes:
        - The function loops over all memories and conditions, assigning wells based on memory and condition.
        - It adjusts the layout depending on the test type and number of replicates.
        - Conditions are plated along the columns and the replicates along the rows. Memories are plated clockewise in the order they appear in the CSV file.
    """
    plate_map = {}  # Initialize the mapping of wells to memories and conditions
    plate_rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']  # Define the rows of the 96-well plate (8 rows)
    plate_cols = list(range(1, 13))  # Define the columns of the 96-well plate (12 columns)

    row_index = 0  # Initialize row index
    col_index = 0  # Initialize column index

    # Step 1 Loop through each unique memory being plated
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

            # Step 5: Adjust layout based on test type
            # If test type is 2 (specific layout adjustment), reset column after 8 cols (2 memories * 4 conditions each) to plate the next memories (3 and/or 4) in the bottom half of the plate
            if test == 2 and col_index > 7:
                col_index = 0  # Reset column index
                row_index += num_replicates  # Move to the next row(s) (depending on replicates)

            # For all other test types, reset column after 12 cols (2 memories * 6 conditions each) to plate the next memories (3 and/or 4) in the bottom half of the plate
            if test != 2 and col_index > 11:
                col_index = 0  # Reset column index
                row_index += num_replicates  # Move to the next row(s) (depending on replicates)

    # Step 6: Return the completed plate map
    return plate_map

def plate_buffer(protocol, p300_single_left, p300_single_right, buffer, plate_map):
    """
    Distributes buffer into designated wells in the 96-well plate according to the buffer volumes 
    specified for each memory and condition. The function handles liquid transfers using both pipettes.

    Args:
        protocol (ProtocolContext): The Opentrons protocol context
        p300_single_left (InstrumentContext): The left-mount pipette (P300 single channel) for liquid transfers.
        p300_single_right (InstrumentContext): The right-mount pipette (P300 single channel) for liquid transfers.
        buffer (dict): A dictionary containing buffer volume data for each memory and condition.
        plate_map (dict): A dictionary mapping each memory and condition to specific well destinations, including the liquids in condition and buffer volumes.

    Returns:
        None: The function directly performs the buffer transfer steps using the protocol context.

    Notes:
        - Buffer is distributed to wells in descending order of buffer volume for optimal pipette usage.
        - The function aspirates and dispenses in chunks based on the maximum pipette volume.
        - Buffer levels in the source tube are updated as the liquid is transferred.
    """
    protocol.comment("Starting buffer plating.")
    
    air_gap = 10  # Air gap in µL to prevent dripping
    max_volume = 300 - air_gap  # Max pipette capacity, accounting for air gap

    def find_chunk_end_index(chunk_start_index, volumes, max_volume):
        """
        Helper function to determine how many wells can be filled in a single pipette load.
        It calculates the total volume that can be aspirated without exceeding the max volume.

        Args:
            chunk_start_index (int): The index to start checking volumes.
            volumes (list): List of buffer volumes to dispense into wells.
            max_volume (float): Maximum volume that can be aspirated by the pipette.

        Returns:
            int: The index where the accumulated volume exceeds the max_volume.
        """
        accumulated_volume = 0
        for i in range(chunk_start_index, len(volumes)):
            if accumulated_volume + volumes[i] <= max_volume:
                accumulated_volume += volumes[i]
            else:
                return i     # Return index where the volume exceeds the max_volume
        return len(volumes)  # Return the end if all volumes fit in a single chunk

    # Step 1: Define the source location for the buffer (tube in a rack)
    buffer_source = buffer["source"]

    # Step 2: Calculate the current liquid height in the buffer tube
    buffer_top = buffer_source.bottom(calculate_liquid_height(buffer["total_volume"], "buffer"))

    # Step 3: Prepare the destination wells and buffer volumes from the plate_map
    destinations = []
    volumes = []
    for memory, conditions in plate_map.items():
        for condition, volume in conditions.items():
            wells = plate_map[memory][condition]['destination']  # Get the destination wells for this condition
            buffer_volume = plate_map[memory][condition]['buffer_volume']  # Buffer volume required per well
            
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
    destinations, volumes = zip(*destinations_volumes)  # Unzip into separate lists

    # Step 5: Pick up tips for both pipettes (plate all wells with these tips as there is no need to change them since the buffer is the same)
    p300_single_left.pick_up_tip()
    p300_single_right.pick_up_tip()

    # Mix buffer in the source tube to ensure uniformity
    p300_single_left.mix(3, 250, buffer_top)

    # Step 6: Begin distributing the buffer across the plate in chunks
    left_chunk_start_index = 0
    left_chunk_end_index = 0
    while left_chunk_end_index < len(volumes):
        # Find the next chunk of wells to fill with the left pipette
        left_chunk_end_index = find_chunk_end_index(left_chunk_start_index, volumes, max_volume)
        left_chunk_volumes = volumes[left_chunk_start_index:left_chunk_end_index]
        left_chunk_volume = sum(left_chunk_volumes)
        left_chunk_destinations = destinations[left_chunk_start_index:left_chunk_end_index]

        # Calculate chunk for the right pipette, if there are more wells to fill
        right_chunk_start_index = left_chunk_end_index if left_chunk_end_index != len(volumes) else None
        if right_chunk_start_index:
            right_chunk_end_index = find_chunk_end_index(right_chunk_start_index, volumes, max_volume)
            right_chunk_volumes = volumes[right_chunk_start_index:right_chunk_end_index]
            right_chunk_volume = sum(right_chunk_volumes)
            right_chunk_destinations = destinations[right_chunk_start_index:right_chunk_end_index]

        # Aspirate buffer with the left pipette and dispense into wells
        p300_single_left.aspirate(left_chunk_volume, buffer_top)
        p300_single_left.touch_tip(buffer_source, radius=0.05, speed=80, v_offset=-1)
        p300_single_left.air_gap(air_gap, height=1)
        p300_single_left.move_to(buffer_source.top(10), speed=50)

        # Update buffer volume and recalculate the liquid height in the buffer tube
        new_buffer_volume = buffer["total_volume"] - left_chunk_volume
        buffer["total_volume"] = new_buffer_volume
        buffer_top = buffer_source.bottom(calculate_liquid_height(new_buffer_volume, "buffer"))

        # Aspirate and dispense buffer with the right pipette, if necessary
        if right_chunk_start_index:
            p300_single_right.aspirate(right_chunk_volume, buffer_top)
            p300_single_right.touch_tip(buffer_source, radius=0.05, speed=80, v_offset=-1)
            p300_single_right.air_gap(air_gap, height=1)
            p300_single_right.move_to(buffer_source.top(10), speed=50)

            # Update buffer volume and liquid height for the right pipette as well
            new_buffer_volume = buffer["total_volume"] - right_chunk_volume
            buffer["total_volume"] = new_buffer_volume
            buffer_top = buffer_source.bottom(calculate_liquid_height(new_buffer_volume, "buffer"))

        # Dispense buffer into wells for the left pipette
        for volume, well in zip(left_chunk_volumes, left_chunk_destinations):
            p300_single_left.move_to(well.top(1))
            p300_single_left.move_to(well.bottom(1))
            p300_single_left.dispense(air_gap + volume, well.bottom(1))
            p300_single_left.touch_tip(well, radius=0.05, speed=80, v_offset=-1)
            protocol.delay(1)  # Short delay after dispensing
            if well != left_chunk_destinations[-1]:  # Air gap before moving to the next well
                p300_single_left.air_gap(air_gap, height=1)
        p300_single_left.blow_out(left_chunk_destinations[-1].top(1))

        # Dispense buffer into wells for the right pipette, if used
        if right_chunk_start_index:
            for volume, well in zip(right_chunk_volumes, right_chunk_destinations):
                p300_single_right.move_to(well.top(1))
                p300_single_right.move_to(well.bottom(1))
                p300_single_right.dispense(air_gap + volume, well.bottom(1))
                p300_single_right.touch_tip(well, radius=0.05, speed=80, v_offset=-1)
                protocol.delay(1)  # Short delay after dispensing
                if well != right_chunk_destinations[-1]:  # Air gap before moving to the next well
                    p300_single_right.air_gap(air_gap, height=1)
            p300_single_right.blow_out(right_chunk_destinations[-1].top(1))

        # If both chunks are done or all wells are filled, exit the loop
        if (right_chunk_start_index and right_chunk_end_index == len(volumes)) or left_chunk_end_index == len(volumes):
            break

        # Move to the next chunk of wells
        left_chunk_start_index = right_chunk_end_index

    # Step 7: Drop tips when the process is complete
    p300_single_left.drop_tip()
    p300_single_right.drop_tip()

    return

def plate_liquids(protocol, liquid_types, p300_single_left, p300_single_right, components, plate_map):
    """
    Distributes the various liquid components to the designated wells in the plate, based on the assay and memory information.
    The function handles pairing liquids for simultaneous distribution using both pipettes.
    
    Plating logic: by a predefined order of liquid type pairs (liquid_types) to avoid early triggering of gates. 
                   The pairs are also in a predefined order to avoid contamination between conditions of a memory 
        
        High level:
                    loop over liquid type pairs (left and right) in liquid_types:
                        loop over each memory:
                            for memory x:
                                - find all conditions left goes into and all that right goes into 
                                - new tip --> plate left --> plate right in order of conditions (replicates before moving to next condition)
                                
                    The max possible volumes needed to distribute a memory's liquid type will never be more than 300 uL given 
                    the current setup of assays, so there is no need to change tips until moving to the next memory. Note that
                    this is also means that for any liquid to be plated, it will only need to be aspirated once.
                                
        
    Args:
        protocol (ProtocolContext): The Opentrons protocol context.
        liquid_types (list of list): List of liquid type pairs (left and right) specifying which liquids to distribute together and order of plating. 
        p300_single_left (InstrumentContext): The left-mount pipette (P300 single channel) used for the plating of the left liquid in the pair.
        p300_single_right (InstrumentContext): The right-mount pipette (P300 single channel) used for the plating of the right liquid in the pair.
        components (dict): A dictionary with the liquid component data for each memory, including source, total volume, and reaction volume.
        plate_map (dict): A dictionary mapping each memory and condition to specific well destinations, including liquids and buffer volumes.

    Returns:
        None: The function performs liquid transfers directly on the Opentrons protocol context.

    Notes:
        - The function distributes liquids by pairs (one liquid type per pipette).
        - Shakes tip, aspirates an air gap, and delays for 1 sec after each dispense to prevent dripping.
        - Predefined order of plating (liquid type list) is used to avoid early triggering of gates.
        - Assays with an odd number of liquid types will use the left pipette for plating the lone liquid type.
        - For Annihilation assay, all memories' components end up in the same wells, i.e. comparing conditions only not memories. 
            (treated as a single memory but need to ensure to change tips between sources)
    """
    protocol.comment("Starting components plating.")
    
    air_gap = 10  # Air gap in µL to prevent dripping during liquid transfers

    # Step 1: Loop through each pair of liquid types to be plated
    for left_right_liquid_types in liquid_types:
        left_liquid_type = left_right_liquid_types[0]  # First liquid type goes to the left pipette
        
        # If there's a second liquid type in the pair, assign it to the right pipette
        if len(left_right_liquid_types) == 2:
            right_liquid_type = left_right_liquid_types[1]

        # Step 2: Iterate over each memory in the components dictionary
        for memory, liquid_info in components.items():
            
            # Get the source and volume for the left liquid type
            left_liquid_source = liquid_info[left_liquid_type]['source']
            left_liquid_top = left_liquid_source.bottom(
                calculate_liquid_height(liquid_info[left_liquid_type]['total_volume'], "liquid")
            )
            left_liquid_reaction_volume = liquid_info[left_liquid_type]['reaction_volume']

            # If there is a right liquid type, get its source and volume
            if len(left_right_liquid_types) == 2:
                right_liquid_source = liquid_info[right_liquid_type]['source']
                right_liquid_top = right_liquid_source.bottom(
                    calculate_liquid_height(liquid_info[right_liquid_type]['total_volume'], "liquid")
                )
                right_liquid_reaction_volume = liquid_info[right_liquid_type]['reaction_volume']

            # Step 3: Collect destination wells and volumes for the left liquid type
            left_destinations = []
            left_volumes = []
            right_destinations = []
            right_volumes = []

            for conditions, mapping in plate_map[memory].items():
                wells = mapping['destination']
                
                # Check if left liquid type is required for the condition
                if left_liquid_type in mapping['liquids']:
                    if isinstance(wells, list):  # Multiple wells (replicates)
                        for well in wells:
                            left_destinations.append(well)
                            left_volumes.append(left_liquid_reaction_volume)
                    else:  # Single well
                        left_destinations.append(wells)
                        left_volumes.append(left_liquid_reaction_volume)
                
                # Check if right liquid type is required for the condition
                if len(left_right_liquid_types) == 2 and right_liquid_type in mapping['liquids']:
                    if isinstance(wells, list):  # Multiple wells (replicates)
                        for well in wells:
                            right_destinations.append(well)
                            right_volumes.append(right_liquid_reaction_volume)
                    else:  # Single well
                        right_destinations.append(wells)
                        right_volumes.append(right_liquid_reaction_volume)

            # Step 4: Pipetting process using the left pipette (and right if applicable)
 
            # Pick up tips for both pipettes
            p300_single_left.pick_up_tip()
            if len(left_right_liquid_types) == 2:
                p300_single_right.pick_up_tip()

            # Mix and aspirate liquid for the left pipette
            p300_single_left.mix(3, 50, left_liquid_top)
            p300_single_left.aspirate(sum(left_volumes), left_liquid_top)
            p300_single_left.touch_tip(left_liquid_source, radius=0.05, speed=80, v_offset=-1)
            p300_single_left.air_gap(air_gap, height=1)

            # Mix and aspirate liquid for the right pipette, if applicable
            if len(left_right_liquid_types) == 2:
                p300_single_right.mix(3, 50, right_liquid_top)
                p300_single_right.aspirate(sum(right_volumes), right_liquid_top)
                p300_single_right.touch_tip(right_liquid_source, radius=0.05, speed=80, v_offset=-1)
                p300_single_right.air_gap(air_gap, height=1)

            # Step 5: Dispense the liquids into the designated wells

            # Dispense using the left pipette
            for volume, well in zip(left_volumes, left_destinations):
                p300_single_left.move_to(well.top(1))
                p300_single_left.move_to(well.bottom(1), speed=50)
                p300_single_left.dispense(air_gap + volume, well.bottom(1))
                p300_single_left.touch_tip(well, radius=0.05, speed=80, v_offset=-1)
                protocol.delay(1)
                if well != left_destinations[-1]:
                    p300_single_left.air_gap(air_gap, height=1)
            p300_single_left.blow_out(left_destinations[-1].top(1))

            # Dispense using the right pipette, if applicable
            if len(left_right_liquid_types) == 2:
                for volume, well in zip(right_volumes, right_destinations):
                    p300_single_right.move_to(well.top(1))
                    p300_single_right.move_to(well.bottom(1), speed=50)
                    p300_single_right.dispense(air_gap + volume, well.bottom(1))
                    p300_single_right.touch_tip(well, radius=0.05, speed=80, v_offset=-1)
                    protocol.delay(1)
                    if well != right_destinations[-1]:
                        p300_single_right.air_gap(air_gap, height=1)
                p300_single_right.blow_out(right_destinations[-1].top(1))

            # Step 6: Drop the tips after the pipetting process is complete for current memory's current liquid types
            p300_single_left.drop_tip()
            if len(left_right_liquid_types) == 2:
                p300_single_right.drop_tip()

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
                - 1: Multiplication (input/trigger: X)
                - 2: Summation (input/trigger: P)
                - 3: Annihilation (input/trigger: S)
                - 4: Restoration (input/trigger: S)

        - `number_of_memories`: Integer representing the number of memories to plate in the protocol.
            - Default: 0 (NONE)
            - Choices: 0 (NONE), 1, 2, 3, or 4 memories.

        - `number_of_replicates`: Integer representing the number of replicates to run for each condition (single or triplicate).
            - Default: 0 (NONE)
            - Choices: 0 (NONE), 1 (Single), 3 (Triplicate)

        - `simulation`: Boolean to indicate whether to run the protocol in simulation mode.
            - Default: False (run on OT-2 robot)
    """

    # Add the 'Test' parameter to allow selection of different test types.
    parameters.add_int(
        variable_name="Test",
        display_name="Test",
        description="Which test to run",
        default=0,  # Default: NONE
        choices=[
            {"display_name": "NONE", "value": 0},            # No test selected
            {"display_name": "Multiplication", "value": 1},  # Multiplication test: input/trigger X
            {"display_name": "Summation", "value": 2},       # Summation test: input/trigger P
            {"display_name": "Annihilation", "value": 3},    # Annihilation test: input/trigger S
            {"display_name": "Restoration", "value": 4},     # Restoration test: input/trigger S
        ]
    )
    
    # Add the 'number_of_memories' parameter to specify how many memories to plate.
    parameters.add_int(
        variable_name="number_of_memories",
        display_name="Number of memories",
        description="Number of memories/fluorophores to plate for",
        default=0,  # Default: NONE
        choices=[
            {"display_name": "NONE", "value": 0},  # No memories selected
            {"display_name": "1", "value": 1},     # 1 memory
            {"display_name": "2", "value": 2},     # 2 memories
            {"display_name": "3", "value": 3},     # 3 memories
            {"display_name": "4", "value": 4},     # 4 memories
        ]
    )
    
    # Add the 'number_of_replicates' parameter to select how many replicates to perform per condition.
    parameters.add_int(
        variable_name="number_of_replicates",
        display_name="Number of replicates",
        description="Number of replicates for each condition (single or triplicate)",
        default=0,  # Default: NONE
        choices=[
            {"display_name": "NONE", "value": 0},        # No replicates selected
            {"display_name": "Single", "value": 1},      # Single replicate
            {"display_name": "Triplicate", "value": 3},  # Triplicate replicates
        ]
    )

    # Add the 'simulation' parameter to determine whether the protocol should run in simulation mode.
    parameters.add_bool(
        variable_name="simulation",
        display_name="Simulation mode",
        description="Turn on if running locally in simulation mode, and not on the OT-2 robot",
        default=False  # Default: run on the OT-2 robot
    )
    
def run(protocol: protocol_api.ProtocolContext):
    """
    Main protocol function for running various DNA Winner-Take-All Neural Network individual memories validation assays on the Opentrons OT-2 robot.
    This function handles the loading of labware, liquid handling instruments, and liquid components based on the 
    test type and other parameters. It performs buffer and liquid plating into a 96-well plate while managing 
    different test conditions, liquid components, and replicates.

    Args:
        protocol (ProtocolContext): The Opentrons protocol context that provides access to the robot, labware, pipettes, and parameters.

    Steps:
        1. Load Labware: Load the racks (1x or 2x 1.5 mL, and 1x 15 mL tubes, 1x 300uL tips), pipettes (2x 300 single channel gen 2), and plate (96 well, 200 uL flat) for the protocol.
        2. Define Conditions and Liquid Types: Set up conditions and liquid components based on the selected test.
        3. Load Liquid Components: Load data from CSV files or bundled data, depending on whether the simulation mode is active.
        4. Map Plate Layout: Map well positions on the 96-well plate based on memories, conditions, and replicates.
        5. Plate Buffer: Distribute buffer into the designated wells.
        6. Plate Liquids: Distribute liquid components into the appropriate wells based on the asssay and its conditions.
    
    Returns:
        None: The function directly controls the robot to execute the experiment.

    Assays Supported:
        - Multiplication: Test 1 (input/trigger: X)
        - Summation: Test 2 (input/trigger: P)
        - Annihilation: Test 3 (input/trigger: S)
        - Restoration: Test 4 (input/trigger: S)

    Raises:
        ValueError: If invalid parameters or configurations are encountered during execution.
    """

    # Step 1: Retrieve the test type parameter
    test = protocol.params.Test
    if test == 0:
        protocol.comment("No test selected. Exiting the protocol.")
        return

    # Step 2: Load labware onto the robot
    liquid_eppendorf_rack_1 = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '1')
    
    # For test 1 (Multiplication), load a second rack of Eppendorf tubes
    if test == 1:
        liquid_eppendorf_rack_2 = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '4')
    
    buffer_tube_rack = protocol.load_labware('opentrons_15_tuberack_falcon_15ml_conical', '5')
    tiprack_300 = protocol.load_labware('opentrons_96_tiprack_300ul', '7')
    plate_96 = protocol.load_labware('nest_96_wellplate_200ul_flat', '3')

    # Load pipettes
    p300_single_left = protocol.load_instrument('p300_single_gen2', 'left', tip_racks=[tiprack_300])
    p300_single_right = protocol.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack_300])

    # Step 3: Define conditions and liquid types for each test type
    if test == 4:  # Restoration test
        conditions = {
            1: ['reporter'],
            2: ['reporter', 'output'],
            3: ['reporter', 'restoration', 'sum'],
            4: ['reporter', 'restoration', 'sum_low'],
            5: ['reporter', 'restoration', 'output_fuel', 'sum'],
            6: ['reporter', 'restoration', 'output_fuel', 'sum_low']
        }
        liquid_types = [
            ["reporter", "restoration"],
            ["output_fuel", "output"],
            ["sum", "sum_low"]
        ]
        # max volume a pipette will handle: when plating reporter and doing 3 reps --> 6 conditions * 10 ul * 3 reps = 180 ul < 300
        # min volume a pipette will handle: when plating output and doing 1 rep --> 1 condition * 10 ul * 1 rep = 10 ul < 300
    
    if test == 3:  # Annihilation test (assumes only 2 memories: x and z)
        conditions = {
            1: ['reporter_x', 'output_x'],
            2: ['reporter_x', 'reporter_z', 'restoration_x', 'restoration_z', 'sum_x'],
            3: ['reporter_x', 'reporter_z', 'restoration_x', 'restoration_z', 'annihilation', 'sum_x', 'sum_lower_z'],
            4: ['reporter_x', 'reporter_z', 'restoration_x', 'restoration_z', 'annihilation', 'sum_x', 'sum_low_z'],
            5: ['reporter_x', 'reporter_z', 'restoration_x', 'restoration_z', 'annihilation', 'sum_z', 'sum_lower_x'],
            6: ['reporter_x', 'reporter_z', 'restoration_x', 'restoration_z', 'annihilation', 'sum_z', 'sum_low_x']
        }
        liquid_types = [
            ['reporter_x', 'reporter_z'],
            ['restoration_x', 'restoration_z'],
            ['annihilation','output_x'],
            ['sum_x', 'sum_lower_z'],
            ['sum_low_z', 'sum_z'],
            ['sum_lower_x','sum_low_x']
        ]
        # max volume a pipette will handle: when plating reporter x and doing 3 reps --> 6 conditions * 10 ul * 3 reps = 180 ul < 300
        # min volume a pipette will handle: when plating output x and doing 1 rep --> 1 condition * 10 ul * 1 rep = 10 ul < 300
    
    if test == 2:  # Summation test
        conditions = {
            1: ['reporter'],
            2: ['reporter', 'output'],
            3: ['reporter', 'restoration', 'output_fuel'],
            4: ['reporter', 'restoration', 'summation', 'output_fuel' , 'product']
        }
        liquid_types = [
            ["reporter", "restoration"],
            ["summation", "output_fuel"],
            ["product", "output"]
        ]
        # max volume a pipette will handle: when plating reporter and doing 3 reps --> 4 conditions * 10 ul * 3 reps = 120 ul < 300
        # min volume a pipette will handle: when plating output and doing 1 rep --> 1 condition * 10 ul * 1 rep = 10 ul < 300
    
    if test == 1:  # Multiplication test (assumes 1 input and 1 weight per memory)
        conditions = {
            1: ['reporter'],
            2: ['reporter', 'output'],
            3: ['reporter', 'restoration', 'summation', 'weight', 'output_fuel', 'input'],
            4: ['reporter', 'restoration', 'summation', 'weight', 'output_fuel', 'input_low'],
            5: ['reporter', 'restoration', 'summation', 'weight', 'output_fuel', 'input_fuel', 'input'],
            6: ['reporter', 'restoration', 'summation', 'weight', 'output_fuel', 'input_fuel', 'input_low']
        }
        liquid_types = [
            ["reporter", "restoration"],
            ["summation", 'weight'],
            ["output_fuel", 'input_fuel'],
            ["input", "input_low"],
            ['output']
        ]
        # max volume a pipette will handle: when plating reporter and doing 3 reps --> 6 conditions * 10 ul * 3 reps = 180 ul < 300
        # min volume a pipette will handle: when plating output and doing 1 rep --> 1 condition * 10 ul * 1 rep = 10 ul < 300

    # Step 4: Load the liquid components depending on simulation mode or actual execution
    if protocol.params.simulation:
        # If in simulation mode, load data from local CSV files
        bundled_data, file_local_name = load_bundled_data()
        if test == 1:
            components, buffer = load_computer(test, conditions, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack, liquid_eppendorf_rack_2=liquid_eppendorf_rack_2, file_local_name=file_local_name)
        else:
            components, buffer = load_computer(test, conditions, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack, file_local_name=file_local_name)
    else:
        # Load bundled data directly when running on the robot
        bundled_data = protocol.bundled_data
        if test == 1:
            components, buffer = load_computer(test, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack, liquid_eppendorf_rack_2=liquid_eppendorf_rack_2)

        else:
            components, buffer = load_computer(test, bundled_data, liquid_eppendorf_rack_1, buffer_tube_rack)

    # Step 5: Map plate layout based on test conditions, memories, and replicates
    num_replicates = protocol.params.number_of_replicates
    memories = components.keys()
    plate_map = map_plate(test, conditions, plate_96, num_replicates, memories, buffer)

    # Step 6: Plate buffer into the designated wells
    plate_buffer(protocol, p300_single_left, p300_single_right, buffer, plate_map)

    # Step 7: Plate the liquids into the designated wells
    plate_liquids(protocol, liquid_types, p300_single_left, p300_single_right, components, plate_map)

    return