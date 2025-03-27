def map_plate(test, num_replicates, num_memories):

    conditions = {
        1: ['reporter'],
        2: ['reporter', 'output'],
        3: ['reporter', 'restoration', 'summation', 'weight', 'output_fuel', 'input'],
        4: ['reporter', 'restoration', 'summation', 'weight', 'output_fuel', 'input_low'],
        5: ['reporter', 'restoration', 'summation', 'weight', 'output_fuel', 'input_fuel', 'input'],
        6: ['reporter', 'restoration', 'summation', 'weight', 'output_fuel', 'input_fuel', 'input_low']
    }

    plate_map = {}  # Initialize the mapping of wells to memories and conditions
    plate_rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']  # Define the rows of the 96-well plate (8 rows)
    plate_cols = list(range(1, 13))  # Define the columns of the 96-well plate (12 columns)

    row_index = 0  # Initialize row index
    col_index = 0  # Initialize column index

    # Step 1 Loop through each unique memory being plated
    for memory in range(1, num_memories + 1):

        # Initialize the entry in the plate map for the current memory
        plate_map[memory] = {}

        # Step 2: Loop through each condition and its corresponding liquids
        for condition, liquids_in_condition in conditions.items():
            # Step 3: Assign wells based on the number of replicates (1 or 3)
            if num_replicates == 1:
                # Single replicate: Assign one well for the condition
                plate_map[memory][condition] = {
                    'destination': f"{plate_rows[row_index]}{plate_cols[col_index]}",
                    'liquids': liquids_in_condition,
                    'buffer_volume': 100
                }
            elif num_replicates == 3:
                # Triple replicates: Assign three consecutive wells in the same column for the condition
                plate_map[memory][condition] = {
                    'destination': [
                        f"{plate_rows[row_index]}{plate_cols[col_index]}",
                        f"{plate_rows[row_index + 1]}{plate_cols[col_index]}",
                        f"{plate_rows[row_index + 2]}{plate_cols[col_index]}" ],
                    'liquids': liquids_in_condition,
                    'buffer_volume': 100
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


plate_map = map_plate(test=1, num_replicates=3, num_memories=4)
print(plate_map)