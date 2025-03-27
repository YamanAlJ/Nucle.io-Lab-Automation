import csv 

def load_computer(liquids_file, buffer_file):
    '''
    Load assay strands/gates information from csv files 
    '''
    
    liquids = {}
    buffer = {}
    
    # Load liquids data
    with open(liquids_file, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            memory_number = int(row["memory_number"])
            liquid_type = row['liquid_type']
            source = row['liquid_source']
            liquid_total_volume = float(row['liquid_total_volume'])
            print(row)
            liquid_reaction_volume = float(row['liquid_reaction_volume'])
            
            if memory_number not in liquids: liquids[memory_number] = {}
            
            liquids[memory_number][liquid_type] = {
                "source": [source],
                "liquid_total_volume": liquid_total_volume,
                "liquid_reaction_volume": liquid_reaction_volume
            }
        
    # Load buffer data
    with open(buffer_file, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            memory_number = int(row["memory_number"])
            condition_number = int(row['condition_number'])
            buffer_well_volume = float(row['buffer_well_volume'])
            
            if memory_number not in buffer: buffer[memory_number] = {}
            
            buffer[memory_number][condition_number] = buffer_well_volume
            
    return liquids, buffer


csv_liquids_file_path = 'protocols\sdr4_liquids.csv'
csv_buffer_file_path = 'protocols\sdr4_buffer.csv'
liquids, buffer = load_computer(csv_liquids_file_path, csv_buffer_file_path)
print(liquids)
print(buffer)