*Pleaese refer to the PDF guide First*

# DNA WTA Validation Assay Guide
This guide provides step-by-step instructions for running two key protocols involved in the validation of DNA WTA (Weighted Transform Arithmetic) computers using Opentrons OT-2. The two protocols described are for:

    1. Individual Memory Validation Assays (SDR)
    2. Assembly Validation Assay (WTA)

The following sections outline the experimental setup, required labware, how to configure the protocols, and specific instructions for running these assays.

## Contents

    1. Contents
    2. Overview
    3. Labware and Setup
    4. Individual Memory Validation (SDR)
    5. Protocol Setup
    6. Running the SDR Protocol
    7. Assembly Validation Assay (WTA)
    8. Protocol Setup
    9. Running the WTA Protocol
    10. Assay Templates and Sample Files

## Overview

The validation of DNA WTA involves two core assays:

    1. SDR (Signal-Domain Recognition) for validating the behavior of individual memories.
    2. WTA (Weighted Transform Arithmetic) for validating the correct assembly and functioning of the complete DNA computer.

## Assay Types

    - Multiplication Assay: Validates multiplication logic across individual memories.
    - Summation Assay: Ensures correct summation of memory states.
    - Annihilation Assay: Tests annihilation behavior in the circuit.
    - Restoration Assay: Focuses on memory restoration.
    - WTA Validation Assay: Ensures correct functioning of the assembled DNA circuit.

## Labware and Setup

    Required Labware:
        Pipettes:
            P300 Single Gen2 (Left)
            P300 Single Gen2 (Right)
        Tip Racks:
            Opentrons 96 Tip Rack 300µL
        Tube Racks:
            Opentrons 24 Tube Rack Eppendorf 1.5mL
            Opentrons 15 Tube Rack Falcon 15mL Conical
            Opentrons 6 Tube Rack Falcon 50mL Conical (for WTA)
        Plates:
            Nest 96 Well Plate 200µL Flat (for SDR)
            Corning 48 Well Plate 1.6mL Flat (for WTA)

## Individual Memory Validation (SDR)

This protocol is designed to validate individual memories within the DNA WTA system. Each memory is associated with a different fluorophore and can store and compute specific values.

### WTA Protocol Setup

1. Test Types:

    - Multiplication Assay
    - Summation Assay
    - Annihilation Assay
    - Restoration Assay

2. Input Parameters:

    - Test: Select the desired test (Multiplication, Summation, etc.)
    - number_of_memories: Number of memories to validate (1-4).
    - number_of_replicates: Choose between single or triplicate replicates.
    - simulation: Enable simulation mode for local testing.

3. Labware Placement:

    - Position the 24-tube Eppendorf rack and 15mL Falcon rack in deck slots 1 and 5.
    - Load the tip racks and ensure the pipettes are loaded and calibrated.

### Running the SDR Protocol

1. Load CSV Data: The SDR_computer.csv file should contain all the liquid compositions, volumes, and source locations. Ensure the file follows the correct template format for the selected test.

2. Execution:

    - Initiate the protocol from the Opentrons App.
    - Select the appropriate test type.
    - Follow the on-screen instructions to load labware and confirm placements.
    - The protocol will automatically handle liquid transfers, buffer plating, and memory validation.

3. Data Collection: Once the protocol completes, analyze the fluorescent output data for the individual memories to verify the expected behavior. Do not waste time and be quick!

## Assembly Validation Assay (WTA)

This protocol validates the full assembly and functioning of the DNA WTA computer. The assay checks whether the complete DNA circuit, with multiple memories and connections, works as expected when given different inputs.

### Protocol Setup

    1. Test Types:
        - WTA Validation

    2. Input Parameters:
        - number_of_memories: Number of memories in the DNA WTA circuit.
        - number_of_replicates: Choose between single or triplicate replicates.
        - simulation: Enable simulation mode for local testing.

    3. Labware Placement:
        - Place the Eppendorf racks and Falcon racks according to the deck map.
        - Ensure pipettes and tip racks are loaded and calibrated.

### Running the WTA Protocol

    1. Load CSV Data: The WTA_computer.csv file should define the required liquid volumes, sources, and components for assembling the full DNA WTA system.

    2. Execution:
        - Load the WTA protocol from the Opentrons App.
        - Configure the number of memories and replicates.
        - Follow the app's instructions for placing labware and starting the run.
    
    3. Validation: After the protocol completes, review the output data and confirm that the assembled DNA WTA circuit functions correctly for the given inputs.

## Assay Templates and Sample Files

For running the assays, ensure you have the correct assay templates. You can find the necessary templates in the links below. Each template corresponds to a different assay type:

- Multiplication Assay Template
- Summation Assay Template
- Annihilation Assay Template
- Restoration Assay Template
- WTA Validation Template

### Template Format

Each template contains the following columns:

- memory_number: The memory being tested.
- liquid_type: Type of liquid (e.g., reporter, restoration, output).
- liquid_source: The source tube for each liquid.
- liquid_total_volume: The total volume of liquid to be used.
- reaction_volume: The volume of liquid to be added to each well.

Ensure your assay files are filled out according to the templates before starting the protocol.

## Notes and Best Practices

- Always confirm that your CSV files are correctly formatted and that the data matches the assay you are running.
- For simulation runs, enable the simulation parameter to avoid running the protocol on the actual OT-2 robot.
- If running multiple tests, ensure to label and organize your outputs carefully for easy comparison.

This guide should help you set up and run DNA WTA validation assays on the Opentrons OT-2 platform. Ensure to follow the specific instructions for each protocol to achieve accurate results.
