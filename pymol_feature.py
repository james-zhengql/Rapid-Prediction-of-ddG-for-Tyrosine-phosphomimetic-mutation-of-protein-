from os.path import isfile
import pymol.cmd as cmd
import os

def extract_pymol_features(pdb_file, residue_num):
    """
    Extracts structural features from a PDB file using PyMOL.
    
    Args:
        pdb_file (str): Path to the PDB file.
        residue_num (int): Residue number for feature extraction.
    
    Returns:
        dict: A dictionary containing extracted features.
    """
    if not os.path.isfile(pdb_file):
        print(f"File not found: {pdb_file}")
        return None

    print(f"Processing {pdb_file} for residue {residue_num}")

    cmd.reinitialize()
    cmd.load(pdb_file, "prot")
    cmd.set('dot_solvent', 'on')

    # Calculate total number of residues in the PDB (assumes 1 alpha carbon per residue)
    total_residues = cmd.count_atoms('name CA')

    # Selecting the residue of interest
    try:
        cmd.select('sel', f'resi {residue_num} and chain A')
    except:
        cmd.select('sel', f'resi {residue_num}')

    # Area calculation
    area = cmd.get_area('sel')

    # Number of nearby atoms & residues
    cmd.select('sel2', f'all within 5 of sel')
    natom = cmd.count_atoms('sel2')
    cmd.select('sel2', 'byres sel2 and not sel')
    nres = cmd.count_atoms('sel2 and name CA')
    # ncooh = cmd.count_atoms('sel2 and name CA and (resname GLU or resname ASP)')

    # Secondary structure calculation
    # cmd.select('helix', 'sel and ss h')
    # cmd.select('sheet', 'sel and ss s')
    # ss = 0
    # if cmd.count_atoms('helix') > 0:
    #     ss = 1
    # elif cmd.count_atoms('sheet') > 0:
    #     ss = 2

    # Compute phi and psi angles
    phi, psi = None, None
    chains = cmd.get_chains("all")

    try:
        if len(chains) > 1:
            chain = 'A'  # Default to chain A, adjust if needed
            if cmd.count_atoms(f'resi {residue_num} and chain {chain}') > 0:
                phi = cmd.get_dihedral(f"name C and resi {residue_num-1} and chain {chain}",
                                       f"name N and resi {residue_num} and chain {chain}",
                                       f"name CA and resi {residue_num} and chain {chain}",
                                       f"name C and resi {residue_num} and chain {chain}")
                
                psi = cmd.get_dihedral(f"name N and resi {residue_num} and chain {chain}",
                                       f"name CA and resi {residue_num} and chain {chain}",
                                       f"name C and resi {residue_num} and chain {chain}",
                                       f"name N and resi {residue_num+1} and chain {chain}")
        else:
            if cmd.count_atoms(f'resi {residue_num}') > 0:
                phi = cmd.get_dihedral(f"name C and resi {residue_num-1}",
                                       f"name N and resi {residue_num}",
                                       f"name CA and resi {residue_num}",
                                       f"name C and resi {residue_num}")
                
                psi = cmd.get_dihedral(f"name N and resi {residue_num}",
                                       f"name CA and resi {residue_num}",
                                       f"name C and resi {residue_num}",
                                       f"name N and resi {residue_num+1}")
    except:
        print(f"Error calculating dihedrals for residue {residue_num}")

    # Calculate the mutation residue ratio
    # mut_res_length_ratio = residue_num / total_residues if total_residues else 0

    cmd.delete('all')  # Clean up PyMOL session

    return {
        "area": area,
        "num_residues": nres,
        # "num_atoms": natom,
        # "num_cooh_atoms": ncooh,
        # "residue_num": residue_num,
        # "secondary_structure": ss,
        "phi": phi,
        "psi": psi,
        # "total_residues": total_residues,
        # "mut_res_length_ratio": mut_res_length_ratio
    }

