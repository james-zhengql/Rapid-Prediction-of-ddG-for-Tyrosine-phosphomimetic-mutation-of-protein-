import os
import subprocess
import shutil

def compute_foldx_stability(pdb_path, mutation, foldx_path="./foldx_20251231", output_directory="./FoldX-Repair"):
    """
    Computes FoldX ΔΔG stability for a given PDB file and mutation.

    Args:
        pdb_path (str): Path to the PDB file.
        mutation (str): Mutation string (e.g., 'A50G' for Alanine-50 to Glycine mutation).
        foldx_path (str): Path to the FoldX executable.
        output_directory (str): Directory to store output files.

    Returns:
        float: Computed ΔΔG stability score, or None if computation fails.
    """
    if not os.path.isfile(pdb_path):
        print(f"File not found: {pdb_path}")
        return None

    os.makedirs(output_directory, exist_ok=True)

    # Generate paths for repaired structure
    pdb_filename = os.path.basename(pdb_path)
    pdb_dir = os.path.dirname(pdb_path)
    base_filename = pdb_filename.replace('.pdb', '')
    repaired_pdb = os.path.join(output_directory, f"{base_filename}_Repair.pdb")
    repaired_pdb_filename = f"{base_filename}_Repair.pdb"
    # Step 1: Repair the PDB structure
    if not os.path.isfile(repaired_pdb):
        repair_command = [
            foldx_path, "--command=RepairPDB", "--water=IGNORE",
            f"--output-dir={output_directory}", f"--pdb-dir={pdb_dir}", f"--pdb={pdb_filename}"
        ]
        # print(repair_command)
        print(f"Repairing {pdb_filename} using FoldX...")
        process = subprocess.run(repair_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if process.returncode != 0:
            print(f"Error repairing PDB: {process.stderr}")
            return None
        
        # # Move repaired file to output directory
        # shutil.move(f"{base_filename}_Repair.pdb", repaired_pdb)

    # Step 2: Create mutation file
    mutant_file = os.path.join(output_directory, "individual_list.txt")
    with open(mutant_file, "w") as f:
        f.write(f"{mutation};\n")

    # Step 3: Generate mutant structure
    build_mutant_command = [
        foldx_path, "--command=BuildModel", "--water=IGNORE",
        f"--mutant-file={mutant_file}", f"--output-dir={output_directory}",
        f"--pdb-dir={output_directory}",f"--pdb={repaired_pdb_filename}"
    ]
    print(f"Building mutant {mutation} using FoldX...")
    process = subprocess.run(build_mutant_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if process.returncode != 0:
        print(f"Error building mutant: {process.stderr}")
        return None

    # shutil.move(f"Optimized_{base_filename}_Repair_1.pdb", mutant_pdb)

    # Step 4: Extract ΔΔG stability score from FoldX output
    foldx_output_file = os.path.join(output_directory, f"Dif_{base_filename}_Repair.fxout")
    try:
        with open(foldx_output_file, "r") as file:
         for line in file:
            parts = line.strip().split('\t')
            if parts[0].endswith('.pdb'):
                entry =parts[1]
                return float(entry)
    except Exception as e:
        print(f"Error extracting FoldX stability score: {e}")
        return None
