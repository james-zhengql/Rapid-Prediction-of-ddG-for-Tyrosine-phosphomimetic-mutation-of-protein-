import os
import shutil
import subprocess

def compute_evoef_stability(pdb_path, mutation, evoef_path="../EvoEF_m/EvoEF", output_directory="./EvoEF-Repair"):
    """
    Computes EvoEF ΔΔG stability for a given PDB file and mutation.

    Args:
        pdb_path (str): Path to the PDB file.
        mutation (str): Mutation string (e.g., 'A50G' for Ala-50 to Gly mutation).
        evoef_path (str): Path to EvoEF executable.
        output_directory (str): Directory to store output files.

    Returns:
        float: Computed ΔΔG stability score, or None if computation fails.
    """
    if not os.path.isfile(pdb_path):
        print(f"File not found: {pdb_path}")
        return None

    os.makedirs(output_directory, exist_ok=True)

    # Generate paths for repaired structure
    base_filename = os.path.basename(pdb_path).replace('.pdb', '')
    repaired_pdb = os.path.join(output_directory, f"{base_filename}_Repair.pdb")

    # Step 1: Repair the PDB structure
    if not os.path.isfile(repaired_pdb):
        repair_command = [evoef_path, "--command=RepairStructure", f"--pdb={pdb_path}"]
        print(f"Repairing {pdb_path}...")
        process = subprocess.run(repair_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if process.returncode != 0:
            print(f"Error in repairing PDB: {process.stderr}")
            return None
        shutil.move(f"{base_filename}_Repair.pdb", repaired_pdb)

    # Step 2: Prepare mutant specification file
    mutant_file = os.path.join(output_directory, "individual.txt")
    with open(mutant_file, "w") as f:
        f.write(f"{mutation};\n")

    # Step 3: Generate mutant structure
    mutant_pdb = os.path.join(output_directory, f"{base_filename}_{mutation}.pdb")
    build_mutant_command = [evoef_path, "--command=BuildMutant", f"--pdb={repaired_pdb}", f"--mutant_file={mutant_file}"]
    print(f"Building mutant {mutation} for {repaired_pdb}...")
    process = subprocess.run(build_mutant_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if process.returncode != 0:
        print(f"Error in mutating PDB: {process.stderr}")
        return None

    shutil.move(f"{base_filename}_Repair_Model_0001.pdb", mutant_pdb)

    # Step 4: Compute stability for wild-type and mutant structures
    compute_stability_ref_command = [evoef_path, "--command=ComputeStability", f"--pdb={repaired_pdb}"]
    compute_stability_mut_command = [evoef_path, "--command=ComputeStability", f"--pdb={mutant_pdb}"]

    try:
        ref_output = subprocess.check_output(compute_stability_ref_command, text=True)
        mut_output = subprocess.check_output(compute_stability_mut_command, text=True)

        ΔGstability_ref = extract_stability_score(ref_output)
        ΔGstability_mut = extract_stability_score(mut_output)
        ΔΔGstability = ΔGstability_mut - ΔGstability_ref

        print(f"ΔΔGstability for {base_filename}, mutation {mutation}: {ΔΔGstability}")
        return ΔΔGstability

    except subprocess.CalledProcessError as e:
        print(f"Error in computing stability: {e}")
        return None

def extract_stability_score(output):
    """Extracts the total stability score from EvoEF output."""
    for line in output.splitlines():
        if "Total" in line:
            return float(line.split()[-1])
    return None
