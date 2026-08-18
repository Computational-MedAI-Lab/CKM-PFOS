import os
import subprocess

receptor_folder = "proteins_clean_pdbqt"  
ligand_folder = "ligand_pdbqt"      
output_folder = "output"                


os.makedirs(output_folder, exist_ok=True)


centers_file = "docking_centers.tsv"


with open(centers_file, "r") as f:
    lines = f.readlines()


for line in lines[1:]:
    cols = line.strip().split("\t")
    receptor_name, ligand_name, center_x, center_y, center_z = cols
    receptor_path = os.path.join(receptor_folder, receptor_name)
    ligand_path = os.path.join(ligand_folder, ligand_name)

    if not os.path.exists(receptor_path):
        print(f"Warning: no file in {receptor_path}")
        continue
    if not os.path.exists(ligand_path):
        print(f"Warning: no file in {ligand_path}")
        continue


    output_file = os.path.join(output_folder, f"{os.path.splitext(receptor_name)[0]}_{os.path.splitext(ligand_name)[0]}_out.pdbqt")
    log_file = os.path.join(output_folder, f"{os.path.splitext(receptor_name)[0]}_{os.path.splitext(ligand_name)[0]}_log.txt")


    docking_command = [
        "vina",
        "--receptor", receptor_path,
        "--ligand", ligand_path,
        "--center_x", center_x,
        "--center_y", center_y,
        "--center_z", center_z,
        "--size_x", "20",
        "--size_y", "20",
        "--size_z", "20",
        "--out", output_file,
        "--exhaustiveness", "16"
    ]


    with open(log_file, "w") as log:
        subprocess.run(docking_command, stdout=log, stderr=log)
    
    print(f"Finished docking: {receptor_name} + {ligand_name}")
