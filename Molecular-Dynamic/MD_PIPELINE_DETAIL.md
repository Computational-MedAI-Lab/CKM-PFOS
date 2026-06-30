# Molecular Dynamics Simulation Pipeline  
## PFOS–MMP9 Complex (GROMACS + AmberTools)

This module provides the **complete molecular dynamics (MD) simulation pipeline**  
used in the study:

The workflow includes:
- Ligand parameterization (AmberTools / ACPYPE)
- Protein preparation (PDB Fixer, GROMACS)
- Complex construction
- Energy minimization
- NVT & NPT equilibration
- Production MD simulation
- Post-simulation trajectory analysis

- ## Requirements

- **AmberTools ≥ 25**
- **GROMACS ≥ 2025**
- **Python ≥ 3.9** (for xvg analysis tools)
- **DuIvyTools v0.6.0**
Tested on Linux (Ubuntu 22.04).

## Quick Start

conda activate AmberTools25

antechamber -i PFOS.pdb -fi pdb -o PFOS.mol2 -fo mol2 -at gaff2 -c bcc -s 2

parmchk2 -i PFOS.mol2 -f mol2 -o PFOS.frcmod

acpype -i PFOS.mol2 -o gmx

pdbfixer MMP9.pdb --output=MMP9_fixed.pdb --add-atoms=all --add-residues --replace-nonstandard

gmx pdb2gmx -f MMP9_fixed.pdb -o MMP9.gro -p topol.top -ignh

gmx editconf -f Complex.gro -o Complex-box.gro -d 1.5 -bt cubic

gmx solvate -cp Complex-box.gro -o Complex-box-solvent.gro -p topol.top

gmx grompp -f em.mdp -c Complex-box-solvent.gro -p topol.top -o ions.tpr -maxwarn 1

gmx genion -s ions.tpr -p topol.top -o Complex-box-solvent-ion.gro -neutral

gmx grompp -f em.mdp -c Complex-box-solvent-ion.gro -r Complex-box-solvent-ion.gro -p topol.top -o em.tpr -maxwarn 1

gmx mdrun -v -deffnm em

gmx energy -f em.edr -o potential.xvg

dit xvg_show -f potential.xvg

gmx make_ndx -f em.gro -o em_index.ndx  
1 | 13  
1 & ! a H*  
13 & ! a H*  
q

gmx genrestr -f PFOS_GMX.gro -o posre_lig_1000.itp -fc 1000 1000 1000 
gmx genrestr -f PFOS_GMX.gro -o posre_lig_100.itp -fc 100 100 100 
gmx genrestr -f PFOS_GMX.gro -o posre_lig_10.itp -fc 10 10 10

modify posre_pro.itp to 1000/100/10 force

#join to topol.top the following parameter  
#ifdef POSRE  
#include "posre_lig_1000.itp"  
#endif

#ifdef POSRENPT  
#include "posre_lig_100.itp"  
#endif

#ifdef POSRENPT1  
#include "posre_lig_10.itp"  
#endif

#ifdef POSRE  
#include "posre_pro_1000.itp"      
#endif

#ifdef POSRENPT  
#include "posre_pro_100.itp"      
#endif

#ifdef POSRENPT1  
#include "posre_pro_10.itp"      
#endif

gmx grompp -f nvt.mdp -c em.gro -r em.gro -n em_index.ndx -p topol.top -o nvt.tpr

gmx mdrun -v -deffnm nvt

gmx energy -f nvt.edr -o temperture.xvg

gmx rms -s nvt.tpr -f nvt.xtc -o rmsd_nvt_pro.xvg -n em_index.ndx

gmx rms -s nvt.tpr -f nvt.xtc -o rmsd_nvt_lig.xvg

dit xvg_compare -f rmsd_nvt_pro.xvg rmsd_nvt_lig.xvg -c 1 1

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -n em_index.ndx -p topol.top -o npt.tpr

gmx mdrun -v -deffnm npt

gmx energy -f npt.edr -o pressure.xvg

gmx energy -f npt.edr -o density.xvg

gmx energy -f npt.edr -o volume.xvg

gmx rms -s npt.tpr -f npt.xtc -o rmsd_npt_pro.xvg -n em_index.ndx

gmx rms -s npt.tpr -f npt.xtc -o rmsd_npt_lig.xvg

gmx grompp -f npt1.mdp -c npt.gro -t npt.cpt -r npt.gro -n em_index.ndx -p topol.top -o npt1.tpr
gmx mdrun -v -deffnm npt1

gmx energy -f npt1.edr -o pressure1.xvg

gmx energy -f npt1.edr -o density1.xvg

gmx energy -f npt1.edr -o volume1.xvg

dit xvg_show -f pressure1.xvg density1.xvg volume1.xvg

gmx rms -s npt.tpr -f npt.xtc -o rmsd_npt1_pro.xvg -n em_index.ndx

gmx rms -s npt.tpr -f npt.xtc -o rmsd_npt1_lig.xvg

dit xvg_compare -f rmsd_npt1_pro.xvg rmsd_npt1_lig.xvg -c 1 1


gmx grompp -f md.mdp -c npt1.gro -t npt1.cpt -r npt1.gro -n em_index.ndx -p topol.top -o md.tpr

gmx mdrun -v -deffnm md

gmx trjconv -s md.tpr -f md.xtc -o md_nojump.xtc -pbc nojump

gmx trjconv -s md.tpr -f md_nojump.xtc -o md_nojump_center.xtc -center -pbc mol -n em_index.ndx

gmx trjconv -s md.tpr -f md_nojump_center.xtc -o md_nojump_center_fit.xtc -fit rot+trans -n em_index.ndx

gmx rms -s md.tpr -f md_nojump_center.xtc -o rmsd_pro.xvg -tu ns 

gmx rms -s md.tpr -f md_nojump_center.xtc -o rmsd_lig.xvg -tu ns 

gmx rms -s md.tpr -f md_nojump_center.xtc -n em_index.ndx -o rmsd_com.xvg -tu ns

gmx rmsf -s md.tpr -f md_nojump_center.xtc -o rmsf.xvg -res

gmx gyrate -s md.tpr -f md_nojump_center_fit.xtc -o gyrate.xvg -tu ns

gmx hbond-legacy -s md.tpr -f md_nojump_center.xtc -num hbond.xvg -tu ns -n em_index.ndx

gmx sasa -f md_nojump_center.xtc -s md.tpr -o sasa.xvg -n em_index.ndx

dit xvg_compare -f rmsd_pro.xvg rmsd_lig.xvg rmsd_com.xvg -c 1 1 1 -t RMSD -l "RMSD_Protein" "RMSD_Ligand" "RMSD_Complex" --legend_location outside

dit xvg_compare -f gyrate.xvg -t Rg -c 1,2,3,4 -l "Rg" "Rgx" "Rgy" "Rgz" --legend_location outside

dit xvg_show -f rmsf.xvg -t RMSF -t RMSF --legend_location outside

dit xvg_compare -f hbond.xvg --legend_location outside -c 1

dit xvg_show -f sasa.xvg
