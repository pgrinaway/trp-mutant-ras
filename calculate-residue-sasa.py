__author__ = 'Patrick B. Grinaway'

import mdtraj as md
import pdbfixer
import operator
from simtk.openmm.app import PDBFile



if __name__=="__main__":
    #arguments = docopt.docopt(__doc__, version='SASA calculator using MDTraj')
    input_filename = '4L8G/4L8G.pdb'
    output_filename = '4L8G/4L8G_trp_mutants.txt'
    fixer = pdbfixer.PDBFixer(filename=input_filename)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(True)
    PDBFile.writeFile(fixer.topology, fixer.positions, open('pdbfixed_4L8G.pdb','w'))
    traj = md.load('pdbfixed_4L8G.pdb')
    sasa = md.shrake_rupley(traj)[0]
    residue_sasa = {}
    for res in traj.top.residues:
        atm_indices = [atm.index for atm in res.atoms]
        residue_sasa[res.name+str(res.resSeq)] =  sasa[atm_indices].sum()
    sorted_sasa = sorted(residue_sasa.items(), key=operator.itemgetter(1))
    top_96 = sorted_sasa[len(sorted_sasa)-96:]
    mutation_list = "\n".join([res[0]+'TRP' for res in top_96])
    outfile = open(output_filename,'w')
    outfile.write(mutation_list)
    outfile.close()






