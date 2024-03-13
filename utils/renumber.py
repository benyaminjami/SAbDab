#!/usr/bin/python
# -*- coding:utf-8 -*-
from anarci import run_anarci
from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1
import numpy as np

def renumber_seq(seq, scheme="imgt"):
    _, numbering, details, _ = run_anarci(
        [("A", seq)], scheme=scheme, allowed_species=["mouse", "human"], hmmerpath='/h/benjami/hmmer/bin'
    )
    numbering = numbering[0]
    fv, position = [], []
    if not numbering:  # not antibody
        return None
    chain_type = details[0][0]["chain_type"]
    # numbering = numbering[0][0]
    mask = np.array([False for _ in range(len(seq))])
    for i, number in enumerate(numbering):
        for pos, res in number[0]:
            if res == "-":
                continue
            fv.append(res)
            position.append((pos[0] + (i * 129), pos[1]))
        temp_mask = np.logical_and(number[2] >= np.arange(len(seq)), number[1] <= np.arange(len(seq)))
        mask = np.logical_or(mask, temp_mask)
    return "".join(fv), position, chain_type, mask


def renumber_pdb(pdb, out_pdb, scheme="imgt", mute=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("anonym", pdb)
    for chain in structure.get_chains():
        seq = []
        for residue in chain:
            hetero_flag, _, _ = residue.get_id()
            if hetero_flag != " ":
                continue
            seq.append(seq1(residue.get_resname()))
        seq = "".join(seq)
        res = renumber_seq(seq, scheme)
        if res is None:
            continue
        fv, position, chain_type, mask = res
        if not mute:
            print(f"chain {chain.id} type: {chain_type}")

        seq_index, pos_index = -1, 0
        for r in list(chain.get_residues()):
            hetero_flag, _, _ = r.get_id()
            if hetero_flag != " ":
                continue
            seq_index += 1
            if not mask[seq_index]:
                chain.__delitem__(r.get_id())
                continue
            assert fv[pos_index] == seq1(r.get_resname()), f"Inconsistent residue in Fv {fv[pos_index]} at {r._id}"
            r._id = (" ", *position[pos_index])
            pos_index += 1
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb)


if __name__ == "__main__":
    import sys

    renumber_pdb('/h/benjami/scrach/dyMEAN/all_structures/imgt/1a0q.pdb', 'test.pdb', 'imgt', 'True')
    infile, outfile, scheme = sys.argv[1:4]
    if len(sys.argv) > 4:
        mute = bool(sys.argv[4])
    else:
        mute = False
    renumber_pdb(infile, outfile, scheme, mute)
