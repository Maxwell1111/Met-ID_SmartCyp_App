#!/usr/bin/env python3
"""
SMARTCyp-Inspired Metabolic Site Predictor

A Python/RDKit implementation inspired by SMARTCyp's approach to metabolism prediction.
Uses SMARTS patterns with quantitative scoring to rank metabolic sites.

Based on the SMARTCyp methodology:
    Rydberg P., Gloriam D.E., Zaretzki J., Breneman C., Olsen L.
    SMARTCyp: A 2D Method for Prediction of Cytochrome P450-Mediated Drug Metabolism.
    ACS Med Chem Lett. 2010;1(3):96-100.

Author: Chemistry Models Project
Date: 2025-12-15
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import sys


@dataclass
class MetabolicSite:
    """Data class for a predicted metabolic site"""
    atom_idx: int
    atom_symbol: str
    score: float
    rank: int
    reaction_type: str
    cyp_isoform: str
    energy: float  # Approximate activation energy (kcal/mol)
    accessibility: float  # Steric accessibility score (0-1)


class SMARTCypInspiredPredictor:
    """
    SMARTCyp-inspired metabolism predictor with quantitative scoring.

    Uses SMARTS patterns combined with:
    - Activation energy estimates
    - Steric accessibility
    - Electronic effects
    - CYP isoform specificity
    """

    # Reactivity rules with approximate activation energies (kcal/mol)
    # Lower energy = more reactive = higher metabolic liability
    REACTIVITY_RULES = [
        # N-Dealkylation (CYP3A4, 2D6)
        {
            'name': 'Aromatic N-Methyl',
            'smarts': '[n;R][CH3]',
            'base_energy': 75.0,  # Relatively activated
            'cyp_isoform': 'CYP3A4/2D6',
            'reaction': 'N-Demethylation'
        },
        {
            'name': 'Aliphatic N-Methyl',
            'smarts': '[NX3;!R][CH3]',
            'base_energy': 78.0,
            'cyp_isoform': 'CYP2D6/3A4',
            'reaction': 'N-Demethylation'
        },

        # Benzylic Oxidation (CYP2D6, 3A4)
        {
            'name': 'Benzylic Primary',
            'smarts': '[CH2;!R][c]',
            'base_energy': 73.0,  # Very reactive
            'cyp_isoform': 'CYP2D6/3A4',
            'reaction': 'Benzylic Hydroxylation'
        },
        {
            'name': 'Benzylic Secondary',
            'smarts': '[CH;!R]([*])[c]',
            'base_energy': 74.0,
            'cyp_isoform': 'CYP2D6/3A4',
            'reaction': 'Benzylic Hydroxylation'
        },

        # Allylic Oxidation
        {
            'name': 'Allylic CH2',
            'smarts': '[CH2;!R][C]=[C]',
            'base_energy': 76.0,
            'cyp_isoform': 'CYP3A4',
            'reaction': 'Allylic Hydroxylation'
        },

        # Aromatic Hydroxylation (CYP1A2, 2C9, 2D6, 3A4)
        {
            'name': 'Aromatic C-H (ortho to substituent)',
            'smarts': '[c;H1]',
            'base_energy': 85.0,  # Less reactive than benzylic
            'cyp_isoform': 'CYP1A2/2D6',
            'reaction': 'Aromatic Hydroxylation'
        },

        # Aliphatic Hydroxylation
        {
            'name': 'Cyclohexyl CH',
            'smarts': '[C;R1;H1]1[C;R1][C;R1][C;R1][C;R1][C;R1]1',
            'base_energy': 87.0,
            'cyp_isoform': 'CYP3A4',
            'reaction': 'Aliphatic Hydroxylation'
        },
        {
            'name': 'Cyclohexyl CH2',
            'smarts': '[C;R1;H2]1[C;R1][C;R1][C;R1][C;R1][C;R1]1',
            'base_energy': 86.0,
            'cyp_isoform': 'CYP3A4',
            'reaction': 'Aliphatic Hydroxylation'
        },

        # O-Dealkylation
        {
            'name': 'O-Methyl Ether',
            'smarts': '[OX2][CH3]',
            'base_energy': 77.0,
            'cyp_isoform': 'CYP2D6/3A4',
            'reaction': 'O-Demethylation'
        },

        # Terminal Methyl Oxidation
        {
            'name': 'Isopropyl Methyl',
            'smarts': '[CH3][CH;!R]([CH3])',
            'base_energy': 90.0,  # Less reactive
            'cyp_isoform': 'CYP3A4',
            'reaction': 'ω-Oxidation'
        },

        # Heterocycle Oxidation
        {
            'name': 'Pyridine C-H',
            'smarts': '[c;H1]1[c][n][c][c][c]1',
            'base_energy': 84.0,
            'cyp_isoform': 'CYP3A4',
            'reaction': 'Heteroaromatic Hydroxylation'
        },

        # ========== PHASE II METABOLISM ==========

        {
            'name': 'Phenolic OH',
            'smarts': '[OH]c1ccccc1',
            'base_energy': 65.0,  # Very reactive for conjugation
            'cyp_isoform': 'UGT/SULT',
            'reaction': 'Glucuronidation/Sulfation'
        },
        {
            'name': 'Primary Alcohol',
            'smarts': '[CH2][OH]',
            'base_energy': 70.0,  # High reactivity
            'cyp_isoform': 'UGT',
            'reaction': 'O-Glucuronidation'
        },
        {
            'name': 'Secondary Alcohol',
            'smarts': '[CH;!R][OH]',
            'base_energy': 72.0,
            'cyp_isoform': 'UGT',
            'reaction': 'O-Glucuronidation'
        },
        {
            'name': 'Carboxylic Acid',
            'smarts': '[CX3](=O)[OH]',
            'base_energy': 68.0,
            'cyp_isoform': 'UGT',
            'reaction': 'Acyl Glucuronidation'
        },
        {
            'name': 'Primary Amine',
            'smarts': '[NX3;H2;!R]',
            'base_energy': 74.0,
            'cyp_isoform': 'NAT1/NAT2',
            'reaction': 'N-Acetylation'
        },
        {
            'name': 'Aromatic Amine',
            'smarts': '[NH2]c1ccccc1',
            'base_energy': 67.0,  # Very reactive
            'cyp_isoform': 'NAT1/NAT2',
            'reaction': 'N-Acetylation'
        },
        {
            'name': 'Thiol',
            'smarts': '[SH]',
            'base_energy': 66.0,  # Very reactive
            'cyp_isoform': 'GST',
            'reaction': 'Glutathione Conjugation'
        },
        {
            'name': 'Catechol',
            'smarts': '[OH]c1ccc([OH])cc1',
            'base_energy': 64.0,  # Extremely reactive
            'cyp_isoform': 'COMT/UGT',
            'reaction': 'O-Methylation/Glucuronidation'
        },
    ]

    def __init__(self, smiles: str):
        """Initialize predictor with SMILES string"""
        self.smiles = smiles
        self.mol = self._parse_smiles(smiles)
        self.sites = []

    def _parse_smiles(self, smiles: str) -> Chem.Mol:
        """Parse and prepare molecule"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        # Add hydrogens for accurate pattern matching
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)
        return mol

    def _calculate_steric_accessibility(self, atom_idx: int) -> float:
        """
        Calculate steric accessibility of an atom (0-1 scale)
        Higher = more accessible = more likely to be metabolized
        """
        atom = self.mol.GetAtomWithIdx(atom_idx)

        # Count number of neighboring heavy atoms
        neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
        num_neighbors = len(neighbors)

        # Base accessibility on degree
        if num_neighbors == 0:
            accessibility = 1.0
        elif num_neighbors == 1:
            accessibility = 0.9
        elif num_neighbors == 2:
            accessibility = 0.7
        elif num_neighbors == 3:
            accessibility = 0.5
        else:
            accessibility = 0.3

        # Bonus for being in a ring (often less accessible)
        if atom.IsInRing():
            accessibility *= 0.8

        # Bonus for terminal positions
        if num_neighbors == 1:
            accessibility *= 1.1

        return min(1.0, accessibility)

    def _calculate_score(self, base_energy: float, accessibility: float) -> float:
        """
        Calculate metabolic liability score (0-100)
        Higher score = more likely to be metabolized

        Formula inspired by SMARTCyp: Score = f(Energy, Accessibility)
        """
        # Normalize energy to 0-100 scale (lower energy = higher score)
        # Typical range: 70-95 kcal/mol
        energy_score = max(0, 100 - (base_energy - 70) * 2)

        # Combine with accessibility
        final_score = energy_score * accessibility

        return round(final_score, 1)

    def predict(self) -> List[MetabolicSite]:
        """
        Predict metabolic sites with quantitative scoring

        Returns ranked list of metabolic sites
        """
        site_dict = {}  # atom_idx -> best site

        for rule in self.REACTIVITY_RULES:
            pattern = Chem.MolFromSmarts(rule['smarts'])
            if pattern is None:
                continue

            matches = self.mol.GetSubstructMatches(pattern)

            for match in matches:
                # Get the reactive atom (usually first in pattern)
                reactive_atom_idx = match[0]
                atom = self.mol.GetAtomWithIdx(reactive_atom_idx)

                # Calculate accessibility
                accessibility = self._calculate_steric_accessibility(reactive_atom_idx)

                # Calculate score
                score = self._calculate_score(rule['base_energy'], accessibility)

                # Keep best score for each atom
                if reactive_atom_idx not in site_dict or score > site_dict[reactive_atom_idx].score:
                    site_dict[reactive_atom_idx] = MetabolicSite(
                        atom_idx=reactive_atom_idx,
                        atom_symbol=atom.GetSymbol(),
                        score=score,
                        rank=0,  # Will assign after sorting
                        reaction_type=rule['reaction'],
                        cyp_isoform=rule['cyp_isoform'],
                        energy=rule['base_energy'],
                        accessibility=accessibility
                    )

        # Sort by score (descending) and assign ranks
        self.sites = sorted(site_dict.values(), key=lambda x: x.score, reverse=True)
        for i, site in enumerate(self.sites, 1):
            site.rank = i

        return self.sites

    def get_top_sites(self, n: int = 5) -> List[MetabolicSite]:
        """Get top N predicted sites"""
        if not self.sites:
            self.predict()
        return self.sites[:n]

    def generate_report(self) -> str:
        """Generate detailed text report"""
        if not self.sites:
            self.predict()

        report = []
        report.append("=" * 90)
        report.append("SMARTCYP-INSPIRED METABOLIC SITE PREDICTION")
        report.append("=" * 90)
        report.append(f"\nSMILES: {self.smiles}")

        # Remove hydrogens for display
        mol_no_h = Chem.RemoveHs(self.mol)
        report.append(f"Molecular Formula: {rdMolDescriptors.CalcMolFormula(mol_no_h)}")
        report.append(f"Molecular Weight: {rdMolDescriptors.CalcExactMolWt(mol_no_h):.2f} g/mol")
        report.append(f"Heavy Atoms: {mol_no_h.GetNumHeavyAtoms()}")

        report.append("\n" + "=" * 90)
        report.append("TOP METABOLIC SITES (Ranked by Score)")
        report.append("=" * 90)
        report.append(f"\n{'Rank':<6} {'Atom':<10} {'Score':<8} {'Energy':<10} {'Access':<8} {'Reaction':<25} {'CYP'}")
        report.append("-" * 90)

        for site in self.sites[:10]:  # Show top 10
            report.append(
                f"{site.rank:<6} "
                f"#{site.atom_idx:<4} ({site.atom_symbol:<2}) "
                f"{site.score:<8.1f} "
                f"{site.energy:<10.1f} "
                f"{site.accessibility:<8.2f} "
                f"{site.reaction_type:<25} "
                f"{site.cyp_isoform}"
            )

        report.append("\n" + "=" * 90)
        report.append("INTERPRETATION GUIDE:")
        report.append("  Score:   0-100 (higher = more likely to be metabolized)")
        report.append("  Energy:  Activation energy estimate in kcal/mol (lower = more reactive)")
        report.append("  Access:  Steric accessibility 0-1 (higher = more accessible)")
        report.append("=" * 90)

        return "\n".join(report)

    def visualize(self, output_path: str = "smartcyp_prediction.png",
                  width: int = 1600, height: int = 1200, top_n: int = 5):
        """Generate annotated structure with top sites highlighted (without hydrogens)"""
        if not self.sites:
            self.predict()

        # Get top N sites
        top_sites = self.sites[:top_n]

        # Remove hydrogens for clean visualization (like pattern matching tool)
        mol_no_h = Chem.RemoveHs(self.mol)
        AllChem.Compute2DCoords(mol_no_h)

        # Map atom indices from mol_with_h to mol_no_h
        atom_map = {}  # maps H-containing mol idx -> no-H mol idx
        h_idx = 0
        no_h_idx = 0
        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() != 1:  # Not hydrogen
                atom_map[atom.GetIdx()] = no_h_idx
                no_h_idx += 1

        # Map highlighted atoms to no-H indices
        highlight_atoms = []
        atom_colors = {}
        for site in top_sites:
            if site.atom_idx in atom_map:
                mapped_idx = atom_map[site.atom_idx]
                highlight_atoms.append(mapped_idx)

                # Color gradient: red (highest) to orange (lower)
                normalized = site.score / 100.0
                atom_colors[mapped_idx] = (1.0, max(0.2, 1.0 - normalized * 0.8), 0.0)

        # Prepare bond highlighting
        bond_colors = {}
        for bond in mol_no_h.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

            # If both atoms are highlighted, color the bond
            if begin_idx in highlight_atoms and end_idx in highlight_atoms:
                # Use average color of the two atoms
                color1 = atom_colors.get(begin_idx, (1.0, 0.5, 0.0))
                color2 = atom_colors.get(end_idx, (1.0, 0.5, 0.0))
                avg_color = tuple((c1 + c2) / 2 for c1, c2 in zip(color1, color2))
                bond_colors[bond.GetIdx()] = avg_color

        # Create drawer
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().clearBackground = True
        drawer.drawOptions().bondLineWidth = 3
        drawer.drawOptions().atomLabelFontSize = 30
        drawer.drawOptions().legendFontSize = 20
        drawer.drawOptions().padding = 0.05

        # Build legend
        legend_lines = [f"METABOLIC SOFT SPOTS (Top {top_n} - SMARTCyp-Inspired):"]
        for site in top_sites:
            # Map to no-H index for display
            display_idx = atom_map.get(site.atom_idx, site.atom_idx)
            legend_lines.append(
                f"• Rank #{site.rank}: Atom {display_idx} ({site.atom_symbol}) - "
                f"Score {site.score:.1f} - {site.reaction_type} - {site.cyp_isoform}"
            )
        legend_text = "\n".join(legend_lines)

        # Draw molecule without hydrogens
        drawer.DrawMolecule(
            mol_no_h,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=atom_colors,
            highlightBonds=list(bond_colors.keys()),
            highlightBondColors=bond_colors,
            legend=legend_text
        )

        drawer.FinishDrawing()
        drawer.WriteDrawingText(output_path)

        print(f"\n✓ Visualization saved to: {output_path}")
        return output_path


def compare_predictions(smiles: str, output_prefix: str = "comparison"):
    """
    Helper function to compare simple pattern matching vs quantitative scoring
    """
    print("\n" + "=" * 90)
    print("COMPARISON: Simple Pattern Matching vs. SMARTCyp-Inspired Scoring")
    print("=" * 90)

    # Run SMARTCyp-inspired prediction
    predictor = SMARTCypInspiredPredictor(smiles)
    sites = predictor.predict()

    print(f"\nFound {len(sites)} potential metabolic sites")
    print("\nTop 5 Sites:")
    for site in sites[:5]:
        print(f"  Rank {site.rank}: Atom #{site.atom_idx} - Score: {site.score:.1f} - {site.reaction_type}")

    # Generate visualization
    predictor.visualize(f"{output_prefix}_smartcyp_inspired.png")

    return predictor


if __name__ == "__main__":
    # Test with the peptidomimetic
    test_smiles = "CC(C)[C@@H](C(=O)N1Cc2ccccc2[C@H]1C(=O)Nc1cn(C)c(=O)n(C)c1=O)NC(=O)[C@H](CC1CCCCC1)NC(=O)C"

    print("╔" + "=" * 88 + "╗")
    print("║" + " " * 20 + "SMARTCYP-INSPIRED METABOLISM PREDICTOR" + " " * 29 + "║")
    print("║" + " " * 25 + "Quantitative Site Scoring" + " " * 38 + "║")
    print("╚" + "=" * 88 + "╝")

    predictor = SMARTCypInspiredPredictor(test_smiles)
    predictor.predict()

    # Print report
    print(predictor.generate_report())

    # Generate visualization
    predictor.visualize("smartcyp_inspired_prediction.png", top_n=5)

    print("\n✓ Analysis complete!")
