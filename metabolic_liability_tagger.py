#!/usr/bin/env python3
"""
Metabolic Liability Tagger - Automated Soft Spot Analysis Tool

This script performs automated metabolic soft spot analysis on drug molecules
using SMARTS pattern matching to identify sites of metabolic vulnerability.

Author: Senior Python Developer & Chemoinformatics Expert
Date: 2025-12-15
Requirements: rdkit, matplotlib, PIL
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from typing import Dict, List, Tuple, Optional
import sys
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class LiabilityRule:
    """Data class to store metabolic liability rules"""
    name: str
    smarts: str
    color: Tuple[float, float, float]  # RGB tuple (0-1 range)
    description: str
    risk_level: str


class MetabolicLiabilityTagger:
    """
    A tool for identifying and visualizing metabolic soft spots in drug molecules.

    Uses SMARTS pattern matching to detect common sites of metabolic transformation
    including CYP-mediated oxidation, hydrolytic cleavage, and conjugation sites.
    """

    # Define metabolic liability rules with SMARTS patterns
    LIABILITY_RULES = [
        LiabilityRule(
            name="CYP N-Dealkylation",
            smarts="[n][CH3]",  # N-methyl on aromatic nitrogen
            color=(1.0, 0.2, 0.2),  # Red
            description="N-Methyl groups on heterocycles (CYP3A4/2D6)",
            risk_level="High"
        ),
        LiabilityRule(
            name="CYP N-Dealkylation (Secondary)",
            smarts="[NX3;!R][CH3]",  # N-methyl on non-ring nitrogen
            color=(1.0, 0.4, 0.4),  # Light red
            description="N-Methyl on aliphatic amines",
            risk_level="Medium"
        ),
        LiabilityRule(
            name="Benzylic Oxidation",
            smarts="[CH2,CH;!R][c]",  # Benzylic carbon
            color=(1.0, 0.6, 0.0),  # Orange
            description="Benzylic positions (CYP2D6/3A4)",
            risk_level="High"
        ),
        LiabilityRule(
            name="Allylic Oxidation",
            smarts="[CH2,CH;!R][C]=[C]",  # Allylic carbon
            color=(1.0, 0.7, 0.2),  # Light orange
            description="Allylic positions (CYP-mediated)",
            risk_level="Medium"
        ),
        LiabilityRule(
            name="Amide Hydrolysis",
            smarts="[C;!R](=O)[NX3;!R]",  # Non-ring amide
            color=(0.3, 0.6, 1.0),  # Blue
            description="Amide bonds (peptidases/amidases)",
            risk_level="High"
        ),
        LiabilityRule(
            name="Ester Hydrolysis",
            smarts="[C](=O)[OX2][C]",  # Ester
            color=(0.4, 0.7, 1.0),  # Light blue
            description="Ester bonds (esterases)",
            risk_level="Very High"
        ),
        LiabilityRule(
            name="Aromatic Hydroxylation",
            smarts="c1ccccc1",  # Benzene ring
            color=(0.6, 0.3, 0.9),  # Purple
            description="Unsubstituted aromatic rings (CYP-mediated)",
            risk_level="Medium"
        ),
        LiabilityRule(
            name="Aliphatic Hydroxylation (Cyclohexyl)",
            smarts="C1CCCCC1",  # Cyclohexyl ring
            color=(0.9, 0.6, 0.2),  # Brown
            description="Cyclohexyl rings (ω-hydroxylation)",
            risk_level="Medium"
        ),
        LiabilityRule(
            name="O-Dealkylation",
            smarts="[OX2][CH3]",  # O-methyl ether
            color=(0.9, 0.4, 0.5),  # Pink
            description="O-Methyl ethers (CYP-mediated)",
            risk_level="Medium-High"
        ),
        LiabilityRule(
            name="Terminal Methyl (Aliphatic)",
            smarts="[CH3][CH;!R]([CH3])",  # Isopropyl pattern
            color=(0.7, 0.7, 0.3),  # Olive
            description="Branched aliphatic methyls (ω-oxidation)",
            risk_level="Low-Medium"
        ),

        # ========== PHASE II METABOLISM ==========

        LiabilityRule(
            name="Phenol (Glucuronidation/Sulfation)",
            smarts="[OH]c1ccccc1",  # Phenolic hydroxyl
            color=(0.0, 0.8, 0.8),  # Cyan
            description="Phenolic -OH (UGT/SULT substrate)",
            risk_level="Very High"
        ),
        LiabilityRule(
            name="Primary Alcohol (Glucuronidation)",
            smarts="[CH2][OH]",  # Primary alcohol
            color=(0.2, 0.9, 0.9),  # Light cyan
            description="Primary -OH (UGT substrate)",
            risk_level="High"
        ),
        LiabilityRule(
            name="Secondary Alcohol (Glucuronidation)",
            smarts="[CH;!R][OH]",  # Secondary alcohol, non-ring
            color=(0.3, 0.85, 0.85),  # Light cyan
            description="Secondary -OH (UGT substrate)",
            risk_level="Medium-High"
        ),
        LiabilityRule(
            name="Carboxylic Acid (Glucuronidation)",
            smarts="[CX3](=O)[OH]",  # Carboxylic acid
            color=(0.1, 0.7, 0.9),  # Blue-cyan
            description="Carboxylic acid (UGT substrate)",
            risk_level="High"
        ),
        LiabilityRule(
            name="Primary Amine (Acetylation)",
            smarts="[NX3;H2;!R]",  # Primary amine, non-ring
            color=(0.5, 0.0, 0.8),  # Purple
            description="Primary -NH₂ (NAT1/NAT2 substrate)",
            risk_level="Medium-High"
        ),
        LiabilityRule(
            name="Aromatic Amine (Acetylation)",
            smarts="[NH2]c1ccccc1",  # Aniline
            color=(0.6, 0.0, 0.9),  # Bright purple
            description="Aromatic -NH₂ (NAT substrate)",
            risk_level="High"
        ),
        LiabilityRule(
            name="Thiol (Glutathione Conjugation)",
            smarts="[SH]",  # Thiol group
            color=(0.9, 0.9, 0.0),  # Yellow
            description="Thiol -SH (GST substrate)",
            risk_level="Medium"
        ),
        LiabilityRule(
            name="Catechol (COMT Methylation)",
            smarts="[OH]c1ccc([OH])cc1",  # Catechol (ortho-dihydroxy)
            color=(0.0, 0.6, 0.6),  # Dark cyan
            description="Catechol (COMT O-methylation)",
            risk_level="High"
        ),
    ]

    def __init__(self, smiles: str):
        """
        Initialize the tagger with a SMILES string.

        Args:
            smiles: SMILES string of the molecule to analyze

        Raises:
            ValueError: If SMILES string is invalid
        """
        self.smiles = smiles
        self.mol = self._parse_smiles(smiles)
        self.liability_matches = {}

    def _parse_smiles(self, smiles: str) -> Chem.Mol:
        """
        Parse SMILES string and validate.

        Args:
            smiles: SMILES string

        Returns:
            RDKit Mol object

        Raises:
            ValueError: If SMILES is invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # Generate 2D coordinates for visualization
        AllChem.Compute2DCoords(mol)
        return mol

    def analyze(self) -> Dict[str, List[Tuple]]:
        """
        Perform metabolic liability analysis by matching SMARTS patterns.

        Returns:
            Dictionary mapping liability names to lists of atom index tuples
        """
        self.liability_matches = defaultdict(list)

        for rule in self.LIABILITY_RULES:
            pattern = Chem.MolFromSmarts(rule.smarts)

            if pattern is None:
                print(f"Warning: Invalid SMARTS pattern for {rule.name}: {rule.smarts}")
                continue

            matches = self.mol.GetSubstructMatches(pattern)

            if matches:
                self.liability_matches[rule.name] = matches
                print(f"✓ {rule.name}: Found {len(matches)} match(es)")
                print(f"  Risk: {rule.risk_level} | {rule.description}")
                print(f"  Atom indices: {matches}")

        if not self.liability_matches:
            print("No metabolic liabilities detected (or molecule is highly stable!)")

        return dict(self.liability_matches)

    def _prepare_highlights(self) -> Tuple[List[int], Dict[int, Tuple], List[int], Dict[int, Tuple]]:
        """
        Prepare atom and bond highlighting information.

        Returns:
            Tuple of (atom_list, atom_colors, bond_list, bond_colors)
        """
        atom_colors = {}
        bond_colors = {}

        # Create a mapping of rule names to colors
        rule_color_map = {rule.name: rule.color for rule in self.LIABILITY_RULES}

        # Priority order (higher priority overwrites lower)
        priority_order = [
            "Terminal Methyl (Aliphatic)",
            "Aromatic Hydroxylation",
            "Aliphatic Hydroxylation (Cyclohexyl)",
            "Allylic Oxidation",
            "O-Dealkylation",
            "CYP N-Dealkylation (Secondary)",
            "Benzylic Oxidation",
            "Amide Hydrolysis",
            "CYP N-Dealkylation",
            "Ester Hydrolysis",
        ]

        # Apply colors in priority order
        for rule_name in priority_order:
            if rule_name in self.liability_matches:
                color = rule_color_map.get(rule_name)
                for match in self.liability_matches[rule_name]:
                    for atom_idx in match:
                        atom_colors[atom_idx] = color

                    # Highlight bonds within matched substructure
                    for i in range(len(match)):
                        for j in range(i + 1, len(match)):
                            bond = self.mol.GetBondBetweenAtoms(match[i], match[j])
                            if bond:
                                bond_colors[bond.GetIdx()] = color

        atom_list = list(atom_colors.keys())
        bond_list = list(bond_colors.keys())

        return atom_list, atom_colors, bond_list, bond_colors

    def generate_report(self) -> str:
        """
        Generate a text report of findings.

        Returns:
            Formatted text report
        """
        report = []
        report.append("=" * 80)
        report.append("METABOLIC LIABILITY ANALYSIS REPORT")
        report.append("=" * 80)
        report.append(f"\nSMILES: {self.smiles}")
        report.append(f"Molecular Formula: {rdMolDescriptors.CalcMolFormula(self.mol)}")
        report.append(f"Molecular Weight: {rdMolDescriptors.CalcExactMolWt(self.mol):.2f} g/mol")
        report.append(f"Heavy Atoms: {self.mol.GetNumHeavyAtoms()}")

        report.append("\n" + "=" * 80)
        report.append("DETECTED LIABILITIES")
        report.append("=" * 80)

        if not self.liability_matches:
            report.append("\nNo metabolic liabilities detected.")
        else:
            for rule in self.LIABILITY_RULES:
                if rule.name in self.liability_matches:
                    matches = self.liability_matches[rule.name]
                    report.append(f"\n[{rule.risk_level}] {rule.name}")
                    report.append(f"  Description: {rule.description}")
                    report.append(f"  Matches: {len(matches)}")
                    report.append(f"  Atom indices: {matches}")

        report.append("\n" + "=" * 80)
        return "\n".join(report)

    def visualize(self, output_path: str = "soft_spot_analysis.png",
                  width: int = 1600, height: int = 1200) -> str:
        """
        Generate annotated molecular structure image with highlighted soft spots.

        Args:
            output_path: Path to save output image
            width: Image width in pixels
            height: Image height in pixels

        Returns:
            Path to saved image
        """
        if not self.liability_matches:
            print("No liabilities to visualize. Running analysis first...")
            self.analyze()

        # Prepare highlighting
        atom_list, atom_colors, bond_list, bond_colors = self._prepare_highlights()

        # Create drawer
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().clearBackground = True
        drawer.drawOptions().padding = 0.05
        drawer.drawOptions().bondLineWidth = 3
        drawer.drawOptions().atomLabelFontSize = 30
        drawer.drawOptions().legendFontSize = 20

        # Build legend text
        legend_lines = ["METABOLIC SOFT SPOTS:"]
        for rule in self.LIABILITY_RULES:
            if rule.name in self.liability_matches:
                count = len(self.liability_matches[rule.name])
                legend_lines.append(f"• {rule.name} ({count} sites) - {rule.risk_level}")

        legend_text = "\n".join(legend_lines)

        # Draw molecule
        drawer.DrawMolecule(
            self.mol,
            highlightAtoms=atom_list,
            highlightAtomColors=atom_colors,
            highlightBonds=bond_list,
            highlightBondColors=bond_colors,
            legend=legend_text
        )

        drawer.FinishDrawing()
        drawer.WriteDrawingText(output_path)

        print(f"\n✓ Visualization saved to: {output_path}")
        return output_path

    def run_full_analysis(self, output_path: str = "soft_spot_analysis.png") -> str:
        """
        Run complete analysis workflow: analyze, report, visualize.

        Args:
            output_path: Path for output image

        Returns:
            Path to generated image
        """
        print("\n" + "=" * 80)
        print("STARTING METABOLIC LIABILITY ANALYSIS")
        print("=" * 80 + "\n")

        # Step 1: Analyze
        self.analyze()

        # Step 2: Generate report
        report = self.generate_report()
        print("\n" + report)

        # Step 3: Visualize
        image_path = self.visualize(output_path)

        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE")
        print("=" * 80 + "\n")

        return image_path


def analyze_smiles(smiles: str, output_path: str = "soft_spot_analysis.png") -> Optional[str]:
    """
    Convenience function to analyze a SMILES string.

    Args:
        smiles: SMILES string to analyze
        output_path: Path for output image

    Returns:
        Path to generated image, or None if error
    """
    try:
        tagger = MetabolicLiabilityTagger(smiles)
        return tagger.run_full_analysis(output_path)
    except ValueError as e:
        print(f"Error: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return None


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Test SMILES: Ac-Cha-Val-Tic-Dimethyluracil peptidomimetic
    test_smiles = "CC(C)[C@@H](C(=O)N1Cc2ccccc2[C@H]1C(=O)Nc1cn(C)c(=O)n(C)c1=O)NC(=O)[C@H](CC1CCCCC1)NC(=O)C"

    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 20 + "METABOLIC LIABILITY TAGGER" + " " * 32 + "║")
    print("║" + " " * 15 + "Automated Soft Spot Analysis Tool" + " " * 30 + "║")
    print("╚" + "=" * 78 + "╝")

    # Option 1: Using the class directly
    print("\n[Method 1] Using MetabolicLiabilityTagger class:")
    tagger = MetabolicLiabilityTagger(test_smiles)
    tagger.run_full_analysis("soft_spot_analysis.png")

    # Option 2: Using convenience function
    # print("\n[Method 2] Using convenience function:")
    # analyze_smiles(test_smiles, "soft_spot_analysis_v2.png")

    print("\n✓ Analysis complete! Check the output images.")

    # Example: Analyze from command line argument
    # if len(sys.argv) > 1:
    #     user_smiles = sys.argv[1]
    #     output = sys.argv[2] if len(sys.argv) > 2 else "soft_spot_analysis.png"
    #     analyze_smiles(user_smiles, output)
