# Met-ID SmartCyp App

**Advanced Metabolic Liability Prediction Tool with Quantitative Site Ranking**

A hybrid metabolism prediction application that combines simple pattern matching with SMARTCyp-inspired quantitative scoring to identify and rank metabolic soft spots in drug molecules.

![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)
![RDKit](https://img.shields.io/badge/RDKit-Required-green.svg)
![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)

---

## 🎯 Overview

Met-ID SmartCyp App provides three complementary approaches to metabolism prediction:

1. **Simple Pattern Matching** - Fast, qualitative detection of known metabolic liabilities
2. **SMARTCyp-Inspired Scoring** - Quantitative ranking with activation energies and accessibility scores
3. **Hybrid Comparison** - Side-by-side evaluation of both methods

### Key Features

- ✅ **Quantitative Site Ranking** - Metabolic sites scored 0-100 and ranked
- ✅ **Activation Energy Estimates** - DFT-inspired energy calculations
- ✅ **Steric Accessibility** - Automated accessibility analysis
- ✅ **CYP Isoform Specificity** - Predictions for CYP3A4, 2D6, 2C9
- ✅ **Dark Theme GUI** - Modern, user-friendly interface
- ✅ **Method Comparison** - Evaluate simple vs. advanced predictions
- ✅ **No Java Required** - Pure Python/RDKit implementation

---

## 📸 Screenshots

### Hybrid GUI
*(Dark-themed interface for comparing prediction methods)*

### Quantitative Prediction Output
```
TOP METABOLIC SITES (Ranked by Score)

Rank   Atom       Score    Energy     Access   Reaction                  CYP
------------------------------------------------------------------------------------------
1      #0    (C ) 59.4     90.0       0.99     ω-Oxidation               CYP3A4
2      #9    (C ) 39.2     85.0       0.56     Aromatic Hydroxylation    CYP1A2/2D6
8      #20   (N ) 36.0     75.0       0.40     N-Demethylation           CYP3A4/2D6
```

---

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/yourusername/Met-ID_SmartCyp_App.git
cd Met-ID_SmartCyp_App

# Install dependencies
pip install rdkit matplotlib pillow pandas

# Or using conda
conda install -c conda-forge rdkit matplotlib pillow pandas
```

### Launch GUI

```bash
python3 hybrid_metabolism_gui.py
```

### Command Line Usage

```python
from smartcyp_inspired_predictor import SMARTCypInspiredPredictor

# Predict metabolism for a SMILES string
smiles = "c1ccccc1CC"  # Ethylbenzene
predictor = SMARTCypInspiredPredictor(smiles)
sites = predictor.predict()

# View top 5 sites
for site in predictor.get_top_sites(n=5):
    print(f"Rank {site.rank}: Atom {site.atom_idx}")
    print(f"  Score: {site.score:.1f}")
    print(f"  Reaction: {site.reaction_type}")
    print(f"  CYP: {site.cyp_isoform}\n")
```

---

## 📁 Project Structure

```
Met-ID_SmartCyp_App/
├── README.md                          # This file
├── QUICK_START_HYBRID.md             # Quick start guide
├── HYBRID_PREDICTOR_README.md        # Full documentation
│
├── hybrid_metabolism_gui.py          # Main GUI application
├── smartcyp_inspired_predictor.py    # Quantitative predictor
├── metabolic_liability_tagger.py     # Pattern matching tool
│
├── requirements.txt                   # Python dependencies
├── .gitignore                        # Git ignore rules
└── LICENSE                           # Apache 2.0 license
```

---

## 🔬 Methodology

### Pattern Matching Approach
- Uses SMARTS patterns to detect known metabolic liabilities
- Binary detection (present/absent)
- Fast screening for compound libraries

### SMARTCyp-Inspired Quantitative Scoring
Inspired by the peer-reviewed SMARTCyp methodology:

> **Rydberg P., Gloriam D.E., Zaretzki J., Breneman C., Olsen L.**
> SMARTCyp: A 2D Method for Prediction of Cytochrome P450-Mediated Drug Metabolism.
> *ACS Med Chem Lett.* 2010;1(3):96-100.

**Our Implementation:**
- SMARTS pattern matching for site identification
- Approximate activation energies from literature
- Steric accessibility calculations
- Quantitative scoring (0-100)
- Ranked output

**Note:** This is an educational/screening tool. For production use, we recommend the official SMARTCyp server (https://smartcyp.sund.ku.dk) which uses DFT-calculated energies.

---

## 📊 Comparison: Pattern vs. Quantitative

| Feature | Pattern Matching | Quantitative Scoring |
|---------|-----------------|---------------------|
| **Output** | Binary (yes/no) | Ranked scores (0-100) |
| **Speed** | Very Fast ⚡⚡⚡ | Fast ⚡⚡ |
| **Ranking** | None | Top 10 sites |
| **Energies** | Not included | Activation energies |
| **Accessibility** | Not calculated | Steric scores |
| **Use Case** | Filtering | Prioritization |

---

## 💡 Use Cases

### Drug Discovery Pipeline

1. **Hit Identification** - Use pattern matching to filter virtual libraries
2. **Lead Optimization** - Use quantitative scoring to prioritize modifications
3. **SAR Analysis** - Compare analogs and rank metabolic liabilities
4. **Experimental Planning** - Guide microsomal stability assays

### Example Workflow

```python
# Screen library
from metabolic_liability_tagger import MetabolicLiabilityTagger

for compound in compound_library:
    tagger = MetabolicLiabilityTagger(compound.smiles)
    matches = tagger.analyze()

    # Filter: Keep only if < 3 high-risk sites
    high_risk_count = len(matches.get("Ester Hydrolysis", [])) + \
                      len(matches.get("CYP N-Dealkylation", []))

    if high_risk_count < 3:
        selected_compounds.append(compound)

# Rank survivors
from smartcyp_inspired_predictor import SMARTCypInspiredPredictor

for compound in selected_compounds:
    predictor = SMARTCypInspiredPredictor(compound.smiles)
    top_site = predictor.get_top_sites(n=1)[0]
    compound.metabolic_score = top_site.score

# Sort by metabolic stability (lower score = more stable)
selected_compounds.sort(key=lambda x: x.metabolic_score)
```

---

## ⚠️ Limitations

1. **Approximate Energies** - Not from quantum calculations
2. **Pattern-Based** - May miss novel metabolic sites
3. **2D Only** - No 3D conformational effects
4. **Not Validated** - Accuracy not determined experimentally
5. **Screening Tool** - Best for early discovery, not late-stage development

### When to Use More Rigorous Tools

- Late-stage drug development
- Regulatory submissions
- Need validated predictions → Use official SMARTCyp server
- Require experimental metabolite ID → LC-MS/MS studies

---

## 📚 Documentation

- **[QUICK_START_HYBRID.md](QUICK_START_HYBRID.md)** - Get started in 30 seconds
- **[HYBRID_PREDICTOR_README.md](HYBRID_PREDICTOR_README.md)** - Full technical documentation
- **In-code comments** - Detailed explanations in all modules

---

## 🤝 Contributing

Contributions welcome! Areas for improvement:

- [ ] Add more SMARTS patterns
- [ ] Refine activation energy estimates with experimental data
- [ ] Implement 3D conformational analysis
- [ ] Add batch processing GUI
- [ ] Validate against known metabolites
- [ ] Add export to CSV/PDF

---

## 📖 References

### SMARTCyp Publications

1. **Original SMARTCyp**:
   Rydberg P., et al. *ACS Med Chem Lett.* 2010;1(3):96-100.
   DOI: 10.1021/ml100016x

2. **SMARTCyp 3.0**:
   Olsen L., et al. *Bioinformatics* 2019;35(17):3174-3175.
   DOI: 10.1093/bioinformatics/btz037

### Metabolic Stability

3. **Drug-like Properties**:
   Kerns E.H., Di L. Academic Press, 2008.

4. **Metabolism in Drug Design**:
   Smith D.A., et al. *Wiley-VCH*, 2007.

---

## 📧 Support

- **Issues**: Open an issue on GitHub
- **Questions**: Check documentation first
- **SMARTCyp Official**: https://smartcyp.sund.ku.dk

---

## 📄 License

Apache License 2.0 - See [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **SMARTCyp Team** - For the original methodology
- **RDKit Community** - For the excellent cheminformatics toolkit
- **MDStudio Project** - For inspiration from their SMARTCyp wrapper

---

## 📊 Citation

If you use this tool in your research, please cite:

```bibtex
@software{metid_smartcyp_app,
  title = {Met-ID SmartCyp App: Hybrid Metabolic Liability Prediction Tool},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/yourusername/Met-ID_SmartCyp_App}
}
```

And cite the original SMARTCyp methodology:

```bibtex
@article{smartcyp2010,
  title={SMARTCyp: A 2D Method for Prediction of Cytochrome P450-Mediated Drug Metabolism},
  author={Rydberg, Patrik and Gloriam, David E and Zaretzki, Jed and Breneman, Curt and Olsen, Lars},
  journal={ACS Medicinal Chemistry Letters},
  volume={1},
  number={3},
  pages={96--100},
  year={2010},
  publisher={ACS Publications}
}
```

---

**Version**: 1.0.0
**Last Updated**: 2025-12-15
**Status**: Beta - Ready for screening and educational use

---

⭐ **If you find this tool useful, please star the repository!**
