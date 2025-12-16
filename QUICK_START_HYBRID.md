# Quick Start: Hybrid Metabolism Predictor

## ✨ What You Got

I built you a **hybrid metabolism prediction tool** that combines:

1. **Your Original Tool** (Simple pattern matching)
2. **SMARTCyp-Inspired Predictor** (Quantitative scoring with rankings)
3. **Comparison GUI** (See both side-by-side)

---

## 🚀 Launch in 30 Seconds

### Option 1: GUI (Recommended)

```bash
cd /Users/aardeshiri/chemistry_models
python3 hybrid_metabolism_gui.py
```

**What you'll see:**
- Dark-themed interface (like your existing tool)
- Three method options:
  - ⚡ Simple Pattern Matching
  - 🎯 SMARTCyp-Inspired Scoring ⭐ **NEW**
  - 📊 Both Methods (Compare)

**Try it:**
1. Click "LOAD EXAMPLE" (loads your peptidomimetic)
2. Select "Both Methods"
3. Click "RUN PREDICTION"
4. View ranked results with scores!

### Option 2: Command Line (Test Quantitative Predictor)

```bash
python3 smartcyp_inspired_predictor.py
```

Outputs:
- Ranked metabolic sites (1-10)
- Quantitative scores (0-100)
- Activation energies (kcal/mol)
- Steric accessibility (0-1)
- Visualization: `smartcyp_inspired_prediction.png`

---

## 🎯 Key Differences from Your Original Tool

### Your Original Tool:
```
✓ CYP N-Dealkylation: Found 2 match(es)
  Atom indices: ((20, 21), (24, 25))
✓ Benzylic Oxidation: Found 1 match(es)
  Atom indices: ((7, 8),)
```
**Output**: Binary (present/absent), no ranking

### New Quantitative Predictor:
```
Rank   Atom       Score    Energy     Access   Reaction                  CYP
------------------------------------------------------------------------------------------
1      #0    (C ) 59.4     90.0       0.99     ω-Oxidation               CYP3A4
8      #20   (N ) 36.0     75.0       0.40     N-Demethylation           CYP3A4/2D6
9      #24   (N ) 36.0     75.0       0.40     N-Demethylation           CYP3A4/2D6
```
**Output**: Ranked (1st most likely → 10th), with scores and energies

---

## 📊 Understanding the Scores

### Score Interpretation:

| Your Peptidomimetic Results |
|----------------------------|

**Top 3 Sites:**
1. **Atom #0 (59.4)** - Terminal methyl, accessible but LOW energy liability
2. **Atoms #9-#12 (39.2)** - Aromatic hydrogens, medium liability
3. **Atom #20, #24 (36.0)** - N-Methyls, HIGH energy liability despite lower score

**Key Insight:**
- **Higher score** = more accessible to metabolism
- **Lower energy** = more reactive when accessed
- **Both matter!** N-methyls (#20, #24) have low accessibility but HIGH reactivity

---

## 💡 When to Use Each Method

### Use Simple Pattern Matching:
- Quick yes/no: "Does this compound have metabolic liabilities?"
- Filtering large libraries
- Teaching/learning about soft spots

### Use Quantitative Scoring:
- "Which sites should I modify first?"
- Ranking multiple analogs
- Prioritizing SAR modifications
- Hypothesis generation for experiments

### Use Both:
- Comparing prediction approaches
- Validating consistency
- Educational purposes

---

## 🧪 Example Python API

```python
from smartcyp_inspired_predictor import SMARTCypInspiredPredictor

# Your molecule
smiles = "CC(C)[C@@H](C(=O)N1Cc2ccccc2[C@H]1C(=O)Nc1cn(C)c(=O)n(C)c1=O)NC(=O)[C@H](CC1CCCCC1)NC(=O)C"

# Predict
predictor = SMARTCypInspiredPredictor(smiles)
sites = predictor.predict()

# Get top 5
for site in sites[:5]:
    print(f"Rank {site.rank}: Atom {site.atom_idx}")
    print(f"  Score: {site.score} (0-100, higher = more likely)")
    print(f"  Reaction: {site.reaction_type}")
    print(f"  CYP Isoform: {site.cyp_isoform}")
    print()

# Generate visualization
predictor.visualize("my_output.png", top_n=5)
```

---

## 📁 New Files Created

```
/Users/aardeshiri/chemistry_models/

NEW FILES:
├── smartcyp_inspired_predictor.py      # Quantitative scoring engine
├── hybrid_metabolism_gui.py            # GUI for comparison
├── HYBRID_PREDICTOR_README.md          # Full documentation
└── QUICK_START_HYBRID.md              # This file

EXISTING (UNCHANGED):
├── metabolic_liability_tagger.py       # Your original tool
└── soft_spot_analyzer_with_legend.py  # Your original GUI

GENERATED OUTPUT:
├── smartcyp_inspired_prediction.png    # Ranked visualization
├── pattern_matching_result.png         # From original tool
└── quantitative_scoring_result.png     # From hybrid GUI
```

---

## 🔬 What's Different from Real SMARTCyp?

**Real SMARTCyp** (peer-reviewed, 76% accuracy):
- Uses DFT quantum calculations for energies
- Requires Java runtime
- Web server: https://smartcyp.sund.ku.dk

**This Tool** (educational/screening):
- Uses approximate energies from literature
- Pure Python/RDKit (no Java)
- Fast, local, easy to customize

**When to use Real SMARTCyp:**
- Production drug discovery
- Need validated predictions
- Publishing research

**When to use This Tool:**
- Early-stage screening
- Learning about metabolism
- Fast hypothesis generation
- Batch processing

---

## ⚡ Quick Tests

### Test 1: Simple Molecule
```bash
python3 -c "
from smartcyp_inspired_predictor import SMARTCypInspiredPredictor
pred = SMARTCypInspiredPredictor('c1ccccc1CC')  # Ethylbenzene
sites = pred.predict()
print(f'Top site: Atom {sites[0].atom_idx}, Score: {sites[0].score}')
"
```

### Test 2: Compare Methods (GUI)
```bash
python3 hybrid_metabolism_gui.py
# Click "LOAD EXAMPLE" → Select "Both Methods" → "RUN PREDICTION"
```

---

## 🐛 Troubleshooting

**GUI doesn't start:**
```bash
python3 -m tkinter  # Test if tkinter works
```

**Missing RDKit:**
```bash
pip install rdkit
# or
conda install -c conda-forge rdkit
```

**Import errors:**
```bash
# Make sure you're in the correct directory
cd /Users/aardeshiri/chemistry_models
python3 hybrid_metabolism_gui.py
```

---

## 📚 Learn More

- **Full Documentation**: `HYBRID_PREDICTOR_README.md`
- **Original Tool Docs**: `README_METABOLIC_TAGGER.md`
- **SMARTCyp Papers**: See references in `HYBRID_PREDICTOR_README.md`

---

## 🎉 Summary

You now have **3 metabolism prediction tools**:

| Tool | Method | Output | Speed | Use Case |
|------|--------|--------|-------|----------|
| **metabolic_liability_tagger** | Pattern | Binary | ⚡⚡⚡ | Screening |
| **smartcyp_inspired_predictor** | Quantitative | Ranked 0-100 | ⚡⚡ | Prioritization |
| **hybrid_metabolism_gui** | Both | Comparison | ⚡⚡ | Analysis |

**Recommended Workflow:**
1. Use **Pattern Matching** to filter libraries
2. Use **Quantitative Scoring** to rank survivors
3. Validate top-ranked sites experimentally

**Try it now:**
```bash
python3 hybrid_metabolism_gui.py
```

Enjoy! 🚀
