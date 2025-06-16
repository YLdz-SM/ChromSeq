# ChromSeq - Sanger Sequence Editing Tool


## Overview
ChromSeq is a Python application designed to replace ChromasPro for basic Sanger sequence editing tasks. It provides a simple interface for visualizing, editing, and generating consensus sequences from forward and reverse Sanger sequencing files.

## Key Features
- 🧬 Dual-sequence visualization (forward & reverse)
- ⚡ Automatic sequence alignment
- 🔍 Multiple consensus generation algorithms
- 📋 FASTA-formatted output
- 🛠️ Basic sequence editing capabilities

## Installation
1. Ensure you have Python 3.6+ installed
2. Clone this repository or download the source files
3. Install required dependencies:
4. Install Python 3.8+ from python.org
5. Install dependencies:

pip install biopython numpy matplotlib tk



QUICK START
-----------
1. Launch the app:
   python3 main.py
   (or drag main.py into your terminal)

2. Load sequences:
   - Click "Add Files"
   - Select forward/reverse .ab1 or .fasta files

3. Assign roles:
   - Select forward sequence → Click "Forward" (turns blue)
   - Select reverse sequence → Click "Reverse" (turns red)

4. Align sequences:
   - Click the "Align" button
   - Wait for confirmation message

5. Generate consensus:
   - Click "Consensus"
   - Select "Two-Step Consensus" (recommended)
   - Copy the FASTA output:

>YourSequenceName
AGCTAGCTNN...GTACGTAC

6. Verify ambiguous bases:
   - BLAST any N/M/W/K bases
   - Manually correct as needed

FEATURES
--------
- 5 consensus algorithms
- Quality score visualization
- Cross-platform (Win/macOS/Linux)
- FASTA export

TROUBLESHOOTING
---------------
Problem: Sequences won't align
Solution:
- Ensure ≥50bp overlap
- Trim low-quality ends

Problem: Many ambiguous bases
Solution:
- Try "Quality Threshold" algorithm
- Increase threshold to Q30

SUPPORT
-------
Email: your.email@example.com
GitHub: github.com/yourusername/ChromSeq/issues
