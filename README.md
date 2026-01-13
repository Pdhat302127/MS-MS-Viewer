# MS-MS Viewer: Peptide Fragment Ion Annotation Tool

## Overview
**MS-MS Viewer** is a Python-based command-line tool for extracting, visualizing, and annotating MS/MS spectra from mzXML (or gzipped mzXML) files.  
It annotates experimentally observed spectra with theoretical b- and y-ion fragment peaks derived from a given peptide sequence.



---

## Key Features
- Supports mzXML and `.gz` compressed mzXML files  
- Extracts peaks from a specific scan number  
- Calculates theoretical b- and y-ion masses  
- Matches ions using an m/z tolerance  
- Annotates MS/MS spectra with color-coded fragment ions  
- Filters low-intensity noise automatically  
- Produces publication-ready matplotlib plots  

---

## How It Works
1. Parses an mzXML file using streaming XML parsing (memory efficient)
2. Extracts m/z and intensity values from the specified scan
3. Computes theoretical b-ions and y-ions from the peptide sequence
4. Matches observed peaks to theoretical ions within a tolerance
5. Visualizes annotated and unannotated peaks

---

## Example Output
- **Blue peaks** → b-ions  
- **Red peaks** → y-ions  
- **Grey peaks** → unmatched peaks  

Annotated ions are labeled directly on the spectrum.
<img width="1280" height="612" alt="Figure_1" src="https://github.com/user-attachments/assets/4894852b-e40b-4fec-850f-e2c235226890" />

## Installation

### Requirements
- Python 3.9+
- matplotlib

### Clone the repository
```bash
git clone https://github.com/Pdhat302127/MS-MS-Viewer.git
cd MS-MS-Viewer
````

**### Usage**
Run the script from the command line:

```bash
python msms_spectrum_annotator.py <mzxml_file> <scan_number> <peptide_sequence>
````

## Input Arguments

| Argument          | Description                                   |
|------------------|-----------------------------------------------|
| `mzxml_file`      | Path to mzXML or `.gz` file                   |
| `scan_number`     | Scan number to extract                         |
| `peptide_sequence`| Peptide sequence (one-letter amino acid codes)|

---

## Algorithm Details

### Fragment Ion Calculation
- **b-ions**: cumulative N-terminal masses + proton  
- **y-ions**: cumulative C-terminal masses + H₂O + proton  
- Uses monoisotopic amino acid masses  

### Peak Matching
- m/z tolerance: ±0.05  
- Intensity threshold: ≥ 5% of maximum intensity  
- Nearest matching peak within tolerance is selected  

---

## Error Handling
The tool provides clear, user-friendly error messages for:
- Invalid peptide sequences  
- Missing scan numbers  
- File parsing errors  
- Missing or incorrect input arguments  

---

## Future Improvements
- Support for mzML files  
- Charge-state aware ion annotation  
- Neutral loss annotations  
- GUI or web-based interface


