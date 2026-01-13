import sys
import gzip
import os
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from base64 import b64decode
from array import array

# Peak extraction
def peak_extraction(afile, ns, scan):
    if afile.endswith('.gz'):
        open_file = gzip.open(afile, 'rt')
    else:
        open_file = open(afile, 'rt')
    
    try:
        xml_document = ET.iterparse(open_file, events=("end",))
        for event, element in xml_document:
            if event == 'end' and element.tag == ns + 'scan':
                if int(element.attrib.get('num')) == scan:
                    peak = element.find(ns + 'peaks')
                    if peak is not None:
                        peaks = array('f', b64decode(peak.text))
                        if sys.byteorder != 'big':
                            peaks.byteswap()
                        mzs = peaks[::2]
                        ints = peaks[1::2]
                        return mzs, ints
                element.clear()
    except ET.ParseError as e:
        raise ValueError("Error parsing XML file: " + str(e))
    open_file.close()
    print("Scan number " + str(scan) + " not found in file " + afile + ". Please check the scan number and try again.")
    raise ValueError("Scan number " + str(scan) + " not found in file " + afile + ".")

# Calculation of b-ions and y-ions
def calculate_ion_masses(peptide_sequence):
    aa_masses = {
        'A': 71.037, 'C': 103.009, 'D': 115.027, 'E': 129.043, 'F': 147.068,
        'G': 57.021, 'H': 137.059, 'I': 113.084, 'K': 128.095, 'L': 113.084,
        'M': 131.040, 'N': 114.043, 'P': 97.053, 'Q': 128.059, 'R': 156.101,
        'S': 87.032, 'T': 101.048, 'V': 99.068, 'W': 186.079, 'Y': 163.063
    }

    b_ions, y_ions = [], []
    b_mass, y_mass = 0, 18

    for i, aa in enumerate(peptide_sequence):
        if aa not in aa_masses:
            raise ValueError("Error: Invalid amino acid '" + aa + "' in peptide sequence. Please provide a valid peptide sequence.")
        b_mass += aa_masses[aa]
        b_ions.append(b_mass + 1)
        y_mass += aa_masses[peptide_sequence[-(i+1)]]
        y_ions.append(y_mass + 1)

    return b_ions, y_ions

# Finding the nearest peak
def find_nearest_peak(mzs, target_mz, tolerance):
    nearest_position = None
    min_diff = tolerance + 1
    for position, mz in enumerate(mzs):
        diff = abs(mz - target_mz)
        if diff <= tolerance and diff < min_diff:
            min_diff = diff
            nearest_position = position
    return nearest_position

# Plotting spectrum
def plot_peptide_spectrum(mzs, ints, b_ions, y_ions, peptide_sequence):
    plt.figure(figsize=(12, 6))

    max_intensity = max(ints)
    intensity_threshold = 0.05 * max_intensity
    tolerance = 0.05

    annotated_indices = set()
    annotations = []

    # b-ions
    for i, b_ion in enumerate(b_ions):
        ion_label = f"b{i+1}"
        nearest_position = find_nearest_peak(mzs, b_ion, tolerance)

        if nearest_position is not None and ints[nearest_position] >= intensity_threshold:
            annotated_indices.add(nearest_position)
            annotations.append(
                (mzs[nearest_position], ints[nearest_position], ion_label, 'blue')
            )
            print(
                f"{ion_label}: m/z={mzs[nearest_position]}, "
                f"intensity={ints[nearest_position]}, "
                f"Δ={abs(mzs[nearest_position] - b_ion):.4f}"
            )

    #y-ions 
    for i, y_ion in enumerate(y_ions):
        ion_label = f"y{i+1}"
        nearest_position = find_nearest_peak(mzs, y_ion, tolerance)

        if nearest_position is not None and ints[nearest_position] >= intensity_threshold:
            annotated_indices.add(nearest_position)
            annotations.append(
                (mzs[nearest_position], ints[nearest_position], ion_label, 'red')
            )
            print(
                f"{ion_label}: m/z={mzs[nearest_position]}, "
                f"intensity={ints[nearest_position]}, "
                f"Δ={abs(mzs[nearest_position] - y_ion):.4f}"
            )

    # Plot unannotated peaks 
    unannotated_mzs = [mzs[i] for i in range(len(mzs)) if i not in annotated_indices]
    unannotated_ints = [ints[i] for i in range(len(ints)) if i not in annotated_indices]

    plt.stem(
        unannotated_mzs,
        unannotated_ints,
        linefmt='grey',
        markerfmt=' ',
        basefmt=" ",
        label='Unmatched peaks'
    )

    # Plot annotated peaks 
    for mz, intensity, label, color in annotations:
        plt.stem([mz], [intensity], linefmt=color, markerfmt=' ', basefmt=" ")
        plt.text(mz, intensity + 0.03 * max_intensity, label, color=color, ha='center')

    # Legend
    plt.stem([0], [0], linefmt='b-', markerfmt=" ", basefmt=" ", label="b-ions")
    plt.stem([0], [0], linefmt='r-', markerfmt=" ", basefmt=" ", label="y-ions")

    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title(f"Annotated Spectrum for Peptide: {peptide_sequence}")
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.show()


# Main function
try:

    if len(sys.argv) != 4:
        raise ValueError("Usage: python project.py <file_name> <scan_number> <peptide_sequence>")

    file_name = sys.argv[1]
    scan_number = int(sys.argv[2])
    peptide_sequence = sys.argv[3].upper()

    if not os.path.exists(file_name):
        print("FileNotFound: File '" + file_name + "' not found.")
        sys.exit(1)
    
    namespace = '{http://sashimi.sourceforge.net/schema/}'
    mzs, ints = peak_extraction(file_name, namespace, scan_number)

    if not mzs or not ints:
        print("Scan number " + str(scan_number) + " not found in file " + file_name + 
              ". Please check the scan number and try again.")
        sys.exit(1)
    else:
        try:
            b_ions, y_ions = calculate_ion_masses(peptide_sequence)
            if b_ions and y_ions:
                plot_peptide_spectrum(mzs, ints, b_ions, y_ions, peptide_sequence)
        except ValueError as e:
            print(e)
            sys.exit(1)

except Exception as e:
    print("An unexpected error occurred: " + str(e))
    sys.exit(1)
