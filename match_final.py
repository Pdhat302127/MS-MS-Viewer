from PeaksFinal import *
from PeptideIonsFinal import *
import gzip

# Import sys is in peaks.py module
if len(sys.argv) < 3:
    print("Require a scan number and peptide sequence on command-line.", file=sys.stderr)
    sys.exit(1)
    
# Scan number and peptide sequence provided by the user
scan=str(sys.argv[1])

peptide=str(sys.argv[2])
seq=PeptideIons(peptide)

# Open the spectrum file:
try:
    file=("17mix_test2.mzxml.gz")
    afile=gzip.open(file, mode="rb")
except FileNotFoundError:
    print("The xml file could not be opened",file=sys.stderr)
    sys.exit(1)
    
ns='{http://sashimi.sourceforge.net/schema/}'

try:
    mzs, ints = peak_extraction(afile, ns, scan)
except ValueError as e:
    print(e,file=sys.stderr)
    sys.exit(1)

afile.close()

# Function for matching the ions
def match_ions(ions, mzs, ints):
    best_matching_ions={}
    highest_peak=max(ints)
    tolerance=0.5
    for key_ion, mz_peptide in ions.items():
        
        all_possibilities={}
        for i in range(0, len(mzs)):
            mz_peak = mzs[i]  
            intensity = ints[i]  
            difference = abs(mz_peptide - mz_peak)
            all_possibilities["peak"+str(i)] = {'m/z': mz_peak,'intensity': intensity,'difference': difference}
            
        #let's find the best match for this key_ion:
        for peak,values in all_possibilities.items():
            if values['intensity'] >= 0.01 * highest_peak:
                if values['difference']<=tolerance:
                    best_matching_ions[key_ion]={'m/z': values['m/z'],'intensity': values['intensity']}
    return best_matching_ions            


# To plot the peaks in a stem plot
import matplotlib.pyplot as plt

best_match=match_ions(seq.ions(), mzs, ints)

mz_b=[]
intensity_b=[]
ion_labels_b=[]
mz_y=[]
intensity_y=[]
ion_labels_y=[]


for key,value in best_match.items():
    if key.startswith("b"):
        mz_b.append(value["m/z"])
        intensity_b.append(value["intensity"])
        ion_labels_b.append(key)
    else:
        mz_y.append(value["m/z"])
        intensity_y.append(value["intensity"])
        ion_labels_y.append(key)


plt.figure(figsize=(10,6))
plt.stem(mzs, ints, markerfmt=" ", basefmt='k-', linefmt='k-')

# To add labels and colors to each matched peak
for i in range(len(mz_b)):
    plt.text(mz_b[i], intensity_b[i] + 1, ion_labels_b[i], ha='center', va='bottom', fontsize=12)
    plt.stem(mz_b[i], intensity_b[i], linefmt='b-',markerfmt=" ",basefmt=" ")
    
for i in range(len(mz_y)):
    plt.text(mz_y[i], intensity_y[i] + 1, ion_labels_y[i], ha='center', va='bottom', fontsize=12)
    plt.stem(mz_y[i], intensity_y[i], linefmt='c-',markerfmt=" ",basefmt=" ")            
    
# To describe each ion type in a legend
if mz_b and intensity_b:
    plt.stem(mz_b, intensity_b, linefmt='b-', markerfmt=" ", basefmt=" ", label="b Ions")

if mz_y and intensity_y:
    plt.stem(mz_y, intensity_y, linefmt='c-', markerfmt=" ", basefmt=" ", label="y Ions")


# Axis labels and title
plt.xlabel("m/z (Th)")
plt.ylabel("Intensity")
plt.title("MS/MS fragmentation for peptide:" + peptide + " " + "& scan:" + scan)
plt.legend(loc="upper right")

# To show the plot
plt.show()
#plt.savefig("test.png")


