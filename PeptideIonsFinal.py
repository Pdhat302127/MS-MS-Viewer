class PeptideIons:

    def __init__(self, seq=""):
        self.seq=seq.upper()
        
    def amino_acids(self):
        aa_masses = {
            'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,
            'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,
            'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
            'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333
        }
        return aa_masses

    def length(self):
        return len(self.seq)
        
    def compute_b_ions(self): # Compute b ions (N-terminal fragments)
        peptide_length = self.length()
        b_ions = {}
    
        b_ion=1
        for i in range(0,peptide_length-1):
            aa=self.seq[i]
            if aa in self.amino_acids():
                b_ion=self.amino_acids().get(aa)+b_ion
                b_ions["b"+str(i+1)]=b_ion
            else:
                raise ValueError("amino acid:"+aa+" "+"not present in amino acid dictionary")
        return b_ions
    
    def compute_y_ions(self): # Compute y ions (C-terminal fragments)
        peptide_length = self.length()
        y_ions = {}
    
        y_ion=19
        for i in range(peptide_length-1,0,-1):
            aa=self.seq[i]
            if aa in self.amino_acids():
                y_ion=self.amino_acids().get(aa)+y_ion
                y_ions["y"+str(peptide_length-i)]=y_ion
            else:
                raise ValueError("amino acid:"+aa+" "+"not present in amino acid dictionary")
        return y_ions
    
    def ions(self):
        b_ions=self.compute_b_ions()
        y_ions=self.compute_y_ions()

        b_ions.update(y_ions)
        return b_ions





        
        

