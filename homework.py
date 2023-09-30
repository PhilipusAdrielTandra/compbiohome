inp = "ttacga"

def reverse_complement(dna_seq):
    comp = ""
    for c in dna_seq.upper():
        if c == "T":
            comp = comp + "A"
        elif c == "A":
            comp = comp + "T"
        elif c == "C":
            comp = comp + "G"
        elif c == "G":
            comp = comp + "C" 
    return comp

#transcripts dna to rna
def mRNA(dna_seq):
    comp = ""
    for c in dna_seq.upper():
        if c == "T":
            comp = comp + "A"
        elif c == "A":
            comp = comp + "U"
        elif c == "C":
            comp = comp + "G"
        elif c == "G":
            comp = comp + "C" 
    return comp

def translate_codon(cod):
    #codons and their corresponding amino acids
    tc = {
        "UUU" : "Phe - (F)", "UUC" : "Phe - (F)", "UUA" : "Leu - (L)", "UUG" : "Leu - (L)",
        "UCU" : "Ser - (S)", "UCC" : "Ser - (S)", "UCA" : "Ser - (S)", "UCG" : "Ser - (S)",
        "UAU" : "Tyr - (Y)", "UAC" : "Tyr - (Y)", "UAA" : "_", "UAG" : "_", 
        "UGU" : "Cys - (C)", "UGC" : "Cys - (C)", "UGA" : "_", "UGG" : "Trp - (W)", 

        "CUU" : "Leu - (L)", "CUC" : "Leu - (L)", "CUA" : "Leu - (L)", "CUG" : "Leu - (L)", 
        "CCU" : "Pro - (P)", "CCC" : "Pro - (P)", "CCA" : "Pro - (P)", "CCG" : "Pro - (P)", 
        "CAU" : "His - (H)", "CAC" : "His - (H)", "CAA" : "Gln - (Q)", "CAG" : "Gln - (Q)",
        "CGU" : "Arg - (R)", "CGC" : "Arg - (R)", "CGA" : "Arg - (R)", "CGG" : "Arg - (R)", 

        "AUU" : "Ile - (I)", "AUC" : "Ile - (I)", "AUA" : "Ile - (I)", "AUG" : "Met - (M)", 
        "ACU" : "Thr - (T)", "ACC" : "Thr - (T)", "ACA" : "Thr - (T)", "ACG" : "Thr - (T)",
        "AAU" : "Asn - (N)", "AAC" : "Asn - (N)", "AAA" : "Lys - (K)", "AAG" : "Lys - (K)",
        "AGU" : "Ser - (S)", "AGC" : "Ser - (S)", "AGA" : "Arg - (R)", "AGG" : "Arg - (R)", 

        "GUU" : "Val - (V)", "GUC" : "Val - (V)", "GUA" : "Val - (V)", "GUG" : "Val - (V)", 
        "GCU" : "Ala - (A)", "GCC" : "Ala - (A)", "GCA" : "Ala - (A)", "GCG" : "Ala - (A)", 
        "GAU" : "Asp - (D)", "GAC" : "Asp - (D)", "GAA" : "Glu - (E)", "GAG" : "Glu - (E)", 
        "GGU" : "Gly - (G)", "GGC" : "Gly - (G)", "GGA" : "Gly - (G)", "GGG" : "Gly - (G)", 

    }

    #create new array split into threes
    cod = cod.upper()
    new = []
    while cod:
        new.append(cod[:3])
        cod = cod[3:]
    translation = []
    # we can simply check if the three characters are valid in the dictionary above with their amino acids
    for codon in new:
        if codon in tc:
            translation.append(tc[codon])
        else:
            translation.append("Unknown")

    return translation

def amino_acid_sequence_to_rna(x):
    #amino acids and their corresponding codons
    codon_table = {
        'A': ['GCU', 'GCC', 'GCA', 'GCG'],
        'C': ['UGU', 'UGC'],
        'D': ['GAU', 'GAC'],
        'E': ['GAA', 'GAG'],
        'F': ['UUU', 'UUC'],
        'G': ['GGU', 'GGC', 'GGA', 'GGG'],
        'H': ['CAU', 'CAC'],
        'I': ['AUU', 'AUC', 'AUA'],
        'K': ['AAA', 'AAG'],
        'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
        'M': ['AUG'],
        'N': ['AAU', 'AAC'],
        'P': ['CCU', 'CCC', 'CCA', 'CCG'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
        'T': ['ACU', 'ACC', 'ACA', 'ACG'],
        'V': ['GUU', 'GUC', 'GUA', 'GUG'],
        'W': ['UGG'],
        'Y': ['UAU', 'UAC'],
        '*': ['UAA', 'UAG', 'UGA'], 
    }

    dna_sequences = ['']

    #append to array if input is valid to the table above
    for i in x:
        if i not in codon_table:
            return ["Invalid amino acid"]

        new_dna_sequences = []
        for codon in codon_table[i]:
            for dna_sequence in dna_sequences:
                new_dna_sequences.append(dna_sequence + codon)
        
        dna_sequences = new_dna_sequences

    return dna_sequences

def freq(seq):
    #similar to the freq function in class, but this time we want to split it into threes to see freq of codons since codons are in threes
    sequence_counts = {}
    #iterate over threes of the string, if already in dictionary then add to the count, if not add a new key and set initial value to one
    for i in range(0, len(seq) - 2, 3):
        sequence = seq[i:i+3]
        sequence_counts[sequence] = sequence_counts.get(sequence, 0) + 1

    return sequence_counts

#printing and calling functions to format
var = mRNA(inp)
print("Input DNA =", inp, "\n")
print("Complement =", reverse_complement(inp))
print("mRNA =", var)
print("Amino Acid =", end=" ")
for i in translate_codon(var):
    print(i, end = " ")
print()

inputa = 'NAN'
print("Input Aminoacid =", inputa, "\n")
dna_sequences = amino_acid_sequence_to_rna(inputa)
print("mRNA =", dna_sequences[0])
result = freq(dna_sequences[0])
for i,count in result.items():
    print(f"{i} = {count}")