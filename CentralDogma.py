## TO EVALUATE DICTIONARY CONCEPTS ##

# seq_dna            # dictionary
# seq_dna.keys()     # only keys
# seq_dna.values()   # only DNA sequences
# seq_dna.items()    # (key, value) pairs


# for seq, dna in seq_dna.items():
#     print(seq, "length:", len(dna))

#TO CHANGE THE KEY NAME

# for seq, dna in seq_dna.items():
#     seq = f"seq_{keynames.index(seq)}_{len(dna)}"
#     print(seq)

    # length = len(dna)
    # print(length)

# for key in seq_dna.keys():
#     print(key)
    # seq_1_len(len(dna))
    #
    #
    # keyname
    # print(keynames[1])

# values = list(seq_dna.values())
# print(values[0])
#
#KEY NAME ORDER
# keynames =list(seq_dna.keys())
# print(keynames[0])

## IMPORTING DATA AS FASTA##

from Bio import SeqIO

fasta_dna = "dna_sequences.fasta"

seq_dna = {}

for record in SeqIO.parse(fasta_dna, "fasta"):
    seq_dna[record.id] = str(record.seq)
print(seq_dna)

##SEQUENCES IN THIS DICTIONARY ARE ACCEPTED AS CODING STRAND.##

# keys_list = list(seq_dna.keys())
# print(keys_list)   #Pulling the key names from the dictionary (Names are incorrect)

## CORRECTING SEQUENCE NAMES ##

keynames = list(seq_dna.keys())
new_seq_dna = {}
for order, dna_strings in enumerate(keynames):
    dna = seq_dna[dna_strings]
    order = order + 1
    new_name = f"seq_{order}_len{len(dna)}"
    new_seq_dna[new_name] = dna

seq_dna = new_seq_dna
print(seq_dna)

## FREQUENCY OF NUCLEOTIDES ##

frequency_nuc = {}
for key, dna in seq_dna.items():
    freq_a = round(dna.count('A')/len(dna), 2)
    freq_t = round(dna.count('T')/len(dna), 2)
    freq_g = round(dna.count('G')/len(dna), 2)
    freq_c = round(dna.count('C')/len(dna), 2)
    GC_Content =  freq_g + freq_c
    keynames = key
    frequency_nuc[key] = f"A:{freq_a:.2f} T:{freq_t:.2f} G:{freq_g:.2f} C:{freq_c:.2f} GC:{GC_Content:.2f}"


print(frequency_nuc)

## CODING STRAND → RNA (T → U) ##

seq_rna = {}
for key, dna in seq_dna.items():
    rna = dna.replace('T', 'U')
    seq_rna[key] = rna

print(seq_rna)


## MAPPING EVERY 3 BASES: RNA → CODONS ##

seq_codon = {}

for key, rna in seq_rna.items():
    codons = []
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        codons.append(codon)
    seq_codon[key] = codons
print(seq_codon)

## CODON → AMINO ACIDS ##

from Bio.Data import CodonTable
table = CodonTable.unambiguous_rna_by_name["Standard"] #Pulling Codon Table into the project
#print(table)

seq_aa = {}
for seq_id, codon_list in seq_codon.items():
    aa_list = []

    # AA CONVERSION CODE 1
    for codon in codon_list:
        if len(codon) != 3:
            continue
        if codon in table.stop_codons:
            break
        aa = table.forward_table.get(codon, "X")
        aa_list.append(aa)

    seq_aa[seq_id] = "".join(aa_list)

print(seq_aa)

### AA CONVERSION CODE 2 ##
#     for codon in codon_list:
#         if len(codon) == 3:
#             aa = table.forward_table.get(codon, "X")
#         if codon in table.stop_codons:
#             break
#         aa_list.append(aa)
#     seq_aa[seq_id] = "".join(aa_list)
#
# print(seq_aa)

### AA CONVERSION CODE 3 ##
#     for codon in codon_list:
#         if len(codon) == 3:
#             codon == aa_list
#         if codon in table.stop_codons:
#             break
#         aa = table.forward_table.get(codon, "X")
#         aa_list.append(aa)
#     seq_aa[seq_id] = "".join(aa_list)
#
# print(seq_aa)


## CREATING TABLE FOR OVERALL DATA WITH PANDAS

import pandas as pd

pd.set_option("display.max_colwidth", None)
pd.set_option("display.width", None)

overall_data = {'ID' : list(seq_dna.keys()),
                'Coding Strand' : list(seq_dna.values()),
                'Nucleotide Frequency' : list(frequency_nuc.values()),
                'Codon Sequences' : list(seq_rna.values()),
                'Amino Acid Sequences' : list(seq_aa.values())
}
df = pd.DataFrame(overall_data)
print(df)
# styled = df.style.background_gradient(cmap="viridis")
# styled.to_html("table.html")

