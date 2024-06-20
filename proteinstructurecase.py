#In this cell, we are installing the libraries that we will need for the first question.
!pip install requests
import requests
from collections import Counter #we defined collections and Counter library; Counter give us count of the sequence

uniprot_api = "https://www.uniprot.org/uniprot/"  #We define the UniProt website URL.
prt_accession_id = "P0DTD1"  #NSP3 non-structural protein ID in UniProt Database
response = requests.get(f"{uniprot_api}{prt_accession_id}.fasta")  
#An HTTP request is made to retrieve the FASTA-formatted sequence of the specified protein using the UniProt URL (uniprot_api) and the protein accession ID (prt_accession_id).

if response.status_code == 200:  # Checks if the HTTP request was successfully completed. If successful, continues processing the protein's data in FASTA format.
    # We create the amino acid sequence from the FASTA file.
    data = response.text.split("\n")  # Splits the data in FASTA format line by line and creates a list named 'data'.
    sequence = ''.join(data[1:])  # Concatenates the part of the data from the second line onwards in FASTA format, creating a string named 'sequence'. This contains the protein's amino acid sequence.

    # We calculate the count of amino acids.
    aa_count = Counter(sequence)  # Uses the Counter class to count the occurrences of each amino acid in the amino acid sequence.

    # We display the amino acid sequence.
    print(f"{prt_accession_id}:")
    for amino_acid, count in aa_count.items():
        print(f"{amino_acid}: {count}")  # Prints the counts of amino acids in the sequence to the screen.

else:  # If the HTTP request is unsuccessful or insufficient, prints an error message.
    print(f"{prt_accession_id}: Data couldn't be retrieved. Status Code: {response.status_code}")

#In this cell, we are installing the libraries that we will need for the second question.
!pip install biopython
from Bio import Entrez

#The NCBI Entrez accession number for the gene is assigned to the gene_id variable
gene_id = "NC_045512.2"  # or gene_id_or_name = "NSP3"  # A specific gene's NCBI Entrez accession number is assigned to the gene_id variable

#An email address is specified for the owner of requests sent to the NCBI Entrez API
Entrez.email = "gulnuruzun20@gmail.com"  # This is the email address that will be used as the owner of requests sent to the NCBI Entrez API. It is important for NCBI authorities to contact you in case of excessive usage.

# The nucleotide sequence of the gene is retrieved
handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")  # An NCBI Entrez API request is made using the Entrez.efetch function to retrieve the data of a specific nucleotide gene. The parameters are as follows:
# db: Type of database (nucleotide)
# id: NCBI accession number of the gene
# rettype: Type of returned data (in GenBank format)
# retmode: Mode of returned data (in text format)

genbank_record = handle.read()

#Only the nucleotide sequence is extracted from the GenBank format data
nucleotide_sequence = ""
reading_sequence = False  # A flag to control the start of the nucleotide sequence

for line in genbank_record.split("\n"):
    if line.startswith("ORIGIN"):
        # The nucleotide sequence begins from the lines under the ORIGIN tag
        reading_sequence = True
    elif reading_sequence and line.strip():  # Skip empty lines at the end of the nucleotide sequence
        # Clear spaces in the line and add to the nucleotide sequence
        nucleotide_sequence += "".join(line.split()[1:])

#Print the obtained nucleotide sequence to the screen
print("Nucleotide Sequence of the Gene:")
print(nucleotide_sequence)

#In this cell, we are installing the libraries that we will need for the third question.
!pip install biopython py3dmol
from Bio.PDB import PDBParser
import py3Dmol
from collections import Counter
import matplotlib.pyplot as plt
import requests

#Downloaded .pdb file.
pdb_file = "ourpath"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)

#Obtaining chain information to be used in subsequent analyses.
amino_acids = []
chain_lengths = []
atom_types = []

for model in structure:
    for chain in model:
        chain_id = chain.get_id()
        print(f"\nChain for analyses {chain_id} information:")

        #We are obtaining chain information.
        for residue in chain:
            for atom in residue:
                #We are obtaining atom information.
                atom_name = atom.get_name()
                residue_name = residue.get_resname()
                atom_coordinates = atom.get_coord()
                print(f"Chain {chain_id}, Residue {residue_name}, Atom {atom_name}: {atom_coordinates}")

                #We are creating the list of amino acids.
                amino_acids.append(residue_name)

        # We are calculate the chain length.
        chain_length = len(list(chain))
        print(f"Zincir {chain_id} Uzunluğu: {chain_length}")
        chain_lengths.append(chain_length)

        #We are add atom types to the list
        for residue in chain:
            for atom in residue:
                atom_types.append(atom.element)

#Visualization stage of analysis results.

#Amino acid composition analysis (Bar plot)
amino_acid_counts = Counter(amino_acids)
plt.figure(figsize=(8, 6))
plt.bar(amino_acid_counts.keys(), amino_acid_counts.values(), color='skyblue')
plt.xlabel('Amino Acid')
plt.ylabel('Count')
plt.title('Amino Acid Composition')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

#Chain length analysis (Histogram)
plt.figure(figsize=(8, 6))
plt.hist(chain_lengths, bins=20, edgecolor='black', color='salmon')
plt.xlabel('Chain Length')
plt.ylabel('Frequency')
plt.title('Distribution of Chain Lengths')
plt.tight_layout()
plt.show()

#Atom type analysis (Bar plot)
plt.figure(figsize=(8, 6))
atom_type_counts = Counter(atom_types)
plt.bar(atom_type_counts.keys(), atom_type_counts.values(), color='lightgreen')
plt.xlabel('Atom Type')
plt.ylabel('Count')
plt.title('Atom Type Composition')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

#Create a Py3Dmol viewer
viewer = py3Dmol.view(width=400, height=400)

#Load the PDB structure and add it to the viewer
pdb_data = open(pdb_file, "r").read()
viewer.addModel(pdb_data, "pdb")

#Show the 3D visualization
viewer.zoomTo()
viewer.show()

#Additionally, visualization to protein structure using different library like os and nglview.
#Visualization with NGLViews Library (optional)
!pip install nglview
import os
import nglview as nv

filepath = os.path.join("PDB_file", "6W9C.pdb")
view = nv.show_file("ourpath") #my path
view

from Bio import PDB
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser

# Downloaded .pdb file.
pdb_file = "C:/Users/gulnu/OneDrive/Desktop/BSB511-HW1/6W9C.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)

# Let's assume there is a variant at a specific amino acid position. UniProt shows 682 different variants.
# We will select 20 variants and visualize them by specifying their positions.
# Specified variant positions
variant_positions = [10, 121, 179, 272, 282, 831, 1005, 1064, 1115, 1364, 1440, 2040, 2347, 2387, 2691, 2925, 3162, 3390, 3746, 4001]

# Create a 3D plot for visualization
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot each atom of the protein
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                ax.scatter(atom.coord[0], atom.coord[1], atom.coord[2], c='b', marker='o')

# Mark specified variant positions
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[1] in variant_positions:
                for atom in residue:
                    ax.scatter(atom.coord[0], atom.coord[1], atom.coord[2], c='r', marker='^', s=100)

# Axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.title('Protein Structure - Specified Variant Positions')
plt.show()