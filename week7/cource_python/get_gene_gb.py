from Bio import SeqIO
records = list(SeqIO.parse("mito.gb", "genbank"))
index = list()
for i in range(len(records[0].features)): 
    if records[0].features[i].type != "gene":
        index.append(i)
index.reverse()
for item in index:
    records[0].features.pop(item)
SeqIO.write(records, "mito_gene.gb", "genbank")
