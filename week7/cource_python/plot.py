import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO
import numpy as np
import sys

record = SeqIO.read("sequence.gb", "genbank")
index = list()
for i in range(len(record.features)):
    if record.features[i].type != sys.argv[1]:
        index.append(i)
index.reverse()
for item in index:
    record.features.pop(item)

fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(15, 3), sharex=True, gridspec_kw={"height_ratios": [10, 1]}
)
# PLOT THE RECORD MAP
graphic_record = BiopythonTranslator().translate_record(record)
graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4, annotate_inline=True)
# PLOT THE LOCAL GC CONTENT (we use 50bp windows)
gc = lambda s: 100.0 * len([c for c in s if c in "GC"]) / 50
xx = np.arange(len(record.seq) - 50)
yy = [gc(record.seq[x : x + 50]) for x in xx]
ax2.fill_between(xx + 25, yy, alpha=0.3)
ax2.set_ylim(bottom=0)
ax2.set_ylabel("GC(%)")
plt.savefig('anno_' + sys.argv[1] + '.pdf')
plt.savefig('anno_' + sys.argv[1] + '.png')
