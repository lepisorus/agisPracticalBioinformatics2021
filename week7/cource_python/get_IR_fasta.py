from Bio import SeqIO

sequence_record = SeqIO.read("sequence.fasta", "fasta")
with open("blastn.txt", "rt") as blastn_in:
    lines = blastn_in.readlines()
    count = 0
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        length = int(fields[3])
        qs = int(fields[6])
        qe = int(fields[7])
        if fields[8] > fields[9]:
            ss = int(fields[9])
            se = int(fields[8])
        else:
            ss = int(fields[8])
            se = int(fields[9])
        if length == len(sequence_record.seq):
            continue
        elif length < 100:
            continue
        if qs < ss:
            start_1 = qs
            end_1 = qe
            start_2 = ss
            end_2 = se
        else:
            start_1 = ss
            end_1 = se
            start_2 = qs
            end_2 = qe
        count += 1
        repeat_1_id = sequence_record.id + "_" + str(start_1) + "_" + str(end_1)
        repeat_1_seq = sequence_record.seq[(start_1-1):end_1]
        repeat_2_id = sequence_record.id + "_" + str(start_1) + "_" + str(end_1)
        repeat_2_seq = sequence_record.seq[(start_2-1):end_2]
        with open("repeat_" + str(count) + "_1_" + str(start_1) + "_" + str(end_1) + ".fasta", "wt") as fasta_1_out:
            print(">" + repeat_1_id, file=fasta_1_out)
            print(repeat_1_seq, file=fasta_1_out)
        with open("repeat_" + str(count) + "_2_" + str(start_2) + "_" + str(end_2) + ".fasta", "wt") as fasta_2_out:
            print(">" + repeat_2_id, file=fasta_2_out)
            print(repeat_2_seq, file=fasta_2_out)

