# Create batches for calculating depth from all contigs

size_map = {}
with open("BGI_contigs_100kb.txt") as contig_ref:
    contig_ref.readline()
    for line in contig_ref:
        line_split = line.split()
        
        contig = line_split[0]
        size = int(line_split[1])
        
        size_map[contig] = size

with open("candidate_contigs.txt") as candidate_file:
    candidate_contigs = candidate_file.read().splitlines()

cum_sum = {x: 0 for x in range(20)}
batch_dict = {x: list() for x in range(20)}

for contig in candidate_contigs:
    idx = min(cum_sum, key = cum_sum.get)
    cum_sum[idx] += size_map[contig]
    batch_dict[idx].append(contig)

for i in range(20):
    with open("depth_batch" + str(i), 'w') as out_file:
        for contig in batch_dict[i]:
            out_file.write(contig + '\n')
