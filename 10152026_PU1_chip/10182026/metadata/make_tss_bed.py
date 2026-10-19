import gzip
import re

gtf = "gencode.v50.annotation.gtf.gz"
tss_records = set()

with gzip.open(gtf, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "transcript":
            continue

        attrs = dict(re.findall(r'(\S+) "([^"]+)"', fields[8]))
        if (attrs.get("gene_type") != "protein_coding"
                or attrs.get("transcript_type") != "protein_coding"):
            continue

        chrom, strand = fields[0], fields[6]
        start, end = int(fields[3]), int(fields[4])
        gene_id = attrs.get("gene_id", "")
        gene_name = attrs.get("gene_name", gene_id)

        # Convert the 1-based GTF TSS to a 0-based BED coordinate
        tss0 = start - 1 if strand == "+" else end - 1
        tss_records.add((chrom, tss0, gene_id, gene_name, strand))

def chrom_key(chrom):
    name = chrom.removeprefix("chr")
    return (0, int(name)) if name.isdigit() else (1, name)

records = sorted(tss_records, key=lambda x: (chrom_key(x[0]), x[1], x[2]))

with open("protein_coding_TSS_hg38.bed", "w") as tss_out, \
     open("protein_coding_TSS_2kb_windows_hg38.bed", "w") as win_out:
    for chrom, tss0, gene_id, gene_name, strand in records:
        label = f"{gene_name}|{gene_id}"
        tss_out.write(f"{chrom}\t{tss0}\t{tss0 + 1}\t{label}\t0\t{strand}\n")

        left = max(0, tss0 - 2000)
        right = tss0 + 2000
        win_out.write(f"{chrom}\t{left}\t{right}\t{label}\t0\t{strand}\n")
