import pysam
import csv
from collections import defaultdict

# Inputs
bam_in_path = "/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/subclustering/subclustering_FAP_02_2025_bigwigs/results/temp-bams-by-cluster-clamp/atac-fibrous.bam"
barcode_to_donor_path = "../barcode_assignments.tsv"
output_prefix = "fibrous_"  # e.g. donor_001.bam

# Read barcode-to-donor map
barcode_to_donor = {}
with open(barcode_to_donor_path) as f:
    reader = csv.reader(f, delimiter='\t')
    for barcode, donor in reader:
        # Barcodes may be stored with "-1" or not. Adjust as needed.
        barcode_to_donor[barcode.strip()] = donor.strip()

# Reverse map: donors → list of barcodes
donor_to_barcodes = defaultdict(set)
for bc, donor in barcode_to_donor.items():
    donor_to_barcodes[donor].add(bc)

# Create BAM writers per donor
bam_in = pysam.AlignmentFile(bam_in_path, "rb")
bam_out_files = {
    donor: pysam.AlignmentFile(f"{output_prefix}{donor}.bam", "wb", template=bam_in)
    for donor in donor_to_barcodes
}

# Iterate through BAM and write to the correct donor BAM file
for read in bam_in.fetch(until_eof=True):
    if read.has_tag("CB"):
        barcode = read.get_tag("CB")
        donor = barcode_to_donor.get(barcode)
        if donor and donor in bam_out_files:
            bam_out_files[donor].write(read)

# Close all files
bam_in.close()
for f in bam_out_files.values():
    f.close()

# Indexing
for donor in bam_out_files:
    pysam.index(f"{output_prefix}{donor}.bam")

