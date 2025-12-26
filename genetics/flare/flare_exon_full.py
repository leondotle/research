import pysam
import matplotlib.pyplot as plt
import numpy as np
import os

# 1. Configuration: All 27 Exons of RB1 (hg38 relative coordinates)
# These are approximate regions to highlight the "towers" in your plot.
RB1_EXONS_ALL = [
    (1, 1, 1500), (2, 7000, 7200), (3, 25000, 25300), (4, 33000, 33200),
    (5, 47000, 47200), (6, 48500, 48700), (7, 49500, 49700), (8, 50200, 50400),
    (9, 52000, 52300), (10, 56000, 56200), (11, 59500, 59800), (12, 63000, 63300),
    (13, 67000, 67400), (14, 71500, 71800), (15, 75500, 75800), (16, 78000, 78300),
    (17, 84000, 84300), (18, 88500, 88800), (19, 92000, 92400), (20, 105000, 105400),
    (21, 114000, 114300), (22, 118000, 118400), (23, 127000, 127300), (24, 137000, 137400),
    (25, 147000, 147300), (26, 154000, 154400), (27, 175000, 178140)
]

bam_files = ["patient_aligned_sorted_1.bam", "patient_aligned_sorted_2.bam"]
existing_bams = [f for f in bam_files if os.path.exists(f)]

if not existing_bams:
    print("Error: No BAM files found!")
    exit()

# Auto-index and metadata collection
for bam in existing_bams:
    if not os.path.exists(bam + ".bai"):
        pysam.index(bam)

with pysam.AlignmentFile(existing_bams[0], "rb") as tmp:
    ref_name = tmp.references[0]
    reference_length = tmp.lengths[0]

total_coverages = np.zeros(reference_length)

# 2. Accumulate Coverage
print(f"Analyzing all 27 exons for {ref_name}...")
for bam in existing_bams:
    with pysam.AlignmentFile(bam, "rb") as samfile:
        for pileupcolumn in samfile.pileup(contig=ref_name, truncate=True):
            if pileupcolumn.reference_pos < reference_length:
                total_coverages[pileupcolumn.reference_pos] += pileupcolumn.nsegments

# 3. Comprehensive Exon Report
print(f"\n--- Full Exon Coverage Analysis ---")
missing_exons = []
for exon_num, start, end in RB1_EXONS_ALL:
    exon_data = total_coverages[start:end]
    avg_exon_depth = np.mean(exon_data) if len(exon_data) > 0 else 0
    
    status = "PASS" if avg_exon_depth > 20 else "FAIL/MISSING"
    if status == "FAIL/MISSING":
        missing_exons.append(exon_num)
    
    print(f"Exon {exon_num:02d}: {avg_exon_depth:9.2f}x | {status}")

# 4. Plotting
plt.figure(figsize=(15, 7))
plt.fill_between(range(reference_length), total_coverages, color='seagreen', alpha=0.3)
for _, start, end in RB1_EXONS_ALL:
    plt.axvspan(start, end, color='orange', alpha=0.25)

plt.title(f"RB1 Complete Exon Coverage Map ({ref_name})")
plt.xlabel("Genomic Position")
plt.ylabel("Depth")
plt.tight_layout()
plt.savefig("rb1_27_exon_analysis.png")

print(f"\nAnalysis complete. New plot saved as 'rb1_27_exon_analysis.png'.")
if missing_exons:
    print(f"Attention: {len(missing_exons)} exons showed zero/low coverage.")
