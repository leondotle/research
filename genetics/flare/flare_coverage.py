import pysam
import matplotlib.pyplot as plt
import numpy as np

# 1. Load the Alignment File
bam_file = "patient_aligned_sorted.bam"
samfile = pysam.AlignmentFile(bam_file, "rb")

# 2. Calculate Coverage across the gene
# Since RB1_ref.fa is small (~180kb), we can map position-by-position
reference_length = 181703 # Based on your wc -c output
coverages = np.zeros(reference_length)

print("Calculating coverage... this may take a moment.")

# .pileup() iterates through every base in the reference
for pileupcolumn in samfile.pileup():
    pos = pileupcolumn.reference_pos
    if pos < reference_length:
        coverages[pos] = pileupcolumn.nsegments

# 3. Basic Stats for Diagnosis
mean_cov = np.mean(coverages)
min_cov = np.min(coverages)
print(f"--- FLARE Report ---")
print(f"Mean Depth: {mean_cov:.2f}x")
print(f"Min Depth:  {min_cov}")

regional_coverages = coverages[49951:54950]
regional_mean_cov = np.mean(regional_coverages)
regional_min_cov = np.min(regional_coverages)
print(f"--- Targeted Region [49951-54950] ---")
print(f"Mean Depth (Region): {regional_mean_cov:.2f}x")
print(f"Min Depth (Region):  {regional_min_cov}")

C_out = mean_cov
C_del = regional_mean_cov
if C_out > 0:
    # Dip Score = 1 - (C_del / C_out)
    dip_score = 1.0 - (C_del / C_out)
else:
    # Avoid division by zero if overall coverage is zero
    dip_score = 0.0
print(f"--- Deletion Evidence ---")
print(f"Dip Score: {dip_score:.4f}")

# 4. Visualization (To show in your project report)
plt.figure(figsize=(10, 4))
plt.plot(coverages, color='blue', alpha=0.6)
plt.title("FLARE Assay: RB1 Gene Coverage")
plt.xlabel("Genomic Position (bp)")
plt.ylabel("Read Depth")
plt.axhline(y=30, color='r', linestyle='--', label='Clinical Threshold (30x)') # Hypothetical threshold
plt.legend()
plt.savefig("rb1_coverage_plot.png")
print("Coverage plot saved to 'rb1_coverage_plot.png'")

samfile.close()
