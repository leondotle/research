import pysam
import matplotlib.pyplot as plt
import numpy as np

# 1. Configuration & Load Alignment File
bam_file = "patient_aligned_sorted.bam"
# Ensure the file exists or replace with your actual path
try:
    samfile = pysam.AlignmentFile(bam_file, "rb")
except FileNotFoundError:
    print(f"Error: Could not find '{bam_file}'. Please check the path.")
    exit()

# Target Deletion Coordinates
del_start = 49951
del_end = 54950

# Analysis Boundaries (Trimmed "Global" Range)
# We trim the start and end of the gene to avoid edge artifacts.
analysis_start = 20000
analysis_end = 160000

# 2. Calculate Coverage across the gene
reference_length = 181703 
coverages = np.zeros(reference_length)

print("Calculating coverage... this may take a moment.")

# .pileup() iterates through every base in the reference
for pileupcolumn in samfile.pileup():
    pos = pileupcolumn.reference_pos
    if pos < reference_length:
        coverages[pos] = pileupcolumn.nsegments

# 3. Basic Stats & Region Extraction
# Extract the specific deletion region and the background
# Background is now defined as the trimmed analysis range EXCLUDING the deletion
left_background = coverages[analysis_start : del_start]
right_background = coverages[del_end : analysis_end]
background = np.concatenate([left_background, right_background])
deletion_region = coverages[del_start : del_end]

# Basic Calculations
mean_cov = np.mean(coverages) # Global mean (untouched)
mu_bg = np.mean(background)   # Background mean (Trimmed Range)
mu_del = np.mean(deletion_region) # Deletion mean

print(f"--- FLARE Report ---")
print(f"Global Mean Depth: {mean_cov:.2f}x")

print(f"--- Targeted Region [{del_start}-{del_end}] ---")
print(f"Mean Depth (Background 20k-160k): {mu_bg:.2f}x")
print(f"Mean Depth (Deletion):            {mu_del:.2f}x")

# --- Original Metric ---
C_out = mean_cov # Using global mean as reference
C_del = mu_del
dip_score = 1.0 - (C_del / C_out) if C_out > 0 else 0.0

# --- Advanced Metrics ---
print(f"\n--- Deletion Evidence Metrics ---")
print(f"1. Standard Dip Score: {dip_score:.4f} (Relative drop)")

# Metric 2: Z-Score (Signal-to-Noise Ratio)
# Measures how many standard deviations the drop is away from the background mean.
# High Z-score (>3) means the drop is significant relative to background noise.
sigma_bg = np.std(background)
z_score = (mu_bg - mu_del) / (sigma_bg + 1e-9) # 1e-9 avoids division by zero
print(f"2. Z-Score:            {z_score:.4f} (Sigmas dropped below background)")

# Metric 3: Cohen's d (Effect Size)
# Accounts for noise inside the deletion as well.
sigma_del = np.std(deletion_region)
n1, n2 = len(background), len(deletion_region)
# Pooled Standard Deviation calculation
s_pooled = np.sqrt(((n1 - 1) * sigma_bg**2 + (n2 - 1) * sigma_del**2) / (n1 + n2 - 2))
cohens_d = (mu_bg - mu_del) / (s_pooled + 1e-9)
print(f"3. Cohen's d:          {cohens_d:.4f} (Separation between distributions)")

# Metric 4: Integrated Deficit
# Sum of missing reads across the whole deletion width. 
# Good for detecting wide, shallow deletions.
deficits = mu_bg - deletion_region
# We only sum positive deficits (where deletion is lower than background)
integrated_deficit = np.sum(deficits[deficits > 0])
print(f"4. Integrated Deficit: {integrated_deficit:.0f} (Total missing coverage area)")


# 4. Visualization
plt.figure(figsize=(12, 6))

# Plot full coverage
plt.plot(coverages, color='lightgray', alpha=0.7, label='Raw Coverage')

# Highlight Background Region (Trimmed Range)
x_left = np.arange(analysis_start, del_start)
x_right = np.arange(del_end, analysis_end)
plt.plot(x_left, left_background, color='green', alpha=0.6, label='Background (20k-160k)')
plt.plot(x_right, right_background, color='green', alpha=0.6)

# Highlight Deletion Region
x_del = np.arange(del_start, del_end)
plt.plot(x_del, deletion_region, color='red', alpha=0.8, label='Deletion Region')

# Add Mean lines for visual comparison
plt.axhline(y=mu_bg, color='green', linestyle='--', alpha=0.8, label=f'Background Mean ({mu_bg:.1f}x)')
plt.axhline(y=mu_del, color='red', linestyle='--', alpha=0.8, label=f'Deletion Mean ({mu_del:.1f}x)')

plt.title(f"FLARE Assay: Deletion Analysis (Z-Score: {z_score:.2f})")
plt.xlabel("Genomic Position (bp)")
plt.ylabel("Read Depth")
plt.legend(loc='upper right')

# Add text box with metrics
stats_text = (
    f"Dip Score: {dip_score:.2f}\n"
    f"Z-Score: {z_score:.2f}\n"
    f"Cohen's d: {cohens_d:.2f}"
)
plt.gca().text(0.02, 0.95, stats_text, transform=plt.gca().transAxes, 
               fontsize=10, verticalalignment='top', 
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig("rb1_coverage_advanced_plot.png")
print("\nCoverage plot saved to 'rb1_coverage_advanced_plot.png'")

samfile.close()
