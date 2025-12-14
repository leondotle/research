# The script assumes all the related package is installed
simulator.py genome     -rg RB1_100.fa     -c human_NA12878_DNA_FAB49712_guppy/training     -o rb1_sim     -n 1000     --fastq     -b guppy
cat rb1_sim_aligned_reads.fastq rb1_sim_unaligned_reads.fastq > patient_sample_raw.fastq
minimap2 -ax map-ont RB1_ref.fa patient_sample_raw.fastq > patient_aligned.sam
samtools view -S -b patient_aligned.sam > patient_aligned.bam
samtools sort patient_aligned.bam -o patient_aligned_sorted.bam
samtools index patient_aligned_sorted.bam
samtools view -S -b patient_aligned.sam > patient_aligned.bam
samtools sort patient_aligned.bam -o patient_aligned_sorted.bam
python flare_coverage.py 
