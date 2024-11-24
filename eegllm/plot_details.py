import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

# Define band-pass filter function
def bandpass_filter(data, lowcut, highcut, fs, order=2):
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

# Function to calculate power spectral density (PSD)
def calculate_psd(signal):
    return np.mean(signal ** 2)

# Sampling rate (assuming 256 Hz)
fs = 256

# Initialize buffers for storing the latest samples
window_size = 128  # Number of samples in the sliding window (e.g., 1 second of data)
tp9_buffer = np.zeros(window_size)
af7_buffer = np.zeros(window_size)
af8_buffer = np.zeros(window_size)
tp10_buffer = np.zeros(window_size)

def process_new_data(new_data):
    global tp9_buffer, af7_buffer, af8_buffer, tp10_buffer
    
    # Update buffers with new data
    tp9_buffer = np.roll(tp9_buffer, -1)
    tp9_buffer[-1] = new_data['TP9']
    
    af7_buffer = np.roll(af7_buffer, -1)
    af7_buffer[-1] = new_data['AF7']
    
    af8_buffer = np.roll(af8_buffer, -1)
    af8_buffer[-1] = new_data['AF8']
    
    tp10_buffer = np.roll(tp10_buffer, -1)
    tp10_buffer[-1] = new_data['TP10']
    
    # Apply band-pass filters to isolate theta, alpha, and beta bands
    theta_tp9 = bandpass_filter(tp9_buffer, 4, 8, fs)
    alpha_tp9 = bandpass_filter(tp9_buffer, 8, 13, fs)
    beta_tp9 = bandpass_filter(tp9_buffer, 13, 30, fs)

    theta_af7 = bandpass_filter(af7_buffer, 4, 8, fs)
    alpha_af7 = bandpass_filter(af7_buffer, 8, 13, fs)
    beta_af7 = bandpass_filter(af7_buffer, 13, 30, fs)

    theta_af8 = bandpass_filter(af8_buffer, 4, 8, fs)
    alpha_af8 = bandpass_filter(af8_buffer, 8, 13, fs)
    beta_af8 = bandpass_filter(af8_buffer, 13, 30, fs)

    theta_tp10 = bandpass_filter(tp10_buffer, 4, 8, fs)
    alpha_tp10 = bandpass_filter(tp10_buffer, 8, 13, fs)
    beta_tp10 = bandpass_filter(tp10_buffer, 13, 30, fs)
    
    # Calculate PSD for each frequency band
    psd_theta_tp9 = calculate_psd(theta_tp9)
    psd_alpha_tp9 = calculate_psd(alpha_tp9)
    psd_beta_tp9 = calculate_psd(beta_tp9)

    psd_theta_af7 = calculate_psd(theta_af7)
    psd_alpha_af7 = calculate_psd(alpha_af7)
    psd_beta_af7 = calculate_psd(beta_af7)

    psd_theta_af8 = calculate_psd(theta_af8)
    psd_alpha_af8 = calculate_psd(alpha_af8)
    psd_beta_af8 = calculate_psd(beta_af8)

    psd_theta_tp10 = calculate_psd(theta_tp10)
    psd_alpha_tp10 = calculate_psd(alpha_tp10)
    psd_beta_tp10 = calculate_psd(beta_tp10)
    
    # Compute Alpha/Theta ratios
    alpha_theta_ratio_tp9 = psd_alpha_tp9 / psd_theta_tp9
    alpha_theta_ratio_af7 = psd_alpha_af7 / psd_theta_af7
    alpha_theta_ratio_af8 = psd_alpha_af8 / psd_theta_af8
    alpha_theta_ratio_tp10 = psd_alpha_tp10 / psd_theta_tp10
    
    # Return the computed values
    return {
        'Theta TP9': psd_theta_tp9,
        'Alpha TP9': psd_alpha_tp9,
        'Beta TP9': psd_beta_tp9,
        'Alpha/Theta Ratio TP9': alpha_theta_ratio_tp9,
        'Theta AF7': psd_theta_af7,
        'Alpha AF7': psd_alpha_af7,
        'Beta AF7': psd_beta_af7,
        'Alpha/Theta Ratio AF7': alpha_theta_ratio_af7,
        'Theta AF8': psd_theta_af8,
        'Alpha AF8': psd_alpha_af8,
        'Beta AF8': psd_beta_af8,
        'Alpha/Theta Ratio AF8': alpha_theta_ratio_af8,
        'Theta TP10': psd_theta_tp10,
        'Alpha TP10': psd_alpha_tp10,
        'Beta TP10': psd_beta_tp10,
        'Alpha/Theta Ratio TP10': alpha_theta_ratio_tp10
    }

def stream_eeg_data(file_path, chunk_size=1):
    """
    Generator function to simulate streaming EEG data from a CSV file.

    :param file_path: Path to the CSV file.
    :param chunk_size: Number of rows to yield at a time.
    :yield: Chunk of data as a DataFrame.
    """
    for chunk in pd.read_csv(file_path, chunksize=chunk_size):
        for _, row in chunk.iterrows():
            yield row.to_dict()

# Initialize plot
fig, ax = plt.subplots(4, 4, figsize=(20, 16))
ax = ax.flatten()

# Initialize data containers for plotting
psd_data = {
    'Theta TP9': [],
    'Alpha TP9': [],
    'Beta TP9': [],
    'Alpha/Theta Ratio TP9': [],
    'Theta AF7': [],
    'Alpha AF7': [],
    'Beta AF7': [],
    'Alpha/Theta Ratio AF7': [],
    'Theta AF8': [],
    'Alpha AF8': [],
    'Beta AF8': [],
    'Alpha/Theta Ratio AF8': [],
    'Theta TP10': [],
    'Alpha TP10': [],
    'Beta TP10': [],
    'Alpha/Theta Ratio TP10': []
}

# Define update function for animation
def update(frame):
    row = frame
    result = process_new_data(row)
    
    # Append new data to the containers
    for key in psd_data.keys():
        psd_data[key].append(result[key])
        # Limit the data to the latest 1024 points
        if len(psd_data[key]) > window_size:
            psd_data[key].pop(0)
    
    # Update plots - arrange all metrics for a specific location in the same row
    location_keys = [['Theta TP9', 'Alpha TP9', 'Beta TP9', 'Alpha/Theta Ratio TP9'],
                     ['Theta AF7', 'Alpha AF7', 'Beta AF7', 'Alpha/Theta Ratio AF7'],
                     ['Theta AF8', 'Alpha AF8', 'Beta AF8', 'Alpha/Theta Ratio AF8'],
                     ['Theta TP10', 'Alpha TP10', 'Beta TP10', 'Alpha/Theta Ratio TP10']]
    
    for row_idx, keys in enumerate(location_keys):
        for col_idx, key in enumerate(keys):
            ax[row_idx * 4 + col_idx].cla()
            ax[row_idx * 4 + col_idx].plot(psd_data[key], label=key)
            ax[row_idx * 4 + col_idx].legend(loc='upper right')

# Stream data and animate
ani = animation.FuncAnimation(fig, update, frames=stream_eeg_data(file_path='eeg_raw.csv'), interval=100)

plt.tight_layout()
plt.show()
