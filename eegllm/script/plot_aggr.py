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

# Function to average PSDs from different electrodes
def average_psd(psd_values):
    return np.mean(psd_values, axis=0)

# Sampling rate (assuming 256 Hz)
fs = 256

# Initialize buffers for storing the latest samples
window_size = 1024  # Number of samples in the sliding window (e.g., 1 second of data)
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
    
    # Apply band-pass filters to isolate delta, theta, alpha, and beta bands
    delta_tp9 = bandpass_filter(tp9_buffer, 0.5, 4, fs)
    theta_tp9 = bandpass_filter(tp9_buffer, 4, 8, fs)
    alpha_tp9 = bandpass_filter(tp9_buffer, 8, 13, fs)
    beta_tp9 = bandpass_filter(tp9_buffer, 13, 30, fs)

    delta_af7 = bandpass_filter(af7_buffer, 0.5, 4, fs)
    theta_af7 = bandpass_filter(af7_buffer, 4, 8, fs)
    alpha_af7 = bandpass_filter(af7_buffer, 8, 13, fs)
    beta_af7 = bandpass_filter(af7_buffer, 13, 30, fs)

    delta_af8 = bandpass_filter(af8_buffer, 0.5, 4, fs)
    theta_af8 = bandpass_filter(af8_buffer, 4, 8, fs)
    alpha_af8 = bandpass_filter(af8_buffer, 8, 13, fs)
    beta_af8 = bandpass_filter(af8_buffer, 13, 30, fs)

    delta_tp10 = bandpass_filter(tp10_buffer, 0.5, 4, fs)
    theta_tp10 = bandpass_filter(tp10_buffer, 4, 8, fs)
    alpha_tp10 = bandpass_filter(tp10_buffer, 8, 13, fs)
    beta_tp10 = bandpass_filter(tp10_buffer, 13, 30, fs)
    
    # Calculate PSD for each frequency band
    psd_delta_tp9 = calculate_psd(delta_tp9)
    psd_theta_tp9 = calculate_psd(theta_tp9)
    psd_alpha_tp9 = calculate_psd(alpha_tp9)
    psd_beta_tp9 = calculate_psd(beta_tp9)

    psd_delta_af7 = calculate_psd(delta_af7)
    psd_theta_af7 = calculate_psd(theta_af7)
    psd_alpha_af7 = calculate_psd(alpha_af7)
    psd_beta_af7 = calculate_psd(beta_af7)

    psd_delta_af8 = calculate_psd(delta_af8)
    psd_theta_af8 = calculate_psd(theta_af8)
    psd_alpha_af8 = calculate_psd(alpha_af8)
    psd_beta_af8 = calculate_psd(beta_af8)

    psd_delta_tp10 = calculate_psd(delta_tp10)
    psd_theta_tp10 = calculate_psd(theta_tp10)
    psd_alpha_tp10 = calculate_psd(alpha_tp10)
    psd_beta_tp10 = calculate_psd(beta_tp10)
    
    # Combine PSDs using averaging
    delta_psd_combined = average_psd([psd_delta_tp9, psd_delta_af7, psd_delta_af8, psd_delta_tp10])
    theta_psd_combined = average_psd([psd_theta_tp9, psd_theta_af7, psd_theta_af8, psd_theta_tp10])
    alpha_psd_combined = average_psd([psd_alpha_tp9, psd_alpha_af7, psd_alpha_af8, psd_alpha_tp10])
    beta_psd_combined = average_psd([psd_beta_tp9, psd_beta_af7, psd_beta_af8, psd_beta_tp10])
    
    # Compute Alpha/Theta ratios
    alpha_theta_ratio_combined = alpha_psd_combined / theta_psd_combined
    
    # Return the computed values
    return {
        'Delta Power Combined': delta_psd_combined,
        'Theta Power Combined': theta_psd_combined,
        'Alpha Power Combined': alpha_psd_combined,
        'Beta Power Combined': beta_psd_combined,
        'Alpha/Theta Ratio Combined': alpha_theta_ratio_combined
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
fig, ax = plt.subplots(5, 1, figsize=(10, 12))

# Initialize data containers for plotting
delta_power = []
theta_power = []
alpha_power = []
beta_power = []
alpha_theta_ratio = []

# Define update function for animation
def update(frame):
    row = frame
    result = process_new_data(row)
    
    # Append new data to the containers
    delta_power.append(result['Delta Power Combined'])
    theta_power.append(result['Theta Power Combined'])
    alpha_power.append(result['Alpha Power Combined'])
    beta_power.append(result['Beta Power Combined'])
    alpha_theta_ratio.append(result['Alpha/Theta Ratio Combined'])
    
    # Update plots
    ax[0].cla()
    ax[0].plot(delta_power, label='Delta Power Combined')
    ax[0].legend(loc='upper right')
    
    ax[1].cla()
    ax[1].plot(theta_power, label='Theta Power Combined')
    ax[1].legend(loc='upper right')
    
    ax[2].cla()
    ax[2].plot(alpha_power, label='Alpha Power Combined')
    ax[2].legend(loc='upper right')
    
    ax[3].cla()
    ax[3].plot(beta_power, label='Beta Power Combined')
    ax[3].legend(loc='upper right')
    
    ax[4].cla()
    ax[4].plot(alpha_theta_ratio, label='Alpha/Theta Ratio Combined')
    ax[4].legend(loc='upper right')

# Stream data and animate
ani = animation.FuncAnimation(fig, update, frames=stream_eeg_data(file_path='eeg_raw.csv'), interval=100)

plt.tight_layout()
plt.show()
