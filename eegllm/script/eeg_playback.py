import time
import csv
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

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
window_size = 256  # Number of samples in the sliding window (e.g., 1 second of data)
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
        'Timestamp': new_data['Timestamp'],
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

def mood_to_color(delta, theta, alpha, beta):
    """
    Map the mood values to a color.

    :param delta: Delta power value
    :param theta: Theta power value
    :param alpha: Alpha power value
    :param beta: Beta power value
    :return: Hex color code
    """
    # Example logic to map mood to color (you can adjust this mapping)
    mood_value = (alpha + beta) / (alpha + beta + theta + 1e-6)  # Avoid division by zero
    # Normalize mood_value to range between 0 and 1
    normalized_mood = max(0.0, min(1.0, mood_value))

    # Map normalized mood to a color gradient from blue to red
    r = int(normalized_mood * 255)
    g = int((1 - normalized_mood) * 255)
    b = 0

    print(f"alpha: {alpha}, beta: {beta}, delta: {delta}, theta: {theta}")
    print(f"mood_value: {mood_value}, normalized_mood: {normalized_mood}, r: {r}, g: {g}, b {b}")

    return f'0x{r:02x}{g:02x}{b:02x}'

# File path to save the results
output_file_path = 'eeg_results.csv'

# Open the CSV file and write the header
with open(output_file_path, mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=[
        'Timestamp', 'Delta Power Combined', 'Theta Power Combined', 
        'Alpha Power Combined', 'Beta Power Combined', 'Alpha/Theta Ratio Combined'
    ])
    writer.writeheader()
    
    # Initialize the Chrome driver with service and options
    service = Service(ChromeDriverManager().install())
    options = webdriver.ChromeOptions()
    driver = webdriver.Chrome(service=service, options=options)

    # Open the webpage initially
    driver.get("https://leondotle.github.io/origami/webvr.html")

    # Simulate streaming data processing
    for row in stream_eeg_data(file_path='eeg_raw.csv'):
        result = process_new_data(row)
        writer.writerow(result)

        # Compute mood and map to color
        color = mood_to_color(result['Delta Power Combined'], result['Theta Power Combined'], 
                              result['Alpha Power Combined'], result['Beta Power Combined'])

        # Update the cube color dynamically using JavaScript
        script = f"updateCubeColor({color});"
        driver.execute_script(script)
        print(f"Updated cube color to: {color}")
        
        time.sleep(3)  # Sleep for 3 seconds to see the new color

    # Close the browser after processing
    driver.quit()
