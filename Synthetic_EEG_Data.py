import numpy as np
import json
import os

def generate_cleaner_eeg_data(num_entries):
    data = []
    for _ in range(num_entries):
        theta_power = np.round(np.random.uniform(0.1, 0.5), 5)
        alpha_power = np.round(np.random.uniform(0.1, 0.4), 5)
        beta_power = np.round(np.random.uniform(0.5, 2.0), 5)
        alpha_over_theta = np.round(alpha_power / theta_power, 5)

        # Cleaner mental state logic without explicit thresholds in notes
        if alpha_over_theta > 1.0:
            mental_state = "relaxed"
            note = "High Alpha over Theta indicates relaxation."
        elif alpha_over_theta < 0.6:
            if beta_power > 1.2:
                mental_state = "concentrating"
                note = "Low Alpha over Theta and high Beta Power indicate concentration."
            else:
                mental_state = "drowsy"
                note = "Low Alpha over Theta and low Beta Power indicate drowsiness."
        else:
            if beta_power > 1.2:
                mental_state = "concentrating"
                note = "Balanced Alpha over Theta and high Beta Power indicate concentration."
            else:
                mental_state = "neutral"
                note = "Balanced EEG patterns indicate a neutral state."

        data.append({
            "Theta Power": theta_power,
            "Alpha Power": alpha_power,
            "Beta Power": beta_power,
            "Alpha_over_Theta": alpha_over_theta,
            "Mental State": mental_state,
            "Note": note
        })
    return data

# Generate the dataset with 10,000 entries
cleaner_synthetic_eeg_data = generate_cleaner_eeg_data(10000)

# Ensure the directory exists
output_directory = "/mnt/data/"
os.makedirs(output_directory, exist_ok=True)

# Save the dataset as a JSON file
cleaner_file_path = os.path.join(output_directory, "Cleaner_Synthetic_EEG_Data.json")
with open(cleaner_file_path, "w") as f:
    json.dump(cleaner_synthetic_eeg_data, f, indent=4)

print(f"Dataset successfully generated and saved to: {cleaner_file_path}")
