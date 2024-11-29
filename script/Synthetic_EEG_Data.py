import numpy as np
import json
import os

def generate_synthetic_eeg_data(num_entries, output_filename):
    data = []
    for _ in range(num_entries):
        # Generate random EEG power values within specified ranges
        theta_power = np.round(np.random.uniform(0.1, 0.5), 5)
        alpha_power = np.round(np.random.uniform(0.1, 0.4), 5)
        beta_power = np.round(np.random.uniform(0.5, 2.0), 5)

        # Calculate the ratio of Alpha to Theta power
        alpha_over_theta = np.round(alpha_power / theta_power, 5)

        # Determine mental state based on Alpha/Theta ratio and Beta power
        if alpha_over_theta > 1.0:
            mental_state = "relaxed"
            note = "High Alpha/Theta ratio (greater than 1.0) suggests a relaxed state."
        elif alpha_over_theta < 0.6:
            if beta_power > 1.2:
                mental_state = "concentrating"
                note = "Low Alpha/Theta ratio (less than 0.6) with high Beta power (greater than 1.2) indicates concentration."
            else:
                mental_state = "drowsy"
                note = "Low Alpha/Theta ratio (less than 0.6) with low Beta power (less than or equal to 1.2) suggests drowsiness."
        else:
            if beta_power > 1.2:
                mental_state = "concentrating"
                note = "Moderate Alpha/Theta ratio (0.6 – 1.0) with high Beta power (greater than 1.2) indicates concentration."
            else:
                mental_state = "neutral"
                note = "Moderate Alpha/Theta ratio (0.6 – 1.0) with low Beta power (less than or equal to 1.2) suggests a neutral state."

        # Append the generated data entry to the list
        data.append({
            "Theta Power": theta_power,
            "Alpha Power": alpha_power,
            "Beta Power": beta_power,
            "Alpha/Theta Ratio": alpha_over_theta,
            "Mental State": mental_state,
            "Note": note
        })

    # Define the output file path
    output_file = os.path.join("/mnt/data", output_filename)

    # Write the generated data to a JSON file
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=4)

    print(f"Synthetic EEG data successfully written to {output_file}")

# Generate another 10,000 synthetic EEG data entries for evaluation
generate_synthetic_eeg_data(10000, "synthetic_eeg_data_10k_eval.json")
