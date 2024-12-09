import pandas as pd
import random
import json

# Load the data from selected.csv
df = pd.read_csv('oasis_select/selected.csv')

# Define the mapping of valence_mean and arousal_mean to 9 regions
def get_suffix(valence_mean, arousal_mean):
    if arousal_mean < 3.5:
        if valence_mean < 4.0:
            return "The picture has a negative valence and low arousal."
        elif valence_mean < 4.5:
            return "The picture has a neutral valence and low arousal."
        else:
            return "The picture has a positive valence and low arousal."
    elif arousal_mean < 4.2:
        if valence_mean < 4.0:
            return "The picture has a negative valence and mild arousal."
        elif valence_mean < 4.5:
            return "The picture has a neutral valence and mild arousal."
        else:
            return "The picture has a positive valence and mild arousal."
    else:
        if valence_mean < 4.0:
            return "The picture has a negative valence and high arousal."
        elif valence_mean < 4.5:
            return "The picture has a neutral valence and high arousal."
        else:
            return "The picture has a positive valence and high arousal."

# Randomly sample indices for training and validation sets
train_indices = random.sample(range(len(df)), 90)
validation_indices = [index for index in range(len(df)) if index not in train_indices]

# Prepare training data
train_data = []
for index in train_indices:
    row = df.iloc[index]
    train_data.append({
        "prefix": row['Theme'],
        "suffix": get_suffix(row['Valence_mean'], row['Arousal_mean']),
        "image": row['Theme'] + ".jpg"
    })

# Prepare validation data
validation_data = []
for index in validation_indices:
    row = df.iloc[index]
    validation_data.append({
        "prefix": row['Theme'],
        "suffix": get_suffix(row['Valence_mean'], row['Arousal_mean']),
        "image": row['Theme'] + ".jpg"
    })

# Write the training data to train.json
with open('train.json', 'w') as f:
    json.dump(train_data, f, indent=4)

# Write the validation data to validation.json
with open('validation.json', 'w') as f:
    json.dump(validation_data, f, indent=4)
