
import json

def convert_json_to_jsonl(input_file, output_file):
    # Read the input JSON file
    with open(input_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Write each entry to the JSONL file
    with open(output_file, 'w', encoding='utf-8') as f:
        for entry in data:
            f.write(json.dumps(entry, ensure_ascii=False) + '\n')
    
    print(f"Converted JSON file has been saved to {output_file}")

# Example usage
if __name__ == "__main__":
    input_file = "image_annotation.json"  # Replace with your JSON file path
    output_file = "image_annotation.jsonl"  # Replace with your desired JSONL file path
    convert_json_to_jsonl(input_file, output_file)
