
import json

def extract_image_descriptions(input_text, output_file):
    # Split the input text into entries
    entries = input_text.strip().split('--------------------')
    output_data = []

    for entry in entries:
        if entry.strip():  # Skip empty sections
            lines = entry.strip().split('\n')
            # Extract image details
            first_line = lines[0]
            image_path = first_line.split(',')[0].split(': ')[1].strip()
            emotion = first_line.split('Emotion: ')[1].strip()

            # Extract description (all lines after the first)
            description = ' '.join(line.strip() for line in lines[1:])

            # Create structured data
            entry_data = {
                "prefix": f"You are a narrator. Describe what you see in this image, using vivid details and evocative language to immerse the viewer in the feeling of {emotion}. Your words should draw the listener deeply into the scene, helping them connect emotionally as they gaze upon the image. Write everything in one paragraph.",
                "suffix": description,
                "image": image_path
            }
            output_data.append(entry_data)

    # Save structured data as JSON file
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, ensure_ascii=False, indent=4)

    print(f"Extracted data has been saved to {output_file}")


# Example usage
if __name__ == "__main__":
    input_text = """
--------------------
Image: oasis/Acorns 1.jpg, Valence: positive, Arousal: low, Emotion: Content
The image depicts a close-up of two acorns resting in the forest, surrounded by a rustic, natural setting.

The acorns sit side by side, their brown, catlike tops facing upward and their light brown bodies lying flat on the forest floor. To the left of the larger acorn, a yellow-brown leaf has broken in half. A thin branch with a dried leaf lies on the right side of the acorn, with a few brown twigs and leaves in the top-left corner adding depth to the scene.

The background reveals a leafy and moss-covered forest floor. The fallen leaves have mounded and pressed the acorns and twigs into the moss, leaving only their tips exposed, as if the forest floor were slowly swallowing the leaves whole. The moss, light green in color, adds a soft and gentle texture to the scene, creating a sense of tranquility and serenity.

The image is reminiscent of the stillness and quiet of the fall season, evoking feelings of comfort and peacefulness in the viewer. The close-up detail of the acorns and surrounding forest floor allows the viewer to immerse themselves in the natural beauty of the scene.<|eot_id|>

--------------------
Image: oasis/Alcohol 7.jpg, Valence: neutral, Arousal: low, Emotion: Calm
As the image unfolds, the visual composition seems to wrap itself around your psyche like a warm, golden embrace. Against the backdrop of a wine rack, softly aglow in a warm amber glow that recalls the golden hour at dusk, a row of empty glass bottles stands guard, their necks neatly capped in dark, matte black caps. Behind the bottles, a small, clear glass bowl lies on the polished wooden countertop, its red wine still suspended inside, awaiting sippers in its translucent bowl, as if anticipating a toast or a solo sip. To the left of the glass, another, larger bottle of wine, similarly capped with black, stands sentry, its own reflective golden surface a promise of rich flavor to come. A few bottles in the middle row also feature a gold and dark label, an extra layer of detail in their rich hues.

On the top row of bottles, a subtle red glow seems to seep through their bottoms, where the glass glows just enough to be seen as they sit in line against the warm, golden background. A glass partition or door divides the counter in the middle, allowing the view from the left to be a reflection of what would have been on the right. In front of the bottles, and to the left and right,

"""

    output_file = "image_descriptions.json"
    extract_image_descriptions(input_text, output_file)
