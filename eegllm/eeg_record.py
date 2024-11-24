import muselsl
from pylsl import StreamInlet, resolve_byprop
import time

def stream_eeg_data():
    # First, check if an EEG stream is available
    streams = resolve_byprop('type', 'EEG', timeout=2)
    if not streams:
        print("No EEG stream available. Make sure your Muse is connected and streaming.")
        return

    # Create an inlet to read from the stream
    inlet = StreamInlet(streams[0])

    print("Starting to read EEG stream. Press Ctrl+C to stop.")
    try:
        while True:
            sample, timestamp = inlet.pull_sample(timeout=1.0)
            if timestamp:
                print(f"Timestamp: {timestamp}, Data: {sample}")
    except KeyboardInterrupt:
        print("Stream stopped.")

if __name__ == '__main__':
    # Ensure the Muse device is paired and streaming data. You may need to start the stream via another method or the muselsl viewer.
    stream_eeg_data()
