# Ghost Spine
A program to visualize the inferences made by Ghost Shark

## Required packages
  - bamnostic
  - streamlit
  - pandas

## Usage
In the directory where GhostSpine is installed, run `streamlit run app.py`. This will launch a browser window with the program running locally.

Assuming you have a BAM file that was run through Ghost Shark, select it inside Ghost Spine and select how many reads to extract, starting from the beginning of the file. All
the reads will be loaded into the read selection dropdown, and the selected read will be rendered in the main window. The canonical sequence, Q-scores, T+U mod scores, and 
the Uracil-containing sequence will be rendered out. The Uracil confidence threshold will determine how high the T+U mod score must be for the T to be converted to a U in the 
Uracil-containing sequence.
