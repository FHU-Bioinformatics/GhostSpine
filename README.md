# Ghost Spine
A program to visualize the inferences made by Ghost Shark

## Usage
In the directory where GhostSpine is installed, run `streamlit run app.py`. This will launch a browser window with the program running locally.

Assuming you have a BAM file that was run through Ghost Shark, select it inside Ghost Spine and choose an analysis mode.

### Specific Read Mode
This mode is for the inspection of individual reads in the BAM file. Select a number of reads to extract from the beginning of the file, and they will be loaded into a dropdown menu. You can either select from this dropdown, or search for an extracted read. Once a read is 
selected, a Uracil threshold can be set. This will show the read length, average Q-score, and Uracil cound, as well as a chart displaying 
the contents of the read itself. The canonical sequence is the original sequence present. The Q-score shows the Q-score for each base. The T+U mod shows the Ghost Shark mod confidence score for each T. The Uracil sequence converts the canonical sequence to a sequence where each T with a T+U Mod score above the uracil threshold is converted to uracil.

### Aggregate Mode
This mode is for viewing aggregate statistics of every read in the BAM file. It currently allows the visualization of the quantities of
uracil at each of the first 50 positions in each read, aggregated into a single line graph.

## Custom Build
1. `streamlit-desktop-app build app.py --name "GhostSpine" --pyinstaller-options --collect-data="bamnostic"`
2. Copy the `.streamlit` directory and paste it in the root directory of `dist/GhostSpine`