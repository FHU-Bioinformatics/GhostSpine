# Ghost Spine
A program to visualize the inferences made by Ghost Shark

## Usage
In the directory where GhostSpine is installed, run `streamlit run app.py`. This will launch a browser window with the program running locally. If you need a compiled version, see the Compiled Versions section of the README.

Assuming you have a BAM file that was run through Ghost Shark, select it inside Ghost Spine and choose an analysis mode.

### Specific Read Mode
This mode is for the inspection of individual reads in the BAM file. Select a read by either inputting its index or name. Once a read is selected, a Uracil threshold can be set. This will show the read length, average Q-score, and Uracil cound, as well as a chart displaying the contents of the read itself. The canonical sequence is the original sequence present. The Q-score shows the Q-score for each base. The T+U mod shows the Ghost Shark mod confidence score for each T. The Uracil sequence converts the canonical sequence to a sequence where each T with a T+U Mod score above the uracil threshold is converted to uracil.

### Aggregate Mode
This mode is for viewing aggregate statistics of every read in the BAM file. It will display the mean and range of read lengths
in the selected file, as well as base proportions and the aggregate quantities of uracil in the first 50 positions.

## Custom Build
1. `streamlit-desktop-app build app.py --name "GhostSpine" --pyinstaller-options --collect-data="bamnostic"`
2. Copy the `.streamlit` directory and paste it in the root directory of `dist/GhostSpine`

## Compiled Versions
Compiled versions of GhostSpine are available, but they are generally inferior to running GhostSpine in a Python environment. Only use the compiled versions of GhostSpine if setting up an environment is impractical.

The downside of compiling a python codebase to an exe is that said exe is flagged by the antivirus and is prevented from running. For Windows, select the downloaded Zip file, right click, and navigate to properties. In the General tab, find security, click unblock, and click apply. Then extract, navigate to the exe, and run it as normal. This step is not required if running GhostSpine in a Python environment.