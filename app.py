import streamlit as st

import tkinter as tk
from tkinter import filedialog
import sys
import bamParsing
import readVisualizer
import readAggregator

st.set_page_config(
    page_title="Ghost Spine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.subheader("Ghost Spine: The Ghost Shark Inference Viewer")

def launch_file_picker():
    if sys.platform == "win32": 
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        
        file_path = filedialog.askopenfilename(master=root, filetypes=[('BAM Files', '*.bam')])
        
        if file_path:
            st.session_state["bam"] = file_path
    elif sys.platform == "darwin":
        script = '''
                set fileSelected to choose file with prompt "Select a BAM file"
                return POSIX path of fileSelected
                '''
        import subprocess
        result = subprocess.run(['osascript', '-e', script], capture_output=True, text=True)
        if result.returncode == 0:
            file_path =  result.stdout.strip()
            if file_path: 
                st.session_state["bam"] = file_path
    else:
        st.sidebar.error("Why you using this OS?")


def on_extract_reads_button_pressed():
    if "bam" not in st.session_state:
        return
    st.session_state["reads"] = bamParsing.get_N_reads(st.session_state["bam"], st.session_state["reads_to_extract"])


#Handles all the sidebar widgets and visualization for the analysis of a specific read
def specific_read_analysis():
    #Get the number of reads to extract from the bam file and extract them
    st.session_state["reads_to_extract"] = st.sidebar.number_input("Number of reads to extract", min_value = 1, max_value = 9999, value=50)
    extract_button = st.sidebar.button(f"Extract {st.session_state["reads_to_extract"]} reads", on_click=on_extract_reads_button_pressed)
    
    if "reads" not in st.session_state : return
    
    #Load the extracted reads into the dropdown and specify one read to view
    selected_read = st.sidebar.selectbox("Select read", st.session_state["reads"], filter_mode="contains")
    st.session_state["current_read"] = selected_read
    
    if "current_read" not in st.session_state : return

    #Set the uracil threshold and visualize the selected read
    uracil_confidence_threshold = st.sidebar.slider("Uracil Threshold", 0, 255, 230)
    readVisualizer.visualize_read(st.session_state["bam"], st.session_state["current_read"], uracil_confidence_threshold)

def read_aggregation_analysis():
    uracil_confidence_threshold = st.sidebar.slider("Uracil Threshold", 0, 255, 230)
    if st.sidebar.button(f"Run Aggregate Analysis", help="Depending on the file size, this may take some time"):
        readAggregator.aggregate_file(st.session_state["bam"], uracil_confidence_threshold)
    

def render_sidebar():
    select_file_button = st.sidebar.button(f"Select BAM File", on_click=launch_file_picker)

    #Bam file selection
    if "bam" in st.session_state:
        st.sidebar.success(st.session_state["bam"])
    else:
        st.sidebar.warning("Please select a BAM file")

    if "bam" not in st.session_state : return

    analysis_mode = st.sidebar.radio("Analysis Mode", ["Specific Read", "Aggregation"])
    
    if analysis_mode == "Specific Read":
        specific_read_analysis()
    else:
        read_aggregation_analysis()

    

render_sidebar()