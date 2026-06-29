
import streamlit as st

import tkinter as tk
from tkinter import filedialog
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

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
st.text("v1.0")


def launch_file_picker():
    #file picker fix for mac
    if sys.platform == "darwin":
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
                
                #Don't try to view reads from an old file if a new one is loaded
                if "current_read" in st.session_state:
                    del st.session_state["current_read"]
                if "reads" in st.session_state:
                    del st.session_state["reads"]
    else:
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        
        file_path = filedialog.askopenfilename(master=root, filetypes=[('BAM Files', '*.bam')])
        
        if file_path:
            st.session_state["bam"] = file_path
            
            #Don't try to view reads from an old file if a new one is loaded
            if "current_read" in st.session_state:
                del st.session_state["current_read"]
            if "reads" in st.session_state:
                del st.session_state["reads"]


#Handles all the sidebar widgets and visualization for the analysis of a specific read
def specific_read_analysis():
    
    selection_mode = st.sidebar.segmented_control(
    "Read Extraction Mode", ["Index", "Name"], selection_mode="single", required=True, default="Index",
    help="Extract a read by its index in the Bam file, or its name"
    )
    
    if selection_mode == "Index":
        index_to_extract = st.sidebar.number_input("Index of read (zero-based)", min_value=0, value=0,
                            help="Extract the read of the index listed in this box. Reads use zero-based indexing, meaning they start at 0 instead of 1.")
        
        if st.sidebar.button(f"Search"):
            with st.spinner(f"Searching for read with index {index_to_extract}..."):
                st.session_state["read"] = bamParsing.get_Nth_read(st.session_state["bam"], index_to_extract)
                st.session_state["read_index"] = index_to_extract
        
    elif selection_mode == "Name":
        name_to_extract = st.sidebar.text_input("Name of read",
                            help="Please input the read name carefully. If you make a mistake, Ghost Spine will search the entire file for a read that does not exist.")
        clean_name = name_to_extract.strip().lower()
        
        if st.sidebar.button("Search"):
            with st.spinner(f"Searching for read with name {clean_name}..."):
                st.session_state["read"], st.session_state["read_index"] = bamParsing.get_data_from_read(st.session_state["bam"], clean_name)
    
    if "read" not in st.session_state : return
    
    uracil_confidence_threshold = st.sidebar.slider("Uracil Threshold", 0, 255, 230)
    
    readVisualizer.visualize_read(st.session_state["read"], st.session_state["read_index"], uracil_confidence_threshold)

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

    analysis_mode = st.sidebar.radio("Analysis Mode", ["Specific Read", "Aggregation"], 
                                     help="View the makeup and stats of a single read (Specific Read), or of the entire file (Aggregation)")
    
    if analysis_mode == "Specific Read":
        specific_read_analysis()
    else:
        read_aggregation_analysis()

    

render_sidebar()