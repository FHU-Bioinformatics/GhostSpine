
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

title, fhu = st.columns([8, 1])

with title:
    st.subheader("Ghost Spine: The Ghost Shark Inference Viewer")

    st.text("v1.1.0 | 6/30/2026")

with fhu:
    st.image("icons/fhu_academics.jpg", width = 200) #max width to 200 because the image gets way too big at low width


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
                
    else:
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        
        file_path = filedialog.askopenfilename(master=root, filetypes=[('BAM Files', '*.bam')])
        
        if file_path:
            st.session_state["bam"] = file_path


#load read and read_index session states based on selection by index
def select_read_by_index():
    with st.sidebar.form("index_form"):
        index_to_extract = st.number_input("Read Index (zero-based)", min_value=0, value=0,
                            help="Extract the read of the index listed in this box. Reads use zero-based indexing, meaning they start at 0 instead of 1.")
        st.form_submit_button("Search")
        
    with st.spinner(f"Searching for read with index {index_to_extract}..."):
        st.session_state["read"] = bamParsing.get_Nth_read(st.session_state["bam"], index_to_extract)
        st.session_state["read_index"] = index_to_extract

#load read and read_index session states based on selection by name
def select_read_by_name():
    
    default_text = ""
    
    with st.sidebar.form("name_form"):
        name_to_extract = st.text_input("Read Name", value=default_text,
                            help="Please input the read name carefully. If you make a mistake, Ghost Spine will search the entire file for a read that does not exist.")
        clean_name = name_to_extract.strip().lower()
        st.form_submit_button("Search")
    
    if name_to_extract == default_text : return #prevent immediately trying to search for read with invalid name
    
    with st.spinner(f"Searching for read with name {clean_name}..."):
        st.session_state["read"], st.session_state["read_index"] = bamParsing.get_data_from_read(st.session_state["bam"], clean_name)
    

#Handles all the sidebar widgets and visualization for the analysis of a specific read
def specific_read_analysis():
    
    selection_mode = st.sidebar.segmented_control(
    "Read Extraction Mode", ["Index", "Name"], selection_mode="single", required=True, default="Index",
    help="Extract a read by its index in the Bam file, or its name"
    )
    
    if selection_mode == "Index":
        select_read_by_index()
    elif selection_mode == "Name":
        select_read_by_name()
    
    if "read" not in st.session_state : return
    
    uracil_confidence_threshold = st.sidebar.slider("Uracil Threshold", 0, 255, 230, help="The T+U mod score required to consider a T as a U")
    
    readVisualizer.visualize_read(st.session_state["read"], st.session_state["read_index"], uracil_confidence_threshold)

def read_aggregation_analysis():
    use_filtration = st.sidebar.checkbox("Apply Read Filtration", help="Exclude reads from aggregate analysis by min and max sequence length")
    
    if use_filtration:
        min_len = st.sidebar.number_input("Ignore reads shorter than:", value = 80, min_value=0)
        max_len = st.sidebar.number_input("Ignore reads longer than:", value = 999999, min_value= 0)
    else:
        min_len = -1
        max_len = -1
    
    uracil_confidence_threshold = st.sidebar.slider("Uracil Threshold", 0, 255, 230, help="The T+U mod score required to consider a T as a U")
    if st.sidebar.button(f"Run Aggregate Analysis", help="Depending on the file size, this may take some time"):
        readAggregator.aggregate_file(st.session_state["bam"], uracil_confidence_threshold, min_len, max_len)
    

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