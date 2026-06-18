import streamlit as st

import tkinter as tk
from tkinter import filedialog

import bamParsing
import readVisualizer

st.set_page_config(
    page_title="Ghost Spine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.subheader("Ghost Spine: Ghost Shark Inference Viewer")


def launch_file_picker():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    
    file_path = filedialog.askopenfilename(master=root, filetypes=[('BAM Files', '*.bam')])
    
    if file_path:
        st.session_state["bam"] = file_path


select_file_button = st.sidebar.button(f"Select BAM File", on_click=launch_file_picker)

if "bam" in st.session_state:
    st.sidebar.success(st.session_state["bam"])
else:
    st.sidebar.warning("Please select a BAM file")



def on_extract_reads_button_pressed():
    if "bam" not in st.session_state:
        print("No BAM File Selected")
        return
    st.session_state["reads"] = bamParsing.get_N_reads(st.session_state["bam"], reads_to_extract)


if "bam" in st.session_state:

    reads_to_extract = st.sidebar.number_input("Number of reads to extract", min_value = 1, max_value = 9999, value=50)
    extract_button = st.sidebar.button(f"Extract {reads_to_extract} reads", on_click=on_extract_reads_button_pressed)
    if "reads" in st.session_state:
        selected_read = st.sidebar.selectbox("Select read", st.session_state["reads"])
        st.session_state["current_read"] = selected_read
    if "current_read" in st.session_state:
        uracil_confidence_threshold = st.sidebar.slider("Uracil Threshold", 0, 255, 230)
        readVisualizer.visualize_read(st.session_state["bam"], st.session_state["current_read"], uracil_confidence_threshold)

