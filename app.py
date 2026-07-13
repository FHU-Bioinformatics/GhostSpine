
import streamlit as st

import sys
import os
import tempfile

import crossfiledialog
from pathlib import Path

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

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    
    return os.path.join(base_path, relative_path)

def build_header(version : str, date : str):
    title, fhu = st.columns([8, 1])

    with title:
        st.subheader("Ghost Spine: The Ghost Shark Inference Viewer")
        st.text(f"{version} | {date}")
    with fhu:
        st.image(resource_path("icons/fhu_academics.jpg"), width = 200) #max width to 200 because the image gets way too big at low width


def launch_file_picker(file_filter, output_location):
    file_path = crossfiledialog.open_file(filter=file_filter)
        
    if file_path:
        st.session_state[output_location] = file_path

def load_filter_file(file_path: str) -> list[str]:
    with Path(file_path).open("r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


#load read and read_index session states based on selection by index
def select_read_by_index():
    with st.sidebar.form("index_form", enter_to_submit=False):
        index_to_extract = st.number_input("Read Index (zero-based)", min_value=0, value=0,
                            help="Extract the read of the index listed in this box. Reads use zero-based indexing, meaning they start at 0 instead of 1.")
        st.form_submit_button("Search")
        
    with st.spinner(f"Searching for read with index {index_to_extract}..."):
        st.session_state["read"] = bamParsing.get_Nth_read(st.session_state["bam"], index_to_extract)
        st.session_state["read_index"] = index_to_extract

#load read and read_index session states based on selection by name
def select_read_by_name():
    default_text = ""
    
    with st.sidebar.form("name_form", enter_to_submit=False):
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
    
    use_base_coloring = st.sidebar.checkbox("Apply Canonical Base Coloring (slow)", value = False, help="Apply coloring to the Canonical Sequence. Depending on the read length, this could take significantly longer than viewing the read with no coloring.")
    
    if use_base_coloring:
        st.session_state["use_base_coloring"] = True
    else:
        st.session_state["use_base_coloring"] = False
    
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
    
    use_filter_file = st.sidebar.checkbox("Use A Filter File", help="Provide a .txt file of read IDs. Only these read IDs will be used in aggregate analysis.")
    
    if use_filter_file:
        ff_button = st.sidebar.button(f"Select A Filter File", on_click=launch_file_picker, args = ("*.txt", "filter_file"))
        if "filter_file" in st.session_state:
            st.sidebar.success(st.session_state["filter_file"])
        else:
            st.sidebar.warning("Please select a filter file")
    
    uracil_confidence_threshold = st.sidebar.slider("Uracil Threshold", 0, 255, 230, help="The T+U mod score required to consider a T as a U")
    
    if st.sidebar.button(f"Run Aggregate Analysis", help="Depending on the file size, this may take some time"):
        
        if use_filter_file:
            filter_list = load_filter_file(st.session_state["filter_file"]) if "filter_file" in st.session_state else []
        else:
            filter_list = []
        
        readAggregator.aggregate_file(st.session_state["bam"], uracil_confidence_threshold, min_len, max_len, filter_list)
    

def render_sidebar():
    select_file_button = st.sidebar.button(f"Select BAM File", on_click=launch_file_picker, args = ("*.bam", "bam"))
    
    
    ## Using streamlits file uploaded if ever moved to server
    # uploaded_file = st.sidebar.file_uploader("Upload a BAM File", type=["bam"], accept_multiple_files=False)
    
    # if uploaded_file is not None:
        
    #     with tempfile.NamedTemporaryFile(delete=False) as temp_file:
    #         temp_file.write(uploaded_file.getvalue())
    #         temp_filepath = temp_file.name
        
    #     st.session_state["bam"] = temp_filepath

    # st.session_state["bam"] = temp_filepath


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

    
build_header("v1.2.1", "7/9/26")
render_sidebar()