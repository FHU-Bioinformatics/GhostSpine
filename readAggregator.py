import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bamParsing

#Make a linechart of U counts at each position
def make_U_count_line_graph(seq) -> None:
    df = pd.DataFrame({'U_per_index': seq})

    fig, ax = plt.subplots(figsize=(15, 5))
    df.plot(kind="line", ax=ax, marker='o',)

    plt.title("Aggregate Uracil Count per Index")
    plt.ylabel("Count")
    # plt.xticks([i + 1 for i in range(20)]) 
    plt.grid(axis='y')

    st.pyplot(fig)


def aggregate_file(bam, U_thresh):
    filtered_reads, num_reads_in_file = bamParsing.get_everything(bam)
    
    st.info(f"Aggregated {len(filtered_reads)} reads of length 80+ from {num_reads_in_file} reads")
    
    list_of_first_twenty_U_seq_reads : list[list[str]] = []
    #Generate the uracil sequence for the first 20 bases of each filtered read
    for read in filtered_reads:
        read.first_twenty_U_seq = read.generate_U_seq(read.sequence, read.mods.copy(), U_thresh)
        list_of_first_twenty_U_seq_reads.append(read.first_twenty_U_seq)
    
    #Count the total number of times U appears at each index
    counts = [sum(seq[i] == "U" for seq in list_of_first_twenty_U_seq_reads)
    for i in range(50)]
    
    make_U_count_line_graph(counts)
    
    

