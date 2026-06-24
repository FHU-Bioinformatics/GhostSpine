import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import bamParsing

#Make a linechart of U counts at each position
def make_U_count_line_graph(seq) -> None:
    df = pd.DataFrame({'U_per_index': seq})

    fig, ax = plt.subplots(figsize=(15, 5))
    df.plot(kind="line", ax=ax, marker='o',)

    plt.title("Aggregate Uracil Count per Index")
    plt.ylabel("Count")
    plt.grid(axis='y')

    st.pyplot(fig)

def get_sum_of_u_positions(filtered_reads, U_thresh, n : int):    
    #Generate a list that has the first N sequences from each filtered read
    list_of_first_U_seq_reads = [
        r.generate_U_seq(r.sequence[:n], r.mods.copy(), U_thresh) for r in filtered_reads
        ]
    
    #Count the total number of times U appears at each index
    counts = [sum(seq[i] == "U" for seq in list_of_first_U_seq_reads)
    for i in range(n)]
    
    return counts
    

def aggregate_file(bam, U_thresh):
    filtered_reads, num_reads_in_file = bamParsing.get_everything(bam)
    st.info(f"Aggregated {len(filtered_reads)} reads of length 80+ from {num_reads_in_file} reads")
    

    counts = get_sum_of_u_positions(filtered_reads, U_thresh, 50)
    
    make_U_count_line_graph(counts)
    