import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
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


def make_aggregate_summary_stats(reads):
    avg_seq_len, min_seq_len, max_seq_len = st.columns(3)
    sequence_lengths = np.array([len(r.sequence) for r in reads])

    avg_seq_len.metric("Mean Read Length", round(sequence_lengths.mean()), border=True)
    min_seq_len.metric("Min Read Length", sequence_lengths.min(), border=True)
    max_seq_len.metric("Max Read Length", sequence_lengths.max(), border=True)


def get_canonical_base_counts(reads) -> Counter:
    base_counts = Counter()
    for r in reads:
        base_counts.update(r.sequence)
    return base_counts
    
def aggregate_file(bam, U_thresh):
    with st.spinner("Extracting reads from file, please wait..."):
        filtered_reads, num_reads_in_file = bamParsing.get_everything(bam)
    
    st.info(f"Aggregated {len(filtered_reads)} reads of length 80+ from {num_reads_in_file} reads")
    
    make_aggregate_summary_stats(filtered_reads)
    
    with st.spinner("Generating canonical base proportions..."):
        cannonical_base_counts = get_canonical_base_counts(filtered_reads)
        prop = pd.DataFrame(cannonical_base_counts.items(), columns=['Base', 'Count'])
        
        #Make bases always appear alphabetically
        prop.sort_values(by='Base', inplace=True)
        st.bar_chart(data = prop, x="Base", y="Count", horizontal=True, color="Base")
    
    counts = get_sum_of_u_positions(filtered_reads, U_thresh, 50)
    make_U_count_line_graph(counts)
    
    