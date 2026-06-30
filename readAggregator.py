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
    
def aggregate_file(bam, U_thresh, min_len, max_len):
    with st.spinner("Extracting reads from file, please wait..."):
        filtered_reads, num_reads_in_file = bamParsing.get_everything(bam, min_len, max_len)
    
    #include length filtration bounds if used
    if min_len == max_len == -1:
        st.info(f"Aggregated {len(filtered_reads)} reads from {num_reads_in_file} reads")
    else:
        st.info(f"Aggregated {len(filtered_reads)} reads of length ∈ [{min_len}, {max_len}] from {num_reads_in_file} reads")
    
    make_aggregate_summary_stats(filtered_reads)
    
    with st.spinner("Generating canonical base proportions..."):
        cannonical_base_counts = get_canonical_base_counts(filtered_reads)
        prop = pd.DataFrame(cannonical_base_counts.items(), columns=['Base', 'Count'])
        
        #Make bases always appear alphabetically
        prop.sort_values(by='Base', inplace=True)
        st.write("Canonical Base Proportion")
        st.bar_chart(data = prop, x="Base", y="Count", horizontal=True, color="Base")
        
        #generate U-containing proportion
        total_u_count = sum(read.get_U_count(U_thresh) for read in filtered_reads)
        u_df = prop.copy()
        u_df.loc[u_df['Base'] == "T", 'Count'] -= total_u_count #subtract total U count from total T count
        u_df.loc[len(u_df)] = ['U', total_u_count] #add total U count to df
        u_df.sort_values(by='Base', inplace=True)
        st.write(f"Base Proportion including Uracil (threshold = {U_thresh})")
        st.bar_chart(data = u_df, x="Base", y="Count", horizontal=True, color="Base")
        
    
    if min_len < 50:
        st.warning('''The minimum read length is < 50, the Aggregate Uracil Count per Index 
                   graph will not be generated. To generate this graph, enable read filtration 
                   and set the minimum sequence length to a value >= 50.''')
    else:
        counts = get_sum_of_u_positions(filtered_reads, U_thresh, 50)
        make_U_count_line_graph(counts)
        
    