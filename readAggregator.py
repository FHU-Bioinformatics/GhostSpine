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
        r.generate_U_sequence(r.sequence[:n], r.mods, U_thresh) for r in filtered_reads
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
    
    if len(filtered_reads) == 0:
        st.error("Could not extract any reads, cannot analyze")
        return
    
    make_aggregate_summary_stats(filtered_reads)
    
    with st.spinner("Generating canonical base proportions..."):
        cannonical_base_counts = get_canonical_base_counts(filtered_reads)
        prop = pd.DataFrame(cannonical_base_counts.items(), columns=['Base', 'Count'])
        
        u_free_list = []
        u_bear_list = []
        for r in filtered_reads:
            useq = r.generate_U_sequence(r.sequence, r.mods, U_thresh)
            
            base_free_mask = r.gen_base_free_mask(useq, 3, "U")
            
            read_free_avg = r.get_avg_qscore_by_mask(base_free_mask)
            read_bear_avg = r.get_avg_qscore_by_mask([not i for i in base_free_mask])
            
            if read_free_avg != 0 and read_bear_avg != 0:
                u_free_list.append(read_free_avg)
                u_bear_list.append(read_bear_avg)
            
            
        
        st.write(f"Total free avg: {sum(u_free_list) / len(u_free_list)}")
        st.write(f"Total bear avg: {sum(u_bear_list) / len(u_bear_list)}")
        
        # st.write(f"U-free avg qual. {read.get_avg_qscore_by_mask(read.gen_base_free_mask(uracil_sequence, 3, "U"))}")
        # st.write(f"U-bearing avg qual. {read.get_avg_qscore_by_mask(read.gen_base_bearing_mask(uracil_sequence, 3, "U"))}")
        
        
        
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
        
    