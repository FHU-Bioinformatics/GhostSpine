import streamlit as st
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from scripts.fullRead import FullRead

st.elements.lib.pandas_styler_utils._use_display_values = lambda df, style: df.astype(str)

# #The current mod list only contains mod values for each T. Extend the list to include 0 for every other base.
def create_full_mod_list(sequence : str, mod_list : list[int]) -> list[int]:
    full_list = []
    i = 0
    for char in sequence:
        if char == "T":
            #should never have out of range exception
            full_list.append(mod_list[i])
            i += 1
        else:
            full_list.append(0)
    
    return full_list


#Generate relevant information and insert into a dataframe for easy read analysis
def make_read_vizualization_dataframe(read : FullRead, thresh):
    
    full_mod_list = create_full_mod_list(read.sequence, read.mods)
    
    uracil_sequence = read.generate_U_sequence(read.sequence, read.mods, thresh)
    
    df = pd.DataFrame({
        "Canonical" : list(read.sequence),
        "Q-Score" : list(read.qualities),
        "T+U Mod" : full_mod_list,
        "Uracil Seq." : uracil_sequence,
        "T-Free" : read.gen_base_free_mask(read.sequence, 3, "T"),
        "U-Free" : read.gen_base_free_mask(uracil_sequence, 3, "U"),
        # "A-Free" : read.gen_base_free_mask(read.sequence, 3, "A")
    })

    return df

def apply_base_highlighting(base):
    
    base_color_dict = {"A" : "mediumblue", "T" : "forestgreen", "C" : "red", "G" : "orange", "U" : "fuchsia"}
    
    # color = 'purple' if base == "U" else ''
    color = base_color_dict[base]
    return f'background-color: {color}'

def apply_U_highlighting(base):
    color = "fuchsia" if base == "U" else ''
    return f'background-color: {color}'

# 1. Define a function that operates on an entire row
def style_target_row(row, target_idx):
    base_color_dict = {"A" : "mediumblue", "T" : "forestgreen", "C" : "red", "G" : "orange", "U" : "fuchsia"}
    
    if row.name != target_idx:
        return [''] * len(row)
    
    # Apply styling logic to the target row values
    return [
        f'background-color: {base_color_dict[base]}' for base in row if base in base_color_dict
    ]

def style_row_with_dict(data, target_idx, mapping):
    # Create an empty styling matrix matching the shape of the data
    style_df = pd.DataFrame('', index=data.index, columns=data.columns)
    
    # Map the dictionary directly to the target row, filling missing keys with a blank string
    styled_row = data.loc[target_idx].map(mapping).fillna('')
    
    # Insert the styled row back into our empty matrix
    style_df.loc[target_idx] = styled_row
    return style_df

def make_U_mod_line_graph(mods, title : str, thresh) -> None:
    df = pd.DataFrame({'Mod Score': mods})

    fig, ax = plt.subplots(figsize=(15, 5))
    df.plot(kind="line", ax=ax, marker='o',)

    plt.title(title)
    plt.ylabel("T+U Mod Score")
    plt.ylim(0, 255)
    plt.grid(axis='y')
    plt.axhline(y=thresh, color='r', linestyle='--', linewidth=2)

    st.pyplot(fig)


def q_score_per_TcallU(sequence, full_mod_list, qualities, title: str, thresh) -> None:
    base_colors = {
        'A': '#2ecc71',   # green
        'T': '#3498db',   # blue
        'U': '#9b59b6',   # purple — T called as U
        'G': '#e67e22',   # orange
        'C': '#e74c3c',   # red
    }

    
    selected_bases = []
    cols = st.columns(5)
    if cols[0].checkbox("A", value=True, key="show_A"):
        selected_bases.append("A")
    if cols[1].checkbox("T", value=True, key="show_T"):
        selected_bases.append("T")
    if cols[2].checkbox("U", value=True, key="show_U"):
        selected_bases.append("U")
    if cols[3].checkbox("G", value=True, key="show_G"):
        selected_bases.append("G")
    if cols[4].checkbox("C", value=True, key="show_C"):
        selected_bases.append("C")
    positions, q_scores, colors, labels = [], [], [], []

    for i, (base, mod, q) in enumerate(zip(sequence, full_mod_list, qualities)):
        called = 'U' if (base == 'T' and mod >= thresh) else base
        if called not in selected_bases:
            continue
        positions.append(i)
        q_scores.append(q)
        colors.append(base_colors.get(called, '#aaaaaa'))
        labels.append(called)

    fig, ax = plt.subplots(figsize=(15, 5))
    ax.plot(positions, q_scores, color='#cccccc', linewidth=0.8, zorder=1)
    ax.scatter(positions, q_scores, c=colors, s=30, zorder=2)

    seen = {}
    for label, color in zip(labels, colors):
        if label not in seen:
            seen[label] = plt.Line2D([0], [0], marker='o', color='w',
                                     markerfacecolor=color, markersize=8, label=label)
    ax.legend(handles=list(seen.values()))

    plt.title(title)
    plt.ylabel("Q-Score")
    plt.xlabel("Position in Read")
    plt.grid(axis='y')
    st.pyplot(fig)

    u_positions = [i for i, mod in enumerate(full_mod_list) if mod >= thresh]

    if not u_positions:
        st.warning("No uracils called at this threshold.")
        return

def calculate_Q_score_dis_to_T(sequence, full_mod_list, qualities, title: str, thresh, max_dist=7) -> None:

    u_positions = [i for i, mod in enumerate(full_mod_list) if mod >= thresh]

    if not u_positions:
        st.warning("No uracils called at this threshold.")
        return

    distance_base_q = {}

    for i, (base, q) in enumerate(zip(sequence, qualities)):
        dist = min(abs(i - u_pos) for u_pos in u_positions)
        if dist == 0 or dist > max_dist:
            continue
        if dist not in distance_base_q:
            distance_base_q[dist] = {}
        if base not in distance_base_q[dist]:
            distance_base_q[dist][base] = []
        distance_base_q[dist][base].append(q)

    rows = []
    for dist in range(1, max_dist + 1):
        row = {"Distance": dist}
        for base in ['A', 'C', 'T', 'G']:
            scores = distance_base_q.get(dist, {}).get(base, [])
            row[base] = round(sum(scores) / len(scores), 1) if scores else "-"
        rows.append(row)

    df = pd.DataFrame(rows).set_index("Distance")
    st.subheader(title)
    st.dataframe(df)
    

def make_summary_stats(read, thresh, suspected_uracils):
    bp_length, qscore, suspected_u = st.columns(3)

    bp_length.metric("Read Length", len(read.sequence), border=True)
    qscore.metric("Average Q-Score", round(sum(read.qualities) / len(read.qualities), 1), border=True)
    suspected_u.metric(f"Uracil Count (>= {thresh})", suspected_uracils, border=True)

#Show visualization df of t-free regions
def build_t_free_analysis(t_free):
    
    if len(t_free) == 0:
        st.warning("The sequence contained no region where a base was 6 or more positions away from a T, nothing to display")
        return
    
    t_free = t_free.drop(columns=["T+U Mod", "Uracil Seq.", "T-Free", "U-Free"]) #Useless columns when no Ts present
        
    st.subheader("T-free Sequence (dist=3)")
    st.dataframe(t_free.T)
    
#Create qscore summary stats of both regions and compare them
def compare_free_and_bearing(free_region, bearing_region, free_name : str, bearing_name : str):
    
    if len(bearing_region) == 0:
        st.warning(f"A {bearing_name} region in this sequence does not exist")
        return
    if len(bearing_region) == 0:
        st.warning(f"A {free_name} region in this sequence does not exist")
        return
    
    free_qscore, bearing_qscore, diff, stat_test = st.columns(4)
    
    free_avg_qscore = round(free_region["Q-Score"].mean(), 1)
    bearing_avg_qscore = round(bearing_region["Q-Score"].mean(), 1)
    qscore_diff = round(free_avg_qscore - bearing_avg_qscore, 1)
    
    free_qscore.metric(f"Mean Q-score of {free_name} Region", free_avg_qscore, f"Length: {len(free_region)}", border=False,
                       delta_color="off", delta_arrow="off")
    
    bearing_qscore.metric(f"Mean Q-score of {bearing_name} Region", bearing_avg_qscore, f"Length: {len(bearing_region)}", border=False,
                          delta_color="off", delta_arrow="off")
    
    diff.metric("Difference of Q-scores", qscore_diff, f"{free_name} avg. - {bearing_name} avg.", border = False,
                delta_color="off", delta_arrow="off")
    
    from scipy import stats
    t_stat, p_val = stats.ttest_ind(free_region["Q-Score"], bearing_region["Q-Score"], equal_var=False)
    stat_test.metric("P-value", f"{p_val:.8f}", f"Welch's t-stat: {t_stat:.8f}", border = False, delta_color="off", delta_arrow="off")

    
#build a relative frequency histogram of all qscores in a df
def make_qscore_distribution_hist(df, title):
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.grid(axis='y', alpha=0.7)
        sns.histplot(data=df, x="Q-Score", stat="percent", binwidth=2, binrange=(0, 50))
        
        ax.set_title(f"Distribution of Q-Scores in {title} Regions (n={len(df)})")
        ax.set_xlim(0, 50)
        ax.set_ylim(0, 100)
        st.pyplot(fig)

def visualize_read(read : FullRead, read_index, thresh) -> None:

    st.info(f"Currently viewing read: {read.name} // Index: {read_index} (zero-based)")
    
    uracil_count = read.get_U_count(thresh)

    make_summary_stats(read, thresh, uracil_count)

    with st.spinner("Building read visualization..."):
        df = make_read_vizualization_dataframe(read, thresh)
        
        
        if st.session_state["use_base_coloring"]:
        
            base_color_dict = {"A" : "background-color: mediumblue",
                            "T" : "background-color: forestgreen",
                            "C" : "background-color: red",
                            "G" : "background-color: orange",
                            "U" : "background-color: fuchsia"}
            
            u_only_dict = {"U" : "background-color: fuchsia"}
            
    
            df_styled = (df.T).style.apply(style_row_with_dict, target_idx="Canonical", mapping=base_color_dict, axis=None).apply(style_row_with_dict, target_idx="Uracil Seq.", mapping=u_only_dict, axis=None)
            st.dataframe(df_styled, hide_index=False, on_select="ignore")
        else:
            st.dataframe(df.T, hide_index=False, on_select="ignore")
        
        
        t_free = df[df["T-Free"] == True].copy()
        t_bearing = df[df["T-Free"] == False].copy()
        
        build_t_free_analysis(t_free)
        
        compare_free_and_bearing(t_free, t_bearing, "T-free", "T-bearing")
        
        u_free = df[df["U-Free"] == True].copy()
        u_bearing = df[df["U-Free"] == False].copy()
        
        compare_free_and_bearing(u_free, u_bearing, "U-free", "U-bearing")
        
        
        
        make_qscore_distribution_hist(u_free, "U-Free")
        if uracil_count == 0:
            st.warning("There are no suspected Uracils in this read, the U-bearing histogram will not be generated")
        else:
            make_qscore_distribution_hist(u_bearing, "U-Bearing")




    make_U_mod_line_graph(read.mods[:100], "First 100 T reads by T+U Mod Score", thresh)
    make_U_mod_line_graph(read.mods[-100:], "Last 100 T reads by T+U Mod Score", thresh)
    
    # q_score_per_TcallU(sequence ,full_mod_list, qualities, "Q-Scores per base with U", thresh)
    # calculate_Q_score_dis_to_T(sequence, full_mod_list, qualities, "Avg Q-Score by Distance from U", thresh)