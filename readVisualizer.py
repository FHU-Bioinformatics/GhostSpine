import streamlit as st
import pandas as pd

import matplotlib.pyplot as plt

import bamParsing

#The current mod list only contains mod values for each T. Extend the list to include 0 for every other base.
def create_full_mod_list(sequence : str, mod_list : list[int]) -> list[int]:
    full_list = []
    for char in sequence:
        if char == "T":
            next_mod = mod_list.pop(0)
            full_list.append(next_mod)
        else:
            full_list.append(0)
    return full_list

#Create a sequence where all T with a mod confidence > thresh is replaced with U
def create_Uracil_sequence(sequence : str, full_mod_list : list[int], thresh : int) -> list[str]:
    U_seq = []
    for i in range(len(sequence)):
        if sequence[i] != "T":
            U_seq.append(sequence[i])
        else:
            if full_mod_list[i] >= thresh:
                U_seq.append("U")
            else:
                U_seq.append("T")
    return U_seq

#Get the count of Ts whose mod score is >= the uracil threshold
def get_num_uracils(mod_list : list[int], thresh) -> int:
    return sum(1 for mod in mod_list if mod >= thresh)

#Create a dataframe to visualize the components of 
def make_read_vizualization_dataframe(sequence, qualities, mods, thresh):
    df = pd.DataFrame([range(len(sequence))])
    df.loc[len(df)] = [char for char in sequence] #Put canonical sequence into df

    df.drop(0, axis=0, inplace=True) #janky workaround to get the df to display everything horizontally
    df = df.reset_index(drop=True)

    df.loc[len(df)] = [q for q in qualities] #insert q scores

    full_mod_list = create_full_mod_list(sequence, mods.copy())

    df.loc[len(df)] = full_mod_list #Put full mod list into df
    df.loc[len(df)] = create_Uracil_sequence(sequence, full_mod_list, thresh) #Make Uracil converted sequence

    df["Index"] = ["Canonical", "Q-Score", "T+U Mod", "Uracil Seq."]
    df.insert(0, 'Index', df.pop('Index')) #make the first column an "index" column

    return df

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

def visualize_read(bam_path, read_name, thresh) -> None:
    read, index = bamParsing.get_data_from_read(bam_path, read_name)
    
    # read = bamParsing.get_Nth_read(bam_path, 0)
    # index = 0
    
    suspected_uracils = get_num_uracils(read.mods, thresh)

    st.info(f"Currently viewing read: {read.name} // Index: {index} (zero-based)")

    make_summary_stats(read, thresh, suspected_uracils)

    with st.spinner("Building read visualization..."):
        df = make_read_vizualization_dataframe(read.sequence, read.qualities, read.mods, thresh)
        st.dataframe(df, hide_index=True, on_select="ignore", )

    make_U_mod_line_graph(read.mods[:100], "First 100 T reads by T+U Mod Score", thresh)
    make_U_mod_line_graph(read.mods[-100:], "Last 100 T reads by T+U Mod Score", thresh)
    
    # q_score_per_TcallU(sequence ,full_mod_list, qualities, "Q-Scores per base with U", thresh)
    # calculate_Q_score_dis_to_T(sequence, full_mod_list, qualities, "Avg Q-Score by Distance from U", thresh)