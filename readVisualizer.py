import streamlit as st
import pandas as pd
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
def get_suspected_uracils(mod_list : list[int], thresh) -> int:
    uracil_sum = 0
    for mod in mod_list:
        if mod > thresh:
            uracil_sum += 1
    return uracil_sum

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

def visualize_read(bam_path, read_name, thresh) -> None:
    mods = bamParsing.get_mods_from_read(bam_path, read_name)
    sequence, qualities = bamParsing.get_seq_and_score_from_read(bam_path, read_name)

    suspected_uracils = get_suspected_uracils(mods, thresh)

    st.warning(f"Currently viewing read: {read_name}")
    st.info(f"{len(sequence)} Bases // Avg Q-Score: {round(sum(qualities) / len(qualities), 1)} // Uracil Count (>= {thresh} threshold): {suspected_uracils}")

    df = make_read_vizualization_dataframe(sequence, qualities, mods, thresh)

    st.dataframe(df, hide_index=True, on_select="ignore", )

