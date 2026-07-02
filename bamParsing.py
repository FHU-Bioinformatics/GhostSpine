from bamnostic import AlignmentFile
import streamlit as st

#A class to hold all attributes of a read for easier aggregate analysis
class FullRead:
    def __init__(self, name, sequence, qualities, mods = []):
        self.name = name
        self.sequence = sequence
        self.qualities = qualities
        self.mods = mods
    
    #Uses the mod scores to convert a sequence to a u containing sequence
    def generate_U_seq(self, seq, mods, U_thresh):
        new_seq = []
        for base in seq:
            if base != "T":
                new_seq.append(base)
            else:
                next_mod = mods.pop(0)
                if next_mod >= U_thresh:
                    new_seq.append("U")
                else:
                    new_seq.append("T")
        return new_seq

    #Uses the mod scores and a threshold to see how many uracils are suspected in a read
    def get_U_count(self, thresh):
        return sum(1 for mod in self.mods if mod >= thresh)
    
                    
#Extract the name, sequence, qscores, and mods from a read and put them in a FullRead object
def get_data_from_read(bam_path, target_read_name) -> tuple[FullRead, int]:
    index = 0
    with AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        for read in bam_file:
            if read.query_name == target_read_name:
                try:

                    r = FullRead(read.query_name,
                                    read.query_sequence,
                                    read.query_qualities,
                                    list(read.get_tag("ML"))
                                    )
                    
                    return r, index

                except:
                    st.warning("Could not extract all data from the selected read. This read may be missing key data.")
                    raise KeyError("Could not extract all data from the selected read. This read may be missing key data.")
            else:
                index += 1
    st.warning(f"Could not find {target_read_name} in the selected file.")
    raise KeyError(f"Could not find {target_read_name} in the selected file.")

#The AlignmentFile object doesn't support bam_file[n], so the Nth read can only be extracted from iteration
#same as get_data_from_read(), but it searches for a given index rather than a given name
def get_Nth_read(bam_path, n : int) -> tuple[FullRead]:
    current_index = 0
    with AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        for read in bam_file:
            if current_index == n:
                try:
                    r = FullRead(read.query_name,
                                read.query_sequence,
                                read.query_qualities,
                                list(read.get_tag("ML"))
                                )
                            
                    return r

                except:
                    st.warning("Could not extract all data from the selected read. This read may be missing key data.")
                    raise KeyError("Could not extract all data from the selected read. This read may be missing key data.")
            else:
                current_index += 1
    st.warning("Could not extract all data from the selected read. This read may be missing key data.")
    raise IndexError("The specified index is higher than the number of reads in the Bam file")

#Determine if this read should be included in aggregation analysis
def is_read_valid_for_aggregation(read : FullRead, min_len, max_len) -> bool:
    seq_len = len(read.sequence)
    if seq_len < min_len or seq_len > max_len:
        return False
    
    return True

#Called in aggregate analysis mode
def get_everything(bam_path, min_len, max_len):
    full_reads : list[FullRead] = []
    num_reads_in_file = 0
    
    with AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        for read in bam_file:
            num_reads_in_file += 1

            #A read isn't guarenteed to have a mod tag for some reason
            try:

                r = FullRead(read.query_name,
                                read.query_sequence,
                                read.query_qualities,
                                list(read.get_tag("ML"))
                                )

                #-1 used as fallback if read filtration isn't enabled
                if min_len != -1 and max_len != -1:
                    if is_read_valid_for_aggregation(r, min_len, max_len):
                        full_reads.append(r)
                else:
                    #append without checking validity
                    full_reads.append(r)

            except:
                pass

    return full_reads, num_reads_in_file