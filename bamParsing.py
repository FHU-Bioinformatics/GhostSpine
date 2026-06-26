import bamnostic

def get_N_reads(bam_path, n : int) -> list[str]:
    reads = []
    with bamnostic.AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        i = 0
        for read in bam_file:
            reads.append(read.query_name)
            i += 1
            if i >= int(n):
                break
        return reads

def get_mods_from_read(bam_path, target_read_name):
    with bamnostic.AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        for read in bam_file:
            if read.query_name == target_read_name:
                return list(read.get_tag("ML"))
                    

def get_seq_and_score_from_read(bam_path, target_read_name):
    with bamnostic.AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        for read in bam_file:
            if read.query_name == target_read_name:
                return read.query_sequence, read.query_qualities

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
    

#Determine if this read should be included in aggregation analysis
def is_read_valid_for_aggregation(read : FullRead) -> bool:
    if len(read.sequence) < 80:
        return False
    
    return True


def get_everything(bam_path):
    full_reads : list[FullRead] = []
    num_reads_in_file = 0
    
    with bamnostic.AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        read = next(bam_file)
        num_reads_in_file += 1
        try:
            r = FullRead(read.query_name,
                            read.query_sequence,                                read.query_qualities,
                            list(read.get_tag("ML")))
            if is_read_valid_for_aggregation(r):
                full_reads.append(r)
        except:
            pass

    return full_reads, num_reads_in_file