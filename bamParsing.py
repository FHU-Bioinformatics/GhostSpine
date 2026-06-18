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