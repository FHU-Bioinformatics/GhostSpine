import pandas as pd

class FullRead:
    def __init__(self, name, sequence, qualities, mods = []):
        self.name = name
        self.sequence = sequence
        self.qualities = qualities
        self.mods = mods
    
    #Uses the mod scores to create a u containing sequence
    def generate_U_sequence(self, seq, mods, U_thresh):
        new_seq = []
        i = 0
        for base in seq:
            if base != "T":
                new_seq.append(base)
            else:
                next_mod = mods[i]
                i += 1
                if next_mod >= U_thresh:
                    new_seq.append("U")
                else:
                    new_seq.append("T")
        return new_seq

    #Uses the mod scores and a threshold to see how many uracils are suspected in a read
    def get_U_count(self, thresh):
        return sum(1 for mod in self.mods if mod >= thresh)
    
    #Returns the a list[bool] with a length equivalent to the length of the sequence
    #The list will be True at all indexes where the base at said index is > n positions away from the specified base
    def gen_base_free_mask(self, sequence, n: int, base : str) -> list[bool]:
        base_free = [True] * len(sequence)

        for i, value in enumerate(sequence):
            if value == base:
                start = max(0, i - n)
                end = min(len(sequence), i + n + 1)  # end is exclusive
                for j in range(start, end):
                    base_free[j] = False

        return base_free
    