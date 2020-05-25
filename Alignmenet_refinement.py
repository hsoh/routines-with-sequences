
def sequence_length(fasta):
    num_site = 0
    seq_r = open(fasta, "r")
    line = seq_r.readline()
    while line != '':
        if line.startswith('>'):
            seq_list = []
            line = seq_r.readline()
            while line != '':
                seq_list.append(line.strip())
                line = seq_r.readline()
                if line.startswith('>'):
                    break
            seq_chars = "".join(seq_list)
            num_site = len(seq_chars)
            break
        else:
            line = seq_r.readline()
    seq_r.close()
    return num_site

def sequence_number(fasta):
    num_seq = 0
    seq_r = open(fasta, "r")
    for line in seq_r:
        if line.startswith('>'):
            num_seq += 1
    seq_r.close()
    return num_seq

import argparse

parser = argparse.ArgumentParser(description = "Input specifications ")
parser.add_argument("--i", dest="input_fasta", required=True, help="input fasta file", type=str)
parser.add_argument("--o", dest="output_fasta", required=True, help="output fatsa file", type=str)
parser.add_argument("--maxgap", dest="max_gap_freq", required=False, help="[0..1] max gap freq, default = 1 (allowing all)", type=float, default = 1)
parser.add_argument("--del_low_cov_entry", dest="want_del_low_cov_entry", required=False, help="use this if you want to remove the entries with low coverage after removing gapped sites by --maxgap", action="store_true", default=False)
parser.add_argument("--entry_maxgap", dest="entry_max_gap", required=False, help="used when --del_low_cov_entry; max gap col freq of a entry allowed; default = 0.5 ", type=float, default=0.5)
args = parser.parse_args()

input_fasta = args.input_fasta
output_fasta = args.output_fasta
max_gap_freq = args.max_gap_freq
want_del_low_cov_entry = args.want_del_low_cov_entry
entry_max_gap = args.entry_max_gap


num_site = sequence_length(input_fasta)
num_seq = sequence_number(input_fasta)
site_gap_count = [0]*num_site

print("Number of sites in input = " + str(num_site))
print("Number of strains in input = " + str(num_seq))

seq_r = open(input_fasta, "r")
line = seq_r.readline()
while line != '':
    if line.startswith('>'):
        seq_id = line.strip()[1:]
        seq_list = []
        line = seq_r.readline()
        while line != '':
            seq_list.append(line.strip())
            line = seq_r.readline()
            if line.startswith('>'):
                break
        
        seq_chars = "".join(seq_list)
        for char_idx in range(num_site):
            aligned_char = seq_chars[char_idx]
            if aligned_char in ['.' , '-']:
                site_gap_count[char_idx] += 1
    else:
        line = seq_r.readline()
seq_r.close()

site_gap_profile = []
dict_siteidx_filterpass = {}
num_filtered_site = 0

for char_idx in range(num_site):
    gap_profile = float(site_gap_count[char_idx])/float(num_site)
    if gap_profile <= max_gap_freq:
        num_filtered_site += 1
        dict_siteidx_filterpass[char_idx] = True
    else:
        dict_siteidx_filterpass[char_idx] = False

print("Number of sites filtered by max gap freq = " + str(num_filtered_site))




# sequences of filtered sites
output_w = open(output_fasta, "w")
seq_r = open(input_fasta, "r")
num_filtered_entry = 0
list_excluded_seqs = []

line = seq_r.readline()
while line != '':
    if line.startswith('>'):
        seq_id = line.strip()[1:]
        seq_list = []
        line = seq_r.readline()
        while line != '':
            seq_list.append(line.strip())
            line = seq_r.readline()
            if line.startswith('>'):
                break
        
        seq_chars = "".join(seq_list)
        filtered_seq = []
        num_gap_site_after_filter = 0
        for char_idx in range(num_site):
            if dict_siteidx_filterpass[char_idx]:
                filtered_seq.append(seq_chars[char_idx])
                if seq_chars[char_idx] in ['-', '.']:
                    num_gap_site_after_filter += 1
        gap_freq_in_filtered_entry = float(num_gap_site_after_filter)/float(num_filtered_site)
        
        to_add_this_seq = False
        if not want_del_low_cov_entry:
            to_add_this_seq = True
        else:
            if gap_freq_in_filtered_entry <= max_gap_freq:
                to_add_this_seq = True
        
        if to_add_this_seq:
            output_w.write(">" + seq_id + "\n" + "".join(filtered_seq) + "\n")
            num_filtered_entry += 1
        else:
            list_excluded_seqs.append(seq_id)
        
    else:
        line = seq_r.readline()
seq_r.close()
output_w.close()

num_excluded = len(list_excluded_seqs)
print("Number of entreis included in the output = " + str(num_filtered_entry))
print("Number of entreis excluded from the output due to high gap content = " + str(num_excluded))

if want_del_low_cov_entry:
    exc_w = open(output_fasta + ".excluded.ids", "w")
    for i in range(num_excluded):
        exc_w.write(list_excluded_seqs[i]+"\n")
    exc_w.close()


