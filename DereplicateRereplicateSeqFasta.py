
def identical_strings(seq_a, seq_b):
    identity = True
    len_a = len(seq_a)
    len_b = len(seq_b)
    if len_a != len_b:
        identity = False
    else:
        for position in range(len_a):
            if seq_a[position] != seq_b[position]:
                identity = False
                break
    return identity


import argparse

parser = argparse.ArgumentParser(description = "Input specifications ")
parser.add_argument("--i", dest="inputFasta", required=True, help="input fasta file", type=str)
parser.add_argument("--o", dest="outputFasta", required=True, help="output fatsa file", type=str)
parser.add_argument("--repmap", dest="repmapFile", required=True, help="mapping tsv between original ID - deplicated rep ID (=is input of --rereplicate, =is output of --dereplicate)", type=str)
parser.add_argument("--do_derep", dest="do_derep", required=False, help="run dereplication", action="store_true", default=False)
parser.add_argument("--do_rerep", dest="do_rerep", required=False, help="run rereplication", action="store_true", default=False)
args = parser.parse_args()

inputFasta = args.inputFasta
outputFasta = args.outputFasta
repmapFile = args.repmapFile
do_derep = args.do_derep
do_rerep = args.do_rerep

if do_derep:
    if do_rerep:
        print("you can do one thing at a time, not both, man.")
        exit

if do_derep:
    # read fasta
    list_id = []
    dict_id_seq = {}
    seqr = open(inputFasta, "r")
    line = seqr.readline()
    while line != '':
        if line.startswith(">"):
            seq_id_field = line.strip()[1:]
            seq_id = seq_id_field
            if seq_id_field.find(" ") > -1:
                seq_id = seq_id_field.split(" ")[0]
            list_id.append(seq_id)
            line = seqr.readline()
            seq_list = []
            while line != '':
                seq_list.append(line.strip())
                line = seqr.readline()
                if line.startswith(">"):
                    break
            seq = "".join(seq_list)
            dict_id_seq[seq_id] = seq
        else:
            line = seqr.readline()    
    seqr.close()

    dict_id_repness = {}    # [id] = True if is rep, [id] = False if is non-rep
    dict_id_rep = {}        # [id] = rep id
    num_id = len(list_id)
    for idx in range(num_id):
        seqid = list_id[idx]
        dict_id_repness[seqid] = True
        dict_id_rep[seqid] = seqid
    print("received " + str(num_id) + " sequences ")

    for repindex in range(num_id -1):
        rep_seqid = list_id[repindex]
        if not dict_id_repness[rep_seqid]:
            continue
        rep_seq = dict_id_seq[rep_seqid]

        for qindex in range(repindex + 1, num_id):
            query_seqid = list_id[qindex]
            if not dict_id_repness[query_seqid]:
                continue
            query_seq = dict_id_seq[query_seqid]
            if identical_strings(query_seq, rep_seq):
                dict_id_repness[query_seqid] = False
                dict_id_rep[query_seqid] = rep_seqid

    num_repseq = 0
    repmap_w = open(repmapFile, "w")
    fasta_w = open(outputFasta, "w")
    for idx in range(num_id):
        seq_id = list_id[idx]
        if dict_id_repness[seq_id]:
            num_repseq += 1
            fasta_w.write(">" + seq_id +"\n" + dict_id_seq[seq_id] + "\n")
        mother_seq_id = dict_id_rep[seq_id]
        repmap_w.write(seq_id + "\t" + mother_seq_id + "\n")
    repmap_w.close()
    fasta_w.close()
    print("dereplicate: " + str(num_repseq) + " unique sequences ")



if do_rerep:
    # read fasta
    list_id = []
    dict_id_seq = {}
    seqr = open(inputFasta, "r")
    line = seqr.readline()
    while line != '':
        if line.startswith(">"):
            seq_id_field = line.strip()[1:]
            seq_id = seq_id_field
            if seq_id_field.find(" ") > -1:
                seq_id = seq_id_field.split(" ")[0]
            list_id.append(seq_id)
            line = seqr.readline()
            seq_list = []
            while line != '':
                seq_list.append(line.strip())
                line = seqr.readline()
                if line.startswith(">"):
                    break
            seq = "".join(seq_list)
            dict_id_seq[seq_id] = seq
        else:
            line = seqr.readline()    
    seqr.close()
    num_id = len(list_id)

    # derplication map
    dict_mother_doughterList = {}
    mapr = open(repmapFile, "r")
    for line in mapr:
        fields = line.strip().split("\t")
        mother = fields[1]
        doughter = fields[0]
        doughter_list = []
        if mother in dict_mother_doughterList:
            doughter_list = dict_mother_doughterList[mother]
        doughter_list.append(doughter)
        dict_mother_doughterList[mother] = doughter_list
    mapr.close()    

    # write re-replicated fasta 
    fasta_w = open(outputFasta, "w")
    for idx in range(num_id):
        seq_id = list_id[idx]
        seq = dict_id_seq[seq_id]
        doughter_list = dict_mother_doughterList[seq_id]
        for doughter in doughter_list:
            fasta_w.write(">" + doughter +"\n" + seq + "\n")
    fasta_w.close()


