# parse input
import argparse
import os
import subprocess

def blastout_to_dict_query_ref(blast_file):
    dict_q_r = {}
    blast_read = open(blast_file, "r")
    for line in blast_read:
        fields = line.strip().split("\t")
        # only if the first hit for the query
        if fields[0] not in dict_q_r:
            dict_q_r[fields[0]] = fields[1]
    blast_read.close()
    return dict_q_r

def blastout_to_dict_ref_query(bdbh_map_file):
    dict_r_q = {}
    blast_read = open(bdbh_map_file, "r")
    for line in blast_read:
        fields = line.strip().split("\t")
        # only if the first hit for the query
        if fields[1] not in dict_r_q:
            dict_r_q[fields[1]] = fields[0]
    blast_read.close()
    return dict_r_q

def list_up_core_gene_ref_protein(ref_gene_genome_freq_file, core_freq_cutoff):
    list_core_ref_protein = []
    tsv_r = open(ref_gene_genome_freq_file, "r")
    for line in tsv_r:
        fields = line.strip().split("\t")
        ref_protein = fields[0]
        freq = float(fields[2])
        if freq >= core_freq_cutoff:
            list_core_ref_protein.append(ref_protein)
    tsv_r.close()
    return list_core_ref_protein


parser = argparse.ArgumentParser(description = "Input specifications ")
parser.add_argument("--run_bdbh", dest="to_run_bdbh", required=False, help="run bdbh (bidirectional blastp) using diamond [REQUIRES 'diamond' in your path env] for each query proteome vs. reference proteome", default = False, action="store_true")
#   from running this command you will create under bdbh_base_dir the followings:
#           bdbh_base_dir/q2r_blastp
#           bdbh_base_dir/r2q_blastp
#           bdbh_base_dir/q2r_blastp/${qaccession}.q2r.blastp
#           bdbh_base_dir/r2q_blastp/${qaccession}.r2q.blastp
parser.add_argument("--bdbh_matrix", dest="to_get_bdbh_matrix", required=False, help="collect bdbh hit matrix", default = False, action="store_true")
#   from running this command you will create under bdbh_base_dir the followings:
#           bdbh_base_dir/bdbh_map
#           bdbh_base_dir/bdbh_map/${qaccession}.bdbh  -- a file per query genome -- which contains bdbh mapping from ref protein to query protein 
#           bdbh_base_dir/bdbh_map.num_core_per_freq.tsv    -- a file containing the table of core genome freqs (1, 0.99, 0.98, ... 0.90) vs. the number of core genes
#           bdbh_base_dir/bdbh_map.ref_gene_genome_freq.tsv    -- a file containing the table of reference protin vs. its frequency in the query genomes
parser.add_argument("--core_seq", dest="to_extract_core_seq", required=False, help="collect core gene sequence fasta per core gene family", default = False, action="store_true")
#   from running this command you will create under bdbh_base_dir the followings:
#           bdbh_base_dir/bdbh_family_${moltype}
#           bdbh_base_dir/bdbh_family_${moltype}/${rprotein}.fasta
parser.add_argument("--bdbh_base_dir", dest="bdbh_base_dir", help="(always required) the output base directory; under which blastp files will be arranged", required=False, type=str, default="-na")
parser.add_argument("--num_thread", dest="num_thread", help="(optional for --run_bdbh; default is 1) num threads to use in diamond blastp/makedb", required=False, type=str, default="1")
parser.add_argument("--ref_faa", dest="ref_faa_file", help="(required for --run_bdbh and --bdbh_matrix) the proteome faa file of refernece genome", required=False, type=str, default="-na")
parser.add_argument("--query_faa_list", dest="query_faa_file_list", help="(required for --run_bdbh) a file containing the list of query proteom faa file paths (1 per line)", required=False, type=str, default="-na")
parser.add_argument("--faa_ext", dest="query_faa_extension", help="(required for --run_bdbh and --core_seq) the terminal extension of query feature seq files; this is to infer the strain name from fasta file path", required=False, type=str, default="-na")
parser.add_argument("--query_acc_list", dest="query_acc_file_list", help="(required for --bdbh_matrix) a file containing the list of query genome accessions (1 per line)", required=False, type=str, default="-na")
parser.add_argument("--query_fseq_list", dest="query_fseq_file_list", help="(required for --core_seq) a file containing the list of query genome feature fasta (whether prot or dna upon your choice) (1 per line)", required=False, type=str, default="-na")
parser.add_argument("--moltype", dest="feature_mol_type", help="(required for --core_seq) (prot or dna) prot or dna sequence of features, are you collecting", required=False, type=str, default="-na")
parser.add_argument("--core_freq", dest="core_freq_cutoff", help="(default 1) (float between 0 - 1) core gene's minimal presence frequency among query genomes", required=False, type=float, default=1)

args = parser.parse_args()

to_run_bdbh = args.to_run_bdbh
to_get_bdbh_matrix = args.to_get_bdbh_matrix
to_extract_core_seq = args.to_extract_core_seq
ref_faa_file = args.ref_faa_file
num_thread = args.num_thread
query_faa_file_list = args.query_faa_file_list
query_faa_extension = args.query_faa_extension
bdbh_base_dir = args.bdbh_base_dir
query_acc_file_list = args.query_acc_file_list
query_fseq_file_list = args.query_fseq_file_list
feature_mol_type = args.feature_mol_type
core_freq_cutoff = args.core_freq_cutoff

if to_run_bdbh:
    if ref_faa_file.startswith("-na"):
        print("since you are doing --run_bdbh you must give input to --ref_faa")
    if query_faa_file_list.startswith("-na"):
        print("since you are doing --run_bdbh you must give input to --query_faa_list")
    if query_faa_extension.startswith("-na"):
        print("since you are doing --run_bdbh you must give input to --faa_ext")
    if bdbh_base_dir.startswith("-na"):
        print("since you are doing --run_bdbh you must give input to --bdbh_base_dir")

if to_get_bdbh_matrix:
    if ref_faa_file.startswith("-na"):
        print("since you are doing --bdbh_matrix you must give input to --ref_faa")
    if query_acc_file_list.startswith("-na"):
        print("since you are doing --bdbh_matrix you must give input to --query_acc_list")
    if bdbh_base_dir.startswith("-na"):
        print("since you are doing --bdbh_matrix you must give input to --bdbh_base_dir")

if to_extract_core_seq:
    if bdbh_base_dir.startswith("-na"):
        print("since you are doing --core_seq you must give input to --bdbh_base_dir")
    if query_fseq_file_list.startswith("-na"):
        print("since you are doing --core_seq you must give input to --query_fseq_list")
    if feature_mol_type.startswith("-na"):
        print("since you are doing --core_seq you must give input to --moltype")
    if query_faa_extension.startswith("-na"):
        print("since you are doing --core_seq you must give input to --faa_ext")

if to_run_bdbh:
    # inputs to use
    #       bdbh_base_dir
    #       ref_faa_file
    #       query_faa_file_list
    #       query_faa_extension
    #you will create under bdbh_base_dir the followings:
    #           bdbh_base_dir/q2r_blastp
    #           bdbh_base_dir/r2q_blastp
    #           bdbh_base_dir/q2r_blastp/${qaccession}.q2r.blastp
    #           bdbh_base_dir/r2q_blastp/${qaccession}.r2q.blastp
    
    # prep dir for output bdbh map per query genome
    if not os.path.exists(bdbh_base_dir):
        os.makedirs(bdbh_base_dir)
    if not os.path.exists(bdbh_base_dir + "/q2r_blastp"):
        os.makedirs(bdbh_base_dir + "/q2r_blastp")
    if not os.path.exists(bdbh_base_dir + "/r2q_blastp"):
        os.makedirs(bdbh_base_dir + "/r2q_blastp")
    
    # reference genome's faa -> make blastp db
    ref_dmnd_file = bdbh_base_dir + "/q2r_blastp/reference.dmnd"
    ex_command = ["diamond", "makedb", "--in", ref_faa_file, "-p", num_thread, "-d", ref_dmnd_file]
    subprocess.run(ex_command)

    # take care of each query faa
    faa_list_read = open(query_faa_file_list, "r")
    for line in faa_list_read:
        query_faa_file = line.strip()
        path_split = query_faa_file.split('/')
        terminal_path = path_split[len(path_split) - 1]
        query_accession = terminal_path[0:(len(terminal_path) - len(query_faa_extension))]
        
        # query faa -> make blastp db
        query_dmnd_file = bdbh_base_dir + "/r2q_blastp/" + query_accession + ".dmnd"
        ex_command = ["diamond", "makedb", "--in", query_faa_file, "-p", num_thread, "-d", query_dmnd_file]
        subprocess.run(ex_command)
        # q2r blastp
        q2r_blastp_output = bdbh_base_dir + "/q2r_blastp/" + query_accession + ".q2r.blastp"
        ex_command = ["diamond", "blastp", "-query", query_faa_file, "-p", num_thread, "-d", ref_dmnd_file, "-f", "6", "-k", "1", "-e", "1e-12", "-o", q2r_blastp_output]
        subprocess.run(ex_command)
        # r2q blastp
        r2q_blastp_output = bdbh_base_dir + "/r2q_blastp/" + query_accession + ".r2q.blastp"
        ex_command = ["diamond", "blastp", "-query", ref_faa_file, "-p", num_thread, "-d", query_dmnd_file, "-f", "6", "-k", "1", "-e", "1e-12", "-o", r2q_blastp_output]
        subprocess.run(ex_command)
        
        # remove query db file
        os.remove(query_dmnd_file)
    faa_list_read.close()

    # remove reference dmnd
    os.remove(ref_dmnd_file)


if to_get_bdbh_matrix:
    #inputs to use
    #   ref_faa_file
    #   query_acc_file_ilst
    #           bdbh_base_dir/q2r_blastp/${qaccession}.q2r.blastp
    #           bdbh_base_dir/r2q_blastp/${qaccession}.r2q.blastp
    #you will create under bdbh_base_dir the followings:
    #           bdbh_base_dir/bdbh_map
    #           bdbh_base_dir/bdbh_map/${qaccession}.bdbh  -- a file per query genome -- which contains bdbh mapping from ref protein to query protein 
    #           bdbh_base_dir/bdbh_map.num_core_per_freq.tsv    -- a file containing the table of core genome freqs (1, 0.99, 0.98, ... 0.90) vs. the number of core genes
    #           bdbh_base_dir/bdbh_map.ref_gene_genome_freq.tsv    -- a file containing the table of reference protin vs. its frequency in the query genomes
    
    # reference proteins
    dict_ref_protein_rpindex = {}
    dict_rpindex_ref_protein = {}
    rpindex_cursor = 0
    rp_read = open(ref_faa_file, "r")
    for line in rp_read:
        line_tr = line.strip()
        if(line_tr.startswith('>')):
            ref_protein = line_tr[1:]
            dict_ref_protein_rpindex[ref_protein] = rpindex_cursor
            dict_rpindex_ref_protein[rpindex_cursor] = ref_protein
            rpindex_cursor += 1
    rp_read.close()
    num_ref_protein = len(dict_rpindex_ref_protein)

    # genome counter per ref protein
    list_num_genome_per_rpindex = [0]*num_ref_protein

    # query genomes
    dict_query_accession_qgindex = {}
    dict_qgindex_query_accession = {}
    qgindex_cursor = 0
    qa_read = open(query_acc_file_list, "r")
    for line in qa_read:
        query_accession = line.strip()
        dict_query_accession_qgindex[query_accession] = qgindex_cursor
        dict_qgindex_query_accession[qgindex_cursor] = query_accession
        qgindex_cursor += 1
    qa_read.close()
    num_query_genome = len(dict_qgindex_query_accession)

    # prep dir for output bdbh map per query genome
    if not os.path.exists(bdbh_base_dir + "/bdbh_map"):
        os.makedirs(bdbh_base_dir + "/bdbh_map")

    # read each query's blastp files doing: (1) update the genome counter (per ref protein) and (2) write bdbh map for each query
    for qgindex in range(num_query_genome):
        query_accession = dict_qgindex_query_accession[qgindex]
        q2r_blastp_file = bdbh_base_dir + "/q2r_blastp/" + query_accession + ".q2r.blastp"
        r2q_blastp_file = bdbh_base_dir + "/r2q_blastp/" + query_accession + ".r2q.blastp"
        dict_ref_query_hit = blastout_to_dict_query_ref(r2q_blastp_file)
        dict_query_ref_hit = blastout_to_dict_query_ref(q2r_blastp_file)

        query_bdbh_map = bdbh_base_dir + "/bdbh_map/" + query_accession + ".bdbh"
        bdbh_map_w = open(query_bdbh_map, "w")
        for rpindex in range(num_ref_protein):
            ref_protein = dict_rpindex_ref_protein[rpindex]
            query_hit = ""
            bdbh_present = True
            if ref_protein in dict_ref_query_hit:
                query_hit = dict_ref_query_hit[ref_protein]
                if query_hit in dict_query_ref_hit:
                    reciprocal_ref = dict_query_ref_hit[query_hit]
                    if reciprocal_ref != ref_protein:
                        bdbh_present = False
                else:
                    bdbh_present = False
            else:
                bdbh_present = False

            if bdbh_present:
                bdbh_map_w.write(ref_protein + "\t" + query_hit + "\n")
                list_num_genome_per_rpindex[rpindex] += 1
        bdbh_map_w.close()
    
    # time to write
    #   - bdbh_base_dir/bdbh_map.ref_gene_genome_freq.tsv    -- a file containing the table of reference protin vs. its frequency in the query genomes
    genome_freq_file = bdbh_base_dir + "/bdbh_map.ref_protein_genome_freq.tsv"
    list_freq_genome_per_rpindex = []
    genome_freq_w = open(genome_freq_file, "w")
    for rpindex in range(num_ref_protein):
        ref_protein = dict_rpindex_ref_protein[rpindex]
        num_genome = list_num_genome_per_rpindex[rpindex]
        freq_genome = float(num_genome)/float(num_query_genome)
        list_freq_genome_per_rpindex.append(freq_genome)
        genome_freq_w.write(ref_protein + "\t" + str(num_genome) + "\t" + str(freq_genome) + "\n")
    genome_freq_w.close()

    # time to count these and write
    #   - bdbh_base_dir/bdbh_map.num_core_per_freq.tsv    -- a file containing the table of core genome freqs (1, 0.99, 0.98, ... 0.90) vs. the number of core genes
    core_freqs = [0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0]
    num_core_gene = [0]*len(core_freqs)
    for cutindex in range(len(core_freqs)):
        core_freq_cut = core_freqs[cutindex]
        num_core_gene_at_cutindex = 0
        for rpindex in range(num_ref_protein):
            freq_genome = list_freq_genome_per_rpindex[rpindex]
            if freq_genome >= core_freq_cut:
                num_core_gene_at_cutindex += 1
        num_core_gene[cutindex] = num_core_gene_at_cutindex
    count_per_cutoff_file = bdbh_base_dir + "/bdbh_map.num_core_per_freq.tsv"
    count_per_cutoff_w = open(count_per_cutoff_file, "w")
    for cutindex in range(len(core_freqs)):
        cutoff = core_freqs[cutindex]
        count = num_core_gene[cutindex]
        count_per_cutoff_w.write(str(cutoff) + "\t" + str(count) + "\n")
    count_per_cutoff_w.close()



if to_extract_core_seq:
    # inputs to be used
    #           bdbh_base_dir 
    #           query_fseq_file_list
    #           feature_mol_type
    #           query_faa_extension
    #           core_freq_cutoff
    # from running this command you will create under bdbh_base_dir the followings:
    #           bdbh_base_dir/bdbh_family_${moltype}
    #           bdbh_base_dir/bdbh_family_${moltype}/${rprotein}.fasta
    ref_gene_genome_freq_file = bdbh_base_dir + "/bdbh_map.ref_protein_genome_freq.tsv"
    list_core_gene_ref_protein = list_up_core_gene_ref_protein(ref_gene_genome_freq_file, core_freq_cutoff)
    num_core_gene = len(list_core_gene_ref_protein)

    # prepare output (sequence fasta files) directory if not exist
    if not os.path.exists(bdbh_base_dir + "/bdbh_family_" + feature_mol_type):
        os.makedirs(bdbh_base_dir + "/bdbh_family_" + feature_mol_type)

    # file to write in the strain (accession) list for which the feature sequences collected
    acc_list_file = bdbh_base_dir + "/bdbh_family_" + feature_mol_type + ".enlisted_accessions"

    # delete any file pre-existing in the fasta to be writen per core gene
    for ref_protein in list_core_gene_ref_protein:
        seq_output = bdbh_base_dir + "/bdbh_family_" + feature_mol_type + "/" + ref_protein + ".fasta"
        if os.path.exists(seq_output):
            os.remove(seq_output)
    
    # add core seq from all query strains to the per-core-gene fasta files using appending writers
    acc_list_writer = open(acc_list_file, "w")
    fseq_list_read = open(query_fseq_file_list, "r")
    for line in fseq_list_read:
        fseq_file = line.strip()
        fseq_file_path_split = fseq_file.split('/')
        fseq_file_endpath = fseq_file_path_split[len(fseq_file_path_split) - 1]
        query_accession = fseq_file_endpath[0:(len(fseq_file_endpath) - len(query_faa_extension))]
        acc_list_writer.write(query_accession + "\n")

        bdbh_map_file = bdbh_base_dir + "/bdbh_map/" + query_accession + ".bdbh"
        dict_query_protein_ref_protein = blastout_to_dict_ref_query(bdbh_map_file)

        fseq_read = open(fseq_file, "r")
        line = fseq_read.readline()
        while line != '':
            if line.startswith('>'):
                seqid = line.strip()[1:].split(" ")[0]
                line = fseq_read.readline()
                seq = []
                while line != '':
                    seq.append(line.strip())
                    line = fseq_read.readline()
                    if line.startswith('>'):
                        break
                # query protein -> ref protein bdbh -> is it core?
                if seqid in dict_query_protein_ref_protein:
                    bdbh_ref_protein = dict_query_protein_ref_protein[seqid]
                    if bdbh_ref_protein in list_core_gene_ref_protein:
                        # yes it is core gene family's ORF
                        seq_output = bdbh_base_dir + "/bdbh_family_" + feature_mol_type + "/" + bdbh_ref_protein + ".fasta"
                        seq_write = open(seq_output, "a")
                        seq_write.write(">" + query_accession + "\n" + "".join(seq) + "\n")
                        seq_write.close()
            else:
                line = fseq_read.readline()
        fseq_read.close()

    fseq_list_read.close()
    acc_list_writer.close()

