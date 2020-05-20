# $1 = input fasta
# $2 = output fasta
input_fasta=$1
output_fasta=$2
python DereplicateRereplicateSeqFasta.py --i ${input_fasta} --o ${input_fasta}.streamline_aln.derep.fa --repmap ${input_fasta}.streamline_aln.derep.map --do_derep
muscle -in ${input_fasta}.streamline_aln.derep.fa -out ${input_fasta}.streamline_aln.derep.aln
python DereplicateRereplicateSeqFasta.py --i ${input_fasta}.streamline_aln.derep.aln --o ${output_fasta} --repmap ${input_fasta}.streamline_aln.derep.map --do_rerep
rm ${input_fasta}.streamline_aln.derep.map
rm ${input_fasta}.streamline_aln.derep.fa
rm ${input_fasta}.streamline_aln.derep.aln

