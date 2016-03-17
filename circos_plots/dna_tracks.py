import csv
import argparse

from pprint import pprint

parser = argparse.ArgumentParser(description='obtiene frecuencias alelicas y tipos de variacion somatica para tablas de TCGA.')
parser.add_argument('--tabla', type=argparse.FileType('r'), required=True, help="tabla de TCGA")
parser.add_argument('--prefix', required=True, help="prefix for output files")
args = parser.parse_args()

reader = csv.reader(args.tabla, delimiter="\t")

header = reader.next()

with open("%s_var_types_tiles.tsv" % args.prefix, 'w') as vartypes:
    vartypes_writer = csv.writer( vartypes, delimiter=" ")

    for v in reader:
        (Hugo_Symbol, Entrez_Gene_Id, Center, Ncbi_Build, Chrom, Start_Position, End_Position, Strand,
         Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
         Dbsnp_Rs, Dbsnp_Val_Status, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode,
         Match_Norm_Seq_Allele1, Match_Norm_Seq_Allele2, Tumor_Validation_Allele1, Tumor_Validation_Allele2,
         Match_Norm_Validation_Allele1, Match_Norm_Validation_Allele2, Verification_Status, Validation_Status,
         Mutation_Status, Sequencing_Phase, Sequence_Source, Validation_Method, Score, Bam_File, Sequencer,
         Tumor_Sample_UUID, Matched_Norm_Sample_UUID, File_Name, Archive_Name, Line_Number) = v

        vartypes_writer.writerow(["hs%s" % Chrom, Start_Position, End_Position, "id=%s" % Variant_Type])
    
