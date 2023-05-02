from Bio import SeqIO
import argparse

def get_uniprot_id(input):
    protein_id =[]
    for seq_record in SeqIO.parse(input, "genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                protein_id.append(seq_feature.qualifiers['protein_id'][0])
    return protein_id

def parse_arguments():

    parser = argparse.ArgumentParser(description='Extract protein_id from genbank')
    parser.add_argument("-gb", '--genbank', default=None, help="Input should be a genbank file")

    return parser


def main():

    parser=parse_arguments()
    args=parser.parse_args()

    protein_id = get_uniprot_id(args.genbank)
    for protein in protein_id:
        print (protein)

    return 0


if __name__=='__main__':
    main()
