import logging
import time
import gzip
import os
import pandas as pd
from joblib.parallel import Parallel, delayed

from dev_settings import LOGGER_NAME
from CustomLogger import initLogging

def process_line(line):
    """Processes a single line into a dictionary data object
    Expecting a line to be formatted as
    CHR START_POS STOP_POS REF;ALT;GENE_SYMBOL;VAR_TYPE
        OR
    CHR START_POS STOP_POS REF;ALT;GENE_SYMBOL;VAR_TYPE;META_DOMAIN_INFO
    """

    # decode if line is gzipped
    if type(line) is bytes:
        line = line.decode('utf-8')

    # split line into tokens
    processed_line = line.strip().split('\t')
    info_tokens = processed_line[3].split(';')

    # craft new data entry
    data_entry = {}
    data_entry['chrom'] = processed_line[0]
    data_entry['pos'] = processed_line[1]
    data_entry['ID'] = info_tokens[0]
    data_entry['ref'] = info_tokens[1]
    data_entry['alt'] = info_tokens[2]
    data_entry['symbol'] = info_tokens[3]
    data_entry['consequence'] = info_tokens[4]

    # Check if there is meta-domain annotation
    if len(info_tokens) == 5:
        data_entry['metadomain_info'] = '-'
    elif len(info_tokens) == 6:
        data_entry['metadomain_info'] = info_tokens[5]
    elif len(info_tokens) > 6:
        data_entry['metadomain_info'] = ",".join(info_tokens[5:])
    else:
        logging.getLogger(LOGGER_NAME).warning("Not enough info tokens found than expected for line '"+str(line)+"'")

    # return data_entry
    return data_entry


def read_and_tokenize_bed_file_line(file_reader, line_counter, report_name):
    """Decorator function that reads and tokenizes a single bed file line with logging service support
    """

    # Read a line of the lifted bed file
    line = file_reader.readline()
    # Check if the end of the file has been reached
    end_of_file = not line

    if not end_of_file:
        # Increment the unlifted line counter
        line_counter += 1
        # Tokenize the line
        processed_line = process_line(line)
    else:
        processed_line = None
        logging.getLogger(LOGGER_NAME).info(
            "Reached end of "+str(report_name)+" file at line " + str(line_counter))

    return line, processed_line, line_counter, end_of_file

def process_line_to_write(processed_line):
    """ Processes a dictionary object of a line to a writeable string format

    Input is a previously processed line from method 'process_line()'
    Output:

    CHR <tab> START_POS <tab> REF <tab> ALT <tab> GENE_SYMBOL <space> VAR_TYPE <space> LIFTED_CHR_POS <space> META_DOMAIN_INFO
    """

    # construct the output string
    output_string = ""
    output_string += str(processed_line['chrom']) + "\t"
    output_string += str(processed_line['pos']) + "\t"
    output_string += str(processed_line['ref']) + "\t"
    output_string += str(processed_line['alt']) + "\t"
    output_string += str(processed_line['symbol']) + " "
    output_string += str(processed_line['consequence'])

    if 'lifted_pos' in processed_line.keys():
        output_string += " " + str(processed_line['lifted_pos'])
        output_string += " " + str(processed_line['metadomain_info'])

    output_string += "\n"

    return output_string



def combine_liftover(unlifted_bed_filename, lifted_bed_filename, output_filename, report_counter=100000):
    """


    # Combine liftover [GChr 38 with ID] +  [hg19 with ID] >>> [GChr38 with ID+hg 19 loc]
    """

    ## unlifted_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz:
    ## chr10   10578   10578   1;G;A;ENSG00000260370;downstream_gene_variant
    ## lifted_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz:
    ## chr18   10774   10774   1;G;A;ENSG00000260370;downstream_gene_variant

    # Read the unlifted bed file
    with gzip.open(unlifted_bed_filename, 'r') as unlifted_f, gzip.open(lifted_bed_filename, 'r') as lifted_f, gzip.open(output_filename, 'wb') as output_f:
        processed_line_unlifted_counter = 0
        processed_line_lifted_counter = 0

        # Read the first line of the unlifted bed file
        _, processed_line_unlifted, processed_line_unlifted_counter, unlifted_end_of_file = read_and_tokenize_bed_file_line(unlifted_f, processed_line_unlifted_counter, "unlifted")
        # Read the first line of the lifted bed file
        _, processed_line_lifted, processed_line_lifted_counter, lifted_end_of_file = read_and_tokenize_bed_file_line(lifted_f, processed_line_lifted_counter, "lifted")

        # Compare the two lines and merge them
        while not unlifted_end_of_file or not lifted_end_of_file:
            # Report the progress
            if processed_line_unlifted_counter % report_counter == 0:
                logging.getLogger(LOGGER_NAME).info("Processed '"+str(processed_line_lifted_counter)+"' lifted lines and merged onto '"+str(processed_line_unlifted_counter)+"' lines")

            if lifted_end_of_file or int(processed_line_unlifted['ID']) <= int(processed_line_lifted['ID']):
                if lifted_end_of_file or int(processed_line_unlifted['ID']) < int(processed_line_lifted['ID']):
                    # Write the unlifted line
                    output_f.write(process_line_to_write(processed_line_unlifted).encode("utf-8"))

                    # Read the next line of the unlifted bed file
                    _, processed_line_unlifted, processed_line_unlifted_counter, unlifted_end_of_file = read_and_tokenize_bed_file_line(
                        unlifted_f, processed_line_unlifted_counter, "unlifted")
                    continue
                elif int(processed_line_unlifted['ID']) == int(processed_line_lifted['ID']):
                    # add the info from the lifted line to the unlifted line
                    processed_line_unlifted['lifted_pos'] = processed_line_lifted['chrom'] + ":" + processed_line_lifted['pos']
                    processed_line_unlifted['metadomain_info'] = processed_line_lifted['metadomain_info']

                    # Write the unlifted line
                    output_f.write(process_line_to_write(processed_line_unlifted).encode("utf-8"))

                    # Read the next lines
                    _, processed_line_unlifted, processed_line_unlifted_counter, unlifted_end_of_file = read_and_tokenize_bed_file_line(
                        unlifted_f, processed_line_unlifted_counter, "unlifted")
                    _, processed_line_lifted, processed_line_lifted_counter, lifted_end_of_file = read_and_tokenize_bed_file_line(
                        lifted_f, processed_line_lifted_counter, "lifted")
                    continue
            else:
                logging.getLogger(LOGGER_NAME).error(
                    "asynchonous error: unlifted (line:" + str(processed_line_unlifted_counter) + ") ID:" + str(
                        processed_line_unlifted['ID']) + " > lifter (line:" + str(
                        processed_line_lifted_counter) + ") ID:" + str(processed_line_lifted['ID']))
                raise Exception("Error: liftover failed")

        logging.getLogger(LOGGER_NAME).info("Finished merging liftover with unlifted bed file. Processed '"+str(processed_line_lifted_counter)+"' lifted lines and merged onto '"+str(processed_line_unlifted_counter)+"' lines")


def chunk_the_chromosomes(source_filename, target_folder, target_filename, report_counter=100000):
    """

    # Chunk per chromosome [GChr38 with ID+hg 19 loc] (split_per_chrom) >>>  [ChrX][GChr38 with ID+hg 19 loc]

    """

    # read the source file
    logging.getLogger(LOGGER_NAME).info("Starting chunking the unlifted bed file into seperate chromosome files")

    # Create the target folder if it doesn't exist
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
    else:
        # Remove the existing files
        for file in os.listdir(target_folder):
            os.remove(os.path.join(target_folder, file))

    # Iterate over the source file
    processed_line_counter = 0
    with gzip.open(source_filename, 'r') as source_f:
        # Read the first line of the source file
        line, processed_line, processed_line_counter, end_of_file = read_and_tokenize_bed_file_line(source_f, processed_line_counter, "chunk_by_chromosome")

        while not end_of_file:
            # Compose the target filename
            target_file = target_folder + processed_line['chrom'] + "_" + target_filename
            # Set the new chromosome the next line is on
            current_chrom = processed_line['chrom']

            logging.getLogger(LOGGER_NAME).info("Writing to '"+target_file+"'")

            # Append the line to the target file
            with open(target_file, 'a') as target_f:
                while not end_of_file and current_chrom == processed_line['chrom']:
                    if processed_line_counter % report_counter == 0:
                        logging.getLogger(LOGGER_NAME).info("Processed '"+str(processed_line_counter)+"' lines")

                    # Append the line to the target file
                    target_f.write(line.decode("utf-8"))
                    line, processed_line, processed_line_counter, end_of_file = read_and_tokenize_bed_file_line(source_f,
                                                                                                          processed_line_counter,
                                                                                                            "chunk_by_chromosome")

    logging.getLogger(LOGGER_NAME).info("Finished chunking the unlifted bed file into seperate chromosome files")

def annotate_metadome(metadome_df, chromosome, chr_pos):
    """

    """

    # Make sure the metadome_df provided contains data from only one chromosome
    chromosomes_in_metadome_df = pd.unique(metadome_df.chrom)

    return_value = {}

    if not chromosome in chromosomes_in_metadome_df:
        logging.getLogger(LOGGER_NAME).error("The metadome_df provided does not contain data from the chromosome '"+chromosome+"'")
        return return_value

    assert chromosome in chromosomes_in_metadome_df, "Error: the metadome_df provided does not contain data from the chromosome '"+chromosome+"'"
    assert len(chromosomes_in_metadome_df) == 1, "Error: the metadome_df provided contains data from multiple chromosomes"

    # check for a metadome hit
    metadome_hit = metadome_df[(metadome_df.pos_start <= int(chr_pos)) & (int(chr_pos) <= metadome_df.pos_stop)]

    # if there is a hit, preprocess the hit
    if len(metadome_hit) > 0:
        # add record
        records = metadome_hit.to_dict('records')
        return_value['metadome_symbol'] = ";".join(set([x['symbol'] for x in records]))
        return_value['metadome_metadomain_positions'] = ";".join(
            set([str(x['domain_id']) + ":" + str(int(x['consensus_pos'])) for x in records]))

    # return the hit
    return return_value

def annotate_chromosome_chunk(chromosome, chunk_filename, metadome_df_chunk, target_folder, target_filename, report_counter):
    """

    """

    logging.getLogger(LOGGER_NAME).info("Starting annotation of MetaDome data to '"+chunk_filename+"'")

    if len(metadome_df_chunk) == 0:
        logging.getLogger(LOGGER_NAME).info("No metadome data for '"+chromosome+"'")

    # Initialize counter
    processed_line_counter = 0
    # Open the source file
    with open(chunk_filename, 'r') as source_f, open(os.path.join(target_folder, chromosome + "_" + target_filename), 'a') as target_f:
        # read the first line
        line, processed_line, processed_line_counter, end_of_file = read_and_tokenize_bed_file_line(source_f, processed_line_counter, "annotate_chunk_" + chromosome)

        while not end_of_file:
            if processed_line_counter % report_counter == 0:
                logging.getLogger(LOGGER_NAME).info("Processed '" + str(processed_line_counter) + "' lines for chromosome '" + chromosome + "'")

            processed_line['metadome_metadomain_positions'] = "-"
            if metadome_df_chunk is not None:
                annotation = annotate_metadome(metadome_df_chunk, processed_line['chrom'], processed_line['pos'])
                if len(annotation) > 0:
                    processed_line['metadome_metadomain_positions'] = annotation['metadome_metadomain_positions']

            # write the line to the target file
            if type(line) is bytes:
                line_to_write = line.decode('utf-8').strip()
            else:
                line_to_write = line.strip()
            line_to_write += ";" + processed_line['metadome_metadomain_positions'] + "\n"

            # write the line to the target file
            target_f.write(line_to_write)

            # read the next line
            line, processed_line, processed_line_counter, end_of_file = read_and_tokenize_bed_file_line(source_f, processed_line_counter, "annotate_chunk_" + chromosome)


def annotate_the_chromosomes(source_folder, target_folder, target_filename, metadome_filename, parallel = True, report_counter=10000):
    """

    # Annotate per chromosome (parallel) [ChrX][GChr38 with ID+hg 19 loc] (domain annotate) >>> [ChrX][GChr38 with ID+hg 19 loc + domain annotation]
    """
    # read the source file
    logging.getLogger(LOGGER_NAME).info("Starting annotation of the chunked chromosome files from '"+source_folder+"'")

    # Create the target folder if it doesn't exist
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
    else:
        # Remove the existing files
        for file in os.listdir(target_folder):
            os.remove(os.path.join(target_folder, file))

    # Check the files in the source folder
    source_files = {}
    for file in os.listdir(source_folder):
        if file.endswith(".bed") and file.startswith("chr"):
            # retrieve the specific chromosome
            chromosome = file.split("_")[0]
            # add the file to the dictionary
            source_files[chromosome] = os.path.join(source_folder, file)

            logging.getLogger(LOGGER_NAME).info("Found chromosome '"+chromosome+"' in '"+source_folder+"' added to the list of files to annotate: '"+source_folder+file+"'")
        else:
            logging.getLogger(LOGGER_NAME).warning("Skipping file '"+file+"' in source folder '"+source_folder+"'")


    # Read the metadome annotation
    metadome_df = pd.read_csv(metadome_filename, index_col=False)

    # Select to only positions that have domain annotation
    metadome_df = metadome_df[(metadome_df.domain_id == metadome_df.domain_id)]

    # Drop duplicates, if any
    metadome_df = metadome_df[['chrom', 'pos_start', 'pos_stop', 'strand', 'symbol', 'domain_id', 'consensus_pos']].drop_duplicates()

    if parallel:
        _ = Parallel(n_jobs=8)(delayed(annotate_chromosome_chunk)(chromosome, source_files[chromosome], metadome_df[metadome_df.chrom == chromosome], target_folder, target_filename, report_counter) for chromosome in source_files.keys())
    else:
        for chromosome in source_files.keys():
            annotate_chromosome_chunk(chromosome, source_files[chromosome], metadome_df[metadome_df.chrom == chromosome], target_folder, target_filename, report_counter)


def combine_annotated_chunks(source_folder, target_filename, report_counter=100000):
    """

    # Combine annotated chunks [ChrX][GChr38 with ID+hg 19 loc + domain annotation] (combine) >>> [GChr38 with ID+hg 19 loc + domain annotation]
    """

    logging.getLogger(LOGGER_NAME).info("Starting combining of the annotated chromosome files from '"+source_folder+"'")

    # Check the files in the source folder
    source_files = {}
    for file in os.listdir(source_folder):
        if file.endswith(".bed") and file.startswith("chr"):
            # retrieve the specific chromosome
            chromosome = file.split("_")[0]
            # add the file to the dictionary
            source_files[chromosome] = os.path.join(source_folder, file)

            logging.getLogger(LOGGER_NAME).info(
                "Found chromosome '" + chromosome + "' in '" + source_folder + "' added to the list of files to annotate: '" + source_folder + file + "'")
        else:
            logging.getLogger(LOGGER_NAME).warning(
                "Skipping file '" + file + "' in source folder '" + source_folder + "'")


    # Open the target file
    processed_line_counter = 0
    with gzip.open(target_filename, 'wb') as target_f:
        for chromosome in source_files.keys():
            # Open the source file
            with open(source_files[chromosome], 'r') as source_f:
                # read the first line
                line, processed_line, processed_line_counter, end_of_file = read_and_tokenize_bed_file_line(source_f, processed_line_counter, "combine_chunk_" + chromosome)

                while not end_of_file:
                    if processed_line_counter % report_counter == 0:
                        logging.getLogger(LOGGER_NAME).info("Written '" + str(processed_line_counter) + "' lines so far and currently at chromosome '" + chromosome + "'")

                    # write the line to the target file
                    if type(line) is bytes:
                        line_to_write = line.decode('utf-8').strip()
                    else:
                        line_to_write = line.strip()
                    line_to_write += "\n"

                    # write the line to the target file
                    target_f.write(line_to_write.encode('utf-8'))

                    # read the next line
                    line, processed_line, processed_line_counter, end_of_file = read_and_tokenize_bed_file_line(source_f, processed_line_counter, "combine_chunk_" + chromosome)

    logging.getLogger(LOGGER_NAME).info("Finished combining of the annotated chromosome files from '"+source_folder+"' and wrote '"+str(processed_line_counter)+"' lines to the target file '"+target_filename+"'")




# Convert to .bed and Add ID [GChr 38] -> [GChr 38 with ID]
## e.g.: genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz | awk '{$3=$3";"$4";"$5";"$6;$6="";$5="";$4="";}1' | sed 's/ /\t/g' | sed -e 's/ \+/\t/g' |  awk -F"\t" '{$2=$2"\t"$2; printf $1"\t"; for(i=2;i<NF;i++){printf $i"\t"}; print ""}' | awk '{$1=$1};1' | awk '/./ {$4 =   NR ";" $4}; {print}' |  sed -e 's/ /\t/g' | gzip > preliftover_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz

# Liftover [GChr 38 with ID] (Liftover)>>> [hg19 with ID]
## e.g.: ./liftOver preliftover_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz hg38ToHg19.over.chain.gz liftover_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz unlifted_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz
## ex: chr10	10578	10578	1;A;G;ENSG00000XXXXXX;variant_type


# Combine liftover [GChr 38 with ID] +  [hg19 with ID] >>> [GChr38 with ID+hg 19 loc]

# Annotate per chromosome (parallel) [ChrX][GChr38 with ID+hg 19 loc] (domain annotate) >>> [ChrX][GChr38 with ID+hg 19 loc + domain annotation]

# Combine annotated chunks [ChrX][GChr38 with ID+hg 19 loc + domain annotation] (combine) >>> [GChr38 with ID+hg 19 loc + domain annotation]



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Annotates a  the spatial clustering of variant locations over cDNA')

    # required arguments
    parser.add_argument('--gene_name', type=str, required=True,
                        help='(Required) Name of the gene of interest, example usage: --gene_name=BRCA1')
    parser.add_argument('--variant_cDNA_locations', type=str, required=True,
                        help='(Required) cDNA based variant locations, example usage: --variant_cDNA_locations=10,50,50,123')
    parser.add_argument('--cDNA_length', type=int, required=True,
                        help='(Required) total cDNA length of the gene (including stop codon), example usage: --cDNA_length=1337')

    # optional arguments
    parser.add_argument('--n_permutations', type=int, required=False, default=100000000,
                        help='(Optional) total nunber of permutations, default=100000000 (1.00E+08), example usage: --n_permutations=100')
    parser.add_argument('--parallel', type=bool, required=False, default=True,
                        help='(Optional) should the algorithm make use of parallel computation?, default=True, example usage: --parallel=True')
    parser.add_argument('--random_seed', type=int, required=False, default=1,
                        help='(Optional) The seed used for initialization of the random permutations, default=1, example usage: --random_seed=1')
    parser.add_argument('--correction', type=int, required=False, default=1,
                        help='(Optional) The number of genes the p-value must be corrected for in a Bonferonni manner, default=1, example usage: --correction=1')

    args = parser.parse_args()
    main(gene_name=args.gene_name, variant_cDNA_locations=args.variant_cDNA_locations.split(','),
         cDNA_length=args.cDNA_length, n_permutations=args.n_permutations, parallel=args.parallel,
         random_seed=args.random_seed, correction=args.correction)



    initLogging(print_to_console=True, logging_level=logging.INFO)


    logging.getLogger(LOGGER_NAME).info("Starting annotation")
    start_time = time.perf_counter()
    

    _metadome_filename = "/usr/data/metadome_data_full_n_transcripts_41772_20220508-011127.tsv.gz"
    # _input_data_filename =  "/usr/data/results/liftover_MultiomicWatershedUpdated-Prediction-AllData_summarizedPosterior-FilteredBygnomAD.tsv"
    # _output_data_filename = "/usr/data/results/metadome_annotated_liftover_MultiomicWatershedUpdated-Prediction-AllData_summarizedPosterior-FilteredBygnomAD.tsv"
    #
    # annotate_dataset_with_metadome(_metadome_filename, _input_data_filename, _output_data_filename, parallel=True)

    _unlifted_input_data_filename = "/usr/data/preliftover_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz"
    _lifted_input_data_filename = "/usr/data/liftover_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz"


    # Chunk the lifted data to chromosome files
    _unlifted_chromosome_chunks_folder = "/usr/data/results/unlifted_chromosome_chunks/"
    _unlifted_chromosome_chunks_name = "unlifted_chromosome_chunks.bed"
    # chunk_the_chromosomes(_unlifted_input_data_filename, _unlifted_chromosome_chunks_folder, _unlifted_chromosome_chunks_name)


    _lifted_chromosome_chunks_folder = "/usr/data/results/lifted_chromosome_chunks/"
    _lifted_chromosome_chunks_name = "lifted_chromosome_chunks.bed"
    # chunk_the_chromosomes(_lifted_input_data_filename, _lifted_chromosome_chunks_folder, _lifted_chromosome_chunks_name)


    # Annotate the chromosome chunks
    # _annotated_unlifted_chromosome_chunks_folder = "/usr/data/results/annotated_unlifted_chromosome_chunks/"
    # _annotated_unlifted_chromosome_chunks_name = "annotated_unlifted_chromosome_chunks.bed"
    # # annotate_the_chromosomes(_unlifted_chromosome_chunks_folder, _annotated_unlifted_chromosome_chunks_folder, _annotated_unlifted_chromosome_chunks_name, _metadome_filename, parallel=True)
    _annotated_lifted_chromosome_chunks_folder = "/usr/data/results/annotated_lifted_chromosome_chunks/"
    _annotated_lifted_chromosome_chunks_name = "annotated_lifted_chromosome_chunks.bed"
    # annotate_the_chromosomes(_lifted_chromosome_chunks_folder, _annotated_lifted_chromosome_chunks_folder, _annotated_lifted_chromosome_chunks_name, _metadome_filename, parallel=True)

    # # Combine the annotated chromosome chunks back to a single file
    # _lifted_annotated_input_data_filename = "/usr/data/annotated_liftover_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz"
    # combine_annotated_chunks(_annotated_lifted_chromosome_chunks_folder, _lifted_annotated_input_data_filename)

    _lifted_annotated_input_data_filename = "/usr/data/sorted_annotated_liftover_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz"


    _unlifted_output_data_filename = "/usr/data/results/final_liftover_annotated_bed_genotype_freeze.6a.pass_only.phased.mesa_1331samples_combined_vep_consequences.txt.gz"

    combine_liftover(_unlifted_input_data_filename, _lifted_annotated_input_data_filename, _unlifted_output_data_filename)

    time_step = time.perf_counter()
    logging.getLogger(LOGGER_NAME).info("Finished annotation in "+str(time_step-start_time)+" seconds")