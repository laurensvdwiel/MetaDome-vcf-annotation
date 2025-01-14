import logging
import time
import gzip
import os
import re
import pandas as pd
from joblib.parallel import Parallel, delayed

from dev_settings import LOGGER_NAME
from CustomLogger import initLogging

def process_vcf_header_line(line):
    """Retrieves the headers of a VCF file
    """
    # decode if line is gzipped
    if type(line) is bytes:
        line = line.decode('utf-8')

    # split line into tokens
    header = line.strip().split('\t')

    return header

def process_vcf_line(line, header):
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

    # check if the number of tokens is as expected
    if len(processed_line) != len(header):
        logging.getLogger(LOGGER_NAME).warning("Not enough tokens found than expected for line '"+str(line)+"'")
        return None
    else:
        # craft new data entry
        data_entry = dict(zip(header, processed_line))

    # return data_entry
    return data_entry

def read_and_tokenize_vcf_file_line(file_reader, line_counter, report_name, header=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']):
    """Decorator function that reads and tokenizes a single vcf file line with logging service support
    """
    # Read a line of the lifted bed file
    line = file_reader.readline()
    # Check if the end of the file has been reached
    end_of_file = not line

    if not end_of_file:
        # Increment the unlifted line counter
        line_counter += 1

        #process the header line
        if line.startswith("#CHROM"):
            header = process_vcf_header_line(line)

        # skip the info header lines
        if line.startswith("#"):
            return line, None, header, line_counter, end_of_file

        # Tokenize the line
        processed_line = process_vcf_line(line, header)
    else:
        processed_line = None
        logging.getLogger(LOGGER_NAME).info(
            "Reached end of "+str(report_name)+" file at line " + str(line_counter))

    return line, processed_line, header, line_counter, end_of_file

def annotate_metadome(metadome_df_grouped, chromosome, chr_pos):
    """ Annotate the metadome data to the processed line
    """
    # Initialize the return value
    return_value = {}

    # Make sure the metadome_df provided contains data from only one chromosome
    chromosomes_in_metadome_df = metadome_df_grouped.groups.keys()

    if not chromosome in chromosomes_in_metadome_df:
        logging.getLogger(LOGGER_NAME).debug("The metadome_df provided does not contain data from the chromosome '"+chromosome+"'")
        return return_value

    # check for a metadome hit
    metadome_df = metadome_df_grouped.get_group(chromosome)
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

def annotate_single_vcf(vcf_filename, metadome_df_grouped, target_folder, report_counter):
    """
    Annotate a single VCF file with metadome data
    """

    logging.getLogger(LOGGER_NAME).info("Starting annotation of MetaDome data to '" + vcf_filename + "'")

    # Initialize counter
    processed_line_counter = 0

    # create source and target file reader/writer contexts based on the file extension
    if vcf_filename.endswith(".vcf.gz"):
        source_f = gzip.open(vcf_filename, 'rt')
        target_f_name = os.path.join(target_folder, "MetaDome_annotated_" + os.path.basename(vcf_filename).replace(".vcf.gz", ".vcf"))
    elif vcf_filename.endswith(".vcf"):
        source_f = open(vcf_filename, 'r')
        target_f_name = os.path.join(target_folder, "MetaDome_annotated_" + os.path.basename(vcf_filename))
    else:
        logging.getLogger(LOGGER_NAME).error("File '" + vcf_filename + "' is not a VCF file")
        return

    # create the target file
    target_f = open(target_f_name, 'a')

    # read the first line
    line, processed_line, header, processed_line_counter, end_of_file = read_and_tokenize_vcf_file_line(source_f, processed_line_counter, "Annotating " + os.path.basename(vcf_filename))

    while not end_of_file:
        if processed_line_counter % report_counter == 0:
            logging.getLogger(LOGGER_NAME).info("Processed '" + str(processed_line_counter) + "' lines for file '" + os.path.basename(vcf_filename) + "'")

        # First check if the line is a header line
        if line.startswith("#"):
            # check if it is the generic column info header line
            if line.startswith("#CHROM"):
                # add the metadome annotation to the header
                if type(line) is bytes:
                    line_to_write = line.decode('utf-8').strip()
                else:
                    line_to_write = line.strip()
                target_f.write("##INFO=<ID=MetaDome,Number=.,Type=String,Description=\"MetaDome domain annotation\">\n")

                line_to_write += "\t" + 'metadome_metadomain_positions' + "\n"

                # write the line to the target file
                target_f.write(line_to_write)
            else:
                target_f.write(line)
            line, processed_line, header, processed_line_counter, end_of_file = read_and_tokenize_vcf_file_line(
                source_f, processed_line_counter, "Annotating " + os.path.basename(vcf_filename), header)
            continue

        # Annotate the line with the metadome data
        processed_line['metadome_metadomain_positions'] = "-"

        if metadome_df_grouped is not None:
            annotation = annotate_metadome(metadome_df_grouped, processed_line['#CHROM'], processed_line['POS'])
            if len(annotation) > 0:
                processed_line['metadome_metadomain_positions'] = annotation['metadome_metadomain_positions']

        # write the line to the target file
        if type(line) is bytes:
            line_to_write = line.decode('utf-8').strip()
        else:
            line_to_write = line.strip()
        line_to_write += "\t" + processed_line['metadome_metadomain_positions'] + "\n"

        # write the line to the target file
        target_f.write(line_to_write)

        # read the next line
        line, processed_line, header, processed_line_counter, end_of_file = read_and_tokenize_vcf_file_line(source_f, processed_line_counter, "Annotating " + os.path.basename(vcf_filename), header)

    target_f.close()
    logging.getLogger(LOGGER_NAME).info("Finished annotation of MetaDome data for '" + vcf_filename + "' to '" + target_f_name + "'")


def annotate_metadome_to_vcf_files(source_folder, target_folder, metadome_filename, n_jobs=None, redo_previous_files = False, parallel = True, report_counter=10000):
    """
    Annotate MetaDome data to the VCF files in the source folder and store the annotated files in the target folder
    Schedules the sequential or parallel processing of the files
    """
    # read the source file
    logging.getLogger(LOGGER_NAME).info("Starting annotation of the chunked chromosome files from '"+source_folder+"'")

    # Create the target folder if it doesn't exist
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # Check the files in the source folder
    source_files = []
    for file in os.listdir(source_folder):
        if file.endswith(".vcf") or file.endswith(".vcf.gz"):
            # add the file to the list of files to annotate
            source_files.append(os.path.join(source_folder, file))

            logging.getLogger(LOGGER_NAME).info("Added '"+source_folder+file+"' to the list of files to annotate")
        else:
            logging.getLogger(LOGGER_NAME).warning("Skipping file '"+file+"' in source folder '"+source_folder+"'")

    # report on the number of files to annotate
    logging.getLogger(LOGGER_NAME).info("Found '"+str(len(source_files))+"' files in the source folder '"+source_folder+"'")
    logging.getLogger(LOGGER_NAME).info("Checking target folder '"+target_folder+"' for already annotated files...")

    # Check if the files have already been annotated
    target_files = [x for x in os.listdir(target_folder) if x.startswith("MetaDome_annotated_")]
    # Count if any of the target files overlap with annotated versions of source file names
    target_files_overlap = [x for x in target_files if any([re.match("MetaDome_annotated_"+os.path.basename(y), x) for y in source_files])]
    logging.getLogger(LOGGER_NAME).info("Found '"+str(len(target_files))+"' files in the target folder '"+target_folder+"' of which '"+str(len(target_files_overlap))+"' overlap with the source files")
    if redo_previous_files:
        # remove files that have already been annotated
        for file in target_files_overlap:
            os.remove(file)
    else:
        # remove the files that have already been annotated from the source files
        source_files = [x for x in source_files if not any([re.match("MetaDome_annotated_"+os.path.basename(x), y) for y in target_files])]
        logging.getLogger(LOGGER_NAME).info("Removed '"+str(len(source_files))+"' files from the source files that have already been annotated")
    
    # Read the metadome annotation
    metadome_df = pd.read_csv(metadome_filename, index_col=False)

    # Select to only positions that have domain annotation
    metadome_df = metadome_df[(metadome_df.domain_id == metadome_df.domain_id)]

    # Drop duplicates, if any
    metadome_df = metadome_df[['chrom', 'pos_start', 'pos_stop', 'strand', 'symbol', 'domain_id', 'consensus_pos']].drop_duplicates()

    # Group the metadome data by chromosome
    metadome_df_grouped = metadome_df.groupby('chrom')

    if parallel:
        if n_jobs is None:
            n_jobs = len(source_files)
        # Report on the number of jobs to be run
        logging.getLogger(LOGGER_NAME).info("Starting annotation of the VCF files in parallel with '"+str(n_jobs)+"' jobs")
        _ = Parallel(n_jobs=n_jobs)(delayed(annotate_single_vcf)(source_file, metadome_df_grouped, target_folder, report_counter) for source_file in source_files)
    else:
        # Report on the number of jobs to be run
        logging.getLogger(LOGGER_NAME).info("Starting annotation of the '"+str(len(source_files))+"' VCF files sequentially")
        for source_file in source_files:
            annotate_single_vcf(source_file, metadome_df_grouped, target_folder, report_counter)

def main(source_vcf_folder, target_vcf_folder, metadome_filename, parallel, n_jobs, redo_previous_files):
    logging.getLogger(LOGGER_NAME).info("Starting annotation")
    start_time = time.perf_counter()

    # Annotate the multiple VCF files with MetaDome data
    annotate_metadome_to_vcf_files(source_vcf_folder, target_vcf_folder, metadome_filename, parallel, n_jobs, redo_previous_files)

    time_step = time.perf_counter()
    logging.getLogger(LOGGER_NAME).info("Finished annotation in " + str(time_step - start_time) + " seconds")

if __name__ == '__main__':
    import argparse

    # Initialize the ArgumentParser
    parser = argparse.ArgumentParser(description=' Running MetaDome annotation on VCF \n Expects a folder with one or more VCF files that contain chr in the filename and the MetaDome annotation file as input')
    # required arguments
    parser.add_argument('--source_vcf_folder', type=str, required=True, help='(Required) Folder with VCF files to annotate')
    parser.add_argument('--target_vcf_folder', type=str, required=True, help='(Required) Folder to store the annotated VCF files')
    parser.add_argument('--metadome_filename', type=str, required=True, help='(Required) MetaDome annotation file')
    parser.add_argument('--parallel', type=bool, required=False, default=True, help='(Optional) should the algorithm make use of parallel computation?, default=True')
    parser.add_argument('--n_jobs', type=int, required=False, default=None, help='(Optional) number of jobs to run in parallel, default=None')
    parser.add_argument('--redo_previous_files', type=bool, required=False, default=False, help='(Optional) should the algorithm redo the annotation of files that have already been annotated?, default=False')
    parser.add_argument('--logging_to_console', type=bool, required=False, default=False, help='(Optional) should the algorithm log to the console?, default=False')

    # parse the arguments
    args = parser.parse_args()

    if args.logging_to_console:
        initLogging(logging_level=logging.DEBUG, print_to_console=True)
    else:
        # Initialize the logger
        initLogging(logging_level=logging.INFO)

    # check if the number of jobs to run in parallel is set but parallel processing is not set
    if args.n_jobs is not None and not args.parallel:
        logging.getLogger(LOGGER_NAME).warning("The number of jobs to run in parallel is set to '" + str(args.n_jobs) + "' but parallel processing is set to False")

    # run the main function
    main(source_vcf_folder=args.source_vcf_folder, target_vcf_folder=args.target_vcf_folder, metadome_filename=args.metadome_filename, parallel=args.parallel, n_jobs=args.n_jobs, redo_previous_files=args.redo_previous_files)