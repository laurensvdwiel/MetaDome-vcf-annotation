import logging
import time
import gzip
import os
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
        return_value['metadome_symbol'] = "MetaDome_Symbol="+",".join(set([x['symbol'] for x in records]))
        return_value['metadome_metadomain_positions'] = "MetaDome_POS="+",".join(set([str(x['domain_id']) + ":" + str(int(x['consensus_pos'])) for x in records]))

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
                # Add the metadome annotation to the header
                target_f.write("##INFO=<ID=MetaDome_POS,Number=.,Type=String,Description=\"MetaDome consensus position in PFAM domain annotation\">\n")
                target_f.write("##INFO=<ID=MetaDome_Symbol,String=.,Type=String,Description=\"MetaDome Havana gene Symbol where PFAM is present\">\n")

                # Write the header line
                if type(line) is bytes:
                    line_to_write = line.decode('utf-8').strip()
                else:
                    line_to_write = line.strip()
                # write the line to the target file
                target_f.write(line_to_write)
            else:
                target_f.write(line)
            line, processed_line, header, processed_line_counter, end_of_file = read_and_tokenize_vcf_file_line(
                source_f, processed_line_counter, "Annotating " + os.path.basename(vcf_filename), header)
            continue

        # Annotate the line with the metadome data
        metadome_annotation = {}
        if metadome_df_grouped is not None:
            metadome_annotation = annotate_metadome(metadome_df_grouped, processed_line['#CHROM'], processed_line['POS'])

        # Add the metadome annotation to the INFO field
        if len(metadome_annotation) > 0:
            processed_line['INFO'] = processed_line['INFO'] + ";" + metadome_annotation['metadome_symbol'] + ";" + metadome_annotation['metadome_metadomain_positions']
            # Update the line to write to the updated processed_line
            line_to_write = "\t".join(str(x) for x in processed_line.values()) + "\n"
        elif type(line) is bytes:
            line_to_write = line.decode('utf-8').strip()
        else:
            line_to_write = line.strip()

        # write the line to the target file
        target_f.write(line_to_write)

        # read the next line
        line, processed_line, header, processed_line_counter, end_of_file = read_and_tokenize_vcf_file_line(source_f, processed_line_counter, "Annotating " + os.path.basename(vcf_filename), header)

    target_f.close()
    logging.getLogger(LOGGER_NAME).info("Finished annotation of MetaDome data for '" + vcf_filename + "' to '" + target_f_name + "'")

def group_and_load_metadome_data(metadome_filename):
    """
    Group the metadome data by chromosome and load the data
    """
    # Check if the metadome file exists
    if not os.path.exists(metadome_filename):
        # exit and report the error
        logging.getLogger(LOGGER_NAME).error("The metadome file '" + metadome_filename + "' does not exist")
        raise FileNotFoundError("The metadome file '" + metadome_filename + "' does not exist")

    # Read the metadome annotation
    metadome_df = pd.read_csv(metadome_filename, index_col=False)

    # Select to only positions that have domain annotation
    metadome_df = metadome_df[(metadome_df.domain_id == metadome_df.domain_id)]

    # Drop duplicates, if any
    metadome_df = metadome_df[['chrom', 'pos_start', 'pos_stop', 'strand', 'symbol', 'domain_id', 'consensus_pos']].drop_duplicates()

    # Group the metadome data by chromosome
    metadome_df_grouped = metadome_df.groupby('chrom')

    return metadome_df_grouped

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
    # Count if any of the target files overlap with annotated versions of source file names, without the gz extension
    target_files_overlap = [x for x in source_files if "MetaDome_annotated_" + os.path.basename(x).replace(".vcf.gz", ".vcf") in target_files]
    logging.getLogger(LOGGER_NAME).info("Found '"+str(len(target_files))+"' files in the target folder '"+target_folder+"' of which '"+str(len(target_files_overlap))+"' overlap with the source files")
    if redo_previous_files:
        # remove files that have already been annotated
        for file in target_files:
            os.remove(os.path.join(target_folder, file))
    else:
        # remove the files that have already been annotated from the source files
        old_source_files = source_files
        source_files = [x for x in source_files if x not in target_files_overlap]
        logging.getLogger(LOGGER_NAME).info("Removed '"+str(len(old_source_files)-len(source_files))+"' files from the initial '"+str(len(old_source_files))+"' source files, due to previous annotation")

    # Check if there are any files to annotate
    if len(source_files) == 0:
        logging.getLogger(LOGGER_NAME).info("No files to annotate")
        return
    elif len(source_files) == 1:
        parallel = False
        logging.getLogger(LOGGER_NAME).info("Only one file to annotate, running in sequential mode")

    # Read the metadome annotation and group by chromosome
    metadome_df_grouped = group_and_load_metadome_data(metadome_filename)

    if parallel:
        # Check if the number of jobs to run in parallel is set and if there if it is greater or equal to the number of files to be annotated
        if n_jobs is None or n_jobs > len(source_files):
            logging.getLogger(LOGGER_NAME).info("Adjusted the user provided number of jobs to run in parallel from '" + str(n_jobs) + "' to the number of '" + str(len(source_files)) + "' files to annotate")
            n_jobs = len(source_files)
        # Report on the number of jobs to be run
        logging.getLogger(LOGGER_NAME).info("Starting annotation of the VCF files in parallel with '"+str(n_jobs)+"' jobs")
        _ = Parallel(n_jobs=n_jobs)(delayed(annotate_single_vcf)(source_file, metadome_df_grouped, target_folder, report_counter) for source_file in source_files)
    else:
        # Report on the number of jobs to be run
        logging.getLogger(LOGGER_NAME).info("Starting annotation of the '"+str(len(source_files))+"' VCF files sequentially")
        for source_file in source_files:
            annotate_single_vcf(source_file, metadome_df_grouped, target_folder, report_counter)

def main_single(source_vcf_file, target_folder, metadome_filename, redo_previous_file, report_counter=10000):
    logging.getLogger(LOGGER_NAME).info("Starting annotation")
    try:
        start_time = time.perf_counter()

        # Check if the source file exists
        if not os.path.exists(source_vcf_file):
            logging.getLogger(LOGGER_NAME).error("The source file '"+source_vcf_file+"' does not exist")
            return

        # Check if the target folder exists
        if not os.path.exists(target_folder):
            # exit and report the error
            logging.getLogger(LOGGER_NAME).error("The target folder '"+target_folder+"' does not exist")
            return

        # Check if the target file already exists
        # bzgipped target file
        bgziped_target_file = os.path.join(target_folder, "MetaDome_annotated_" + os.path.basename(source_vcf_file))
        # non bgzipped target file
        target_file = os.path.join(target_folder, "MetaDome_annotated_" + os.path.basename(source_vcf_file).replace(".vcf.gz", ".vcf"))

        # check if bgzipped file exists
        if os.path.exists(bgziped_target_file) or os.path.exists(target_file):
            if redo_previous_file:
                if os.path.exists(bgziped_target_file):
                    logging.getLogger(LOGGER_NAME).info("Removing previous bgzipped target file '"+bgziped_target_file+"'")
                    os.remove(bgziped_target_file)
                if os.path.exists(target_file):
                    logging.getLogger(LOGGER_NAME).info("Removing previous target file '"+target_file+"'")
                    os.remove(target_file)
            else:
                if os.path.exists(bgziped_target_file):
                    logging.getLogger(LOGGER_NAME).info("The bgzipped target file '"+bgziped_target_file+"' already exists, canceling annotation")
                if os.path.exists(target_file):
                    logging.getLogger(LOGGER_NAME).info("The target file '"+target_file+"' already exists, canceling annotation")
                return

        # Read the metadome annotation and group by chromosome
        metadome_df_grouped = group_and_load_metadome_data(metadome_filename)

        # Annotate the single VCF file with MetaDome data
        annotate_single_vcf(source_vcf_file, metadome_df_grouped, target_folder, report_counter)

        time_step = time.perf_counter()
        logging.getLogger(LOGGER_NAME).info("Finished annotation in " + str(time_step - start_time) + " seconds")
    except Exception as e:
        logging.getLogger(LOGGER_NAME).error("An error occurred during the setup of the annotation process: " + str(e) + "\n" + str(e.with_traceback))
        raise e

def main_multi(source_vcf_folder, target_vcf_folder, metadome_filename, parallel, n_jobs, redo_previous_files):
    logging.getLogger(LOGGER_NAME).info("Starting annotation")
    try:
        start_time = time.perf_counter()

        # Annotate the multiple VCF files with MetaDome data
        annotate_metadome_to_vcf_files(source_folder=source_vcf_folder, target_folder=target_vcf_folder, metadome_filename=metadome_filename, parallel=parallel, n_jobs=n_jobs, redo_previous_files=redo_previous_files)

        time_step = time.perf_counter()
        logging.getLogger(LOGGER_NAME).info("Finished annotation in " + str(time_step - start_time) + " seconds")
    except Exception as e:
        logging.getLogger(LOGGER_NAME).error("An error occurred during the annotation process: " + str(e) + "\n" + str(e.with_traceback))
        raise e

if __name__ == '__main__':
    import argparse

    # Initialize the ArgumentParser
    _parser = argparse.ArgumentParser(description=' Running MetaDome annotation on VCF \n Expects a folder with one or more VCF files that contain chr in the filename and the MetaDome annotation file as input')

    # required arguments determining single or multi mode and logging
    _parser.add_argument("annotation_mode", nargs='?', type=str, help="The function to run, either 'single' or 'multi'",
                        choices=['single', 'multi'], default='single')
    _parser.add_argument('--logging_to_console', type=bool, required=False, default=False,
                        help='(Optional) should the algorithm log to the console?, default=False')

    # choose the single mode
    _args, _sub_args = _parser.parse_known_args()

    # Initialize the logger
    if _args.logging_to_console:
        initLogging(logging_level=logging.DEBUG, print_to_console=True)
    else:
        initLogging(logging_level=logging.INFO)

    # single mode
    if _args.annotation_mode == "single":
        _sub_parser = argparse.ArgumentParser()

        # required arguments
        _sub_parser.add_argument('--source_vcf_file', type=str, required=True, help='(Required) VCF file to annotate')
        _sub_parser.add_argument('--target_folder', type=str, required=True, help='(Required) Folder to store the annotated VCF files')
        _sub_parser.add_argument('--metadome_filename', type=str, required=True, help='(Required) MetaDome annotation file')
        _sub_parser.add_argument('--redo_previous_file', type=bool, required=False, default=False, help='(Optional) should the algorithm redo the annotation of files that have already been annotated?, default=False')

        # parse the arguments
        _mode_args = _sub_parser.parse_args(_sub_args)

        # run the main_single function
        main_single(source_vcf_file=_mode_args.source_vcf_file,
            target_folder=_mode_args.target_folder,
            metadome_filename=_mode_args.metadome_filename,
            redo_previous_file=_mode_args.redo_previous_file)

    # multi mode
    elif _args.annotation_mode == "multi":
        _sub_parser = argparse.ArgumentParser()

        # required arguments
        _sub_parser.add_argument('--source_vcf_folder', type=str, required=True,
                            help='(Required) Folder with VCF files to annotate')
        _sub_parser.add_argument('--target_vcf_folder', type=str, required=True,
                            help='(Required) Folder to store the annotated VCF files')
        _sub_parser.add_argument('--metadome_filename', type=str, required=True, help='(Required) MetaDome annotation file')
        _sub_parser.add_argument('--parallel', type=bool, required=False, default=True,
                            help='(Optional) should the algorithm make use of parallel computation?, default=True')
        _sub_parser.add_argument('--n_jobs', type=int, required=False, default=None,
                            help='(Optional) number of jobs to run in parallel, default=None')
        _sub_parser.add_argument('--redo_previous_files', type=bool, required=False, default=False,
                            help='(Optional) should the algorithm redo the annotation of files that have already been annotated?, default=False')

        # parse the arguments
        _mode_args = _sub_parser.parse_args(_sub_args)

        # check if the number of jobs to run in parallel is set but parallel processing is not set
        if _mode_args.n_jobs is not None and not _mode_args.parallel:
            logging.getLogger(LOGGER_NAME).warning("The number of jobs to run in parallel is set to '" + str(
                _mode_args.n_jobs) + "' but parallel processing is set to False")

        # run the main_multi function
        main_multi(source_vcf_folder=_mode_args.source_vcf_folder, target_vcf_folder=_mode_args.target_vcf_folder,
             metadome_filename=_mode_args.metadome_filename, parallel=_mode_args.parallel, n_jobs=_mode_args.n_jobs,
             redo_previous_files=_mode_args.redo_previous_files)
    else:
        logging.getLogger(LOGGER_NAME).error("The annotation mode provided is not valid: '"+str(_args.annotation_mode)+"'")