import pandas as pd
from joblib.parallel import Parallel, delayed

def annotate_chunk(metadome_df_chunk, data_chunk):
    variants_that_effect_domains = []
    for index, row in data_chunk.iterrows():
        # craft new data entry
        data_entry = row.to_dict()
        # check for a metadome hit
        metadome_hit = metadome_df_chunk[(metadome_df_chunk.pos_start <= int(data_entry['pos_lifted']))& (int(data_entry['pos_lifted']) <= metadome_df_chunk.pos_stop)]
        if len(metadome_hit)>0:
            # add record
            records = metadome_hit.to_dict('records')
            data_entry['metadome_symbol'] = ";".join(set([x['symbol'] for x in records]))
            data_entry['metadome_metadomain_positions'] = ";".join(set([str(x['domain_id'])+":"+str(x['consensus_pos']) for x in records]))
        else:
            # add empty data
            data_entry['metadome_symbol'] = ""
            data_entry['metadome_metadomain_positions'] = ""

        variants_that_effect_domains.append(data_entry)
    return variants_that_effect_domains


def annotate_dataset_with_metadome(metadome_filename, input_data_filename, output_data_filename, parallel=True):
    """Annotates a pandas DataFrame containing with """

    # Read the metadome annotation
    metadome_df = pd.read_csv(metadome_filename, index_col=False)

    # Select to only positions that have domain annotation
    metadome_df = metadome_df[(metadome_df.domain_id==metadome_df.domain_id)]

    # Drop duplicates, if any
    metadome_df = metadome_df[['chrom', 'pos_start', 'pos_stop', 'strand', 'symbol', 'domain_id', 'consensus_pos']].drop_duplicates()

    # Read the metadome annotation
    input_data_df = pd.read_csv(input_data_filename, index_col=False, sep='\t')

    # Divide the work into chromosome chunks
    metadome_chromosomes = pd.unique(metadome_df.chrom)
    output_annotations = []
    if parallel:
        output_annotations = Parallel(n_jobs=len(metadome_chromosomes))(delayed(annotate_chunk)(metadome_df[metadome_df.chrom == chromosome], input_data_df[input_data_df['chrom_lifted'] == chromosome]) for chromosome in metadome_chromosomes)
    else:
        output_annotations = [annotate_chunk(metadome_df[metadome_df.chrom == chromosome], input_data_df[input_data_df['chrom_lifted'] == chromosome]) for chromosome in metadome_chromosomes]

    # combine again
    result_data_df = []
    for output in output_annotations:
        result_data_df += output

    # convert to dataframe and save results
    result_data_df = pd.DataFrame(result_data_df)
    result_data_df.to_csv(output_data_filename, index=False, sep="\t")