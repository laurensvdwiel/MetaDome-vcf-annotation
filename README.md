# MetaDome-vcf-annotation

An implementation to vcf annotations using MetaDome data.

MetaDome web server is available [here](https://stuart.radboudumc.nl/metadome/)
MetaDome data is available [here](https://zenodo.org/record/6625251)

## Citation ##

If you make use of this software for academic purposes, please cite the article [10.1016/j.ajhg.2017.08.004](https://doi.org/10.1016/j.ajhg.2017.08.004).

## How to run ##

### Dependencies ###

Please ensure you have the following software installed on your machine:

	docker

You can get docker [here](https://www.docker.com/get-docker)

## Build Docker image ##

Run the following command in the root directory of this project

    docker build -t metadome_annotation .

### Example run ###

[//]: # (Example for a hypothetical gene)

[//]: # ()
[//]: # (    docker run --rm -v $&#40;pwd&#41;:/app --name metadome_annotation_container metadome_annotation python /app/src/spatial_clustering.py --gene_name=TESTGENE --variant_cDNA_locations=1,2,3,4,5,5,5,10,1000 --cDNA_length=1337 --n_permutations=10 --parallel=True --random_seed=1 --correction=1)

[//]: # ()
[//]: # (This should result in the following output:)

[//]: # (    )
[//]: # (    Computing spatial clustering for gene: TESTGENE)

[//]: # (    cDNA_length: 1337)

[//]: # (    variant_cDNA_locations: ['1', '2', '3', '4', '5', '5', '5', '10', '1000'])

[//]: # (    Settings: random_seed: 1, parallel: True, n_permutations: 10)

[//]: # ()
[//]: # (    Results:)

[//]: # (    gene: TESTGENE, with n variants: 9)

[//]: # (    geometric_mean: 6.898887607777424955370112958229102e-08)

[//]: # (    corrected p-value: 0.09090909090909091 &#40;Bonferroni correction = 1&#41;)

[//]: # ()
[//]: # (The arguments '--n_permutations', '--parallel', `--random_seed`, and, `--correction` are optional and need not be included.)

[//]: # ()
[//]: # (For further information on arguments, please refer to the help file:)

[//]: # ()
[//]: # (    -h, --help            show this help message and exit)

[//]: # (    --gene_name GENE_NAME)

[//]: # (                          &#40;Required&#41; Name of the gene of interest, example)

[//]: # (                          usage: --gene_name=BRCA1)

[//]: # (    --variant_cDNA_locations VARIANT_CDNA_LOCATIONS)

[//]: # (                          &#40;Required&#41; cDNA based variant locations, example)

[//]: # (                          usage: --variant_cDNA_locations=10,50,50,123)

[//]: # (    --cDNA_length CDNA_LENGTH)

[//]: # (                          &#40;Required&#41; total cDNA length of the gene &#40;including)

[//]: # (                          stop codon&#41;, example usage: --cDNA_length=1337)

[//]: # (    --n_permutations N_PERMUTATIONS)

[//]: # (                          &#40;Optional&#41; total nunber of permutations,)

[//]: # (                          default=100000000 &#40;1.00E+08&#41;, example usage:)

[//]: # (                          --n_permutations=100)

[//]: # (    --parallel PARALLEL   &#40;Optional&#41; should the algorithm make use of parallel)

[//]: # (                          computation?, default=True, example usage:)

[//]: # (                          --parallel=True)

[//]: # (    --random_seed RANDOM_SEED)

[//]: # (                          &#40;Optional&#41; The seed used for initialization of the)

[//]: # (                          random permutations, default=1, example usage:)

[//]: # (                          --random_seed=1)

[//]: # (    --correction CORRECTION)

[//]: # (                          &#40;Optional&#41; The number of genes the p-value must be)

[//]: # (                          corrected for in a Bonferonni manner, default=1,)

[//]: # (                          example usage: --correction=1)
