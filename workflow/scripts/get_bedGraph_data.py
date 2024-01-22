import urllib.request

bedGraphs = snakemake.params["bedGraphs"]

for bedGraph in bedGraphs:

    accession_number = bedGraph.split('_')[0]
    urllib.request.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5649nnn/' + accession_number +
                               '/suppl/' + bedGraph + '.bedGraph.gz', 'resources/HG002/' + bedGraph + '.bedGraph.gz')
