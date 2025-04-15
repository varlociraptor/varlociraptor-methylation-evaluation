import urllib.request
import os

# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


bedGraphs = snakemake.params["bedGraphs"]
output_dir = snakemake.params["bedgraph_path"]


max_download_attempts = 3
for bedGraph in bedGraphs:
    accession_number = bedGraph.split("_")[0]

    source_url = (
        "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5649nnn/"
        + accession_number
        + "/suppl/"
        + bedGraph
        + ".bedGraph.gz"
    )
    output_file = os.path.join(output_dir, bedGraph + ".bedGraph.gz")

    download_successful = False
    for attempt in range(1, max_download_attempts + 1):
        try:
            urllib.request.urlretrieve(source_url, output_file)
            download_successful = True
            break  # Break out of the loop if download is successful
        except Exception as e:
            print(
                f"Error downloading {source_url} (attempt {
            attempt} / {max_download_attempts}): {e}"
            )
            os.remove(output_file)

    if not download_successful:
        print(
            f"Failed to download {source_url} after {
        max_download_attempts} attempts.Exiting."
        )
        raise Exception(
            f"Failed to download {source_url} after {
        max_download_attempts} attempts."
        )
