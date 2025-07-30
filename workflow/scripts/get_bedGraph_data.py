import requests
import os
import sys
import time
import shutil
import tempfile

# Redirect stderr to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")

# Parameters from Snakemake
bedGraphs = snakemake.params["bedGraphs"]
output_dir = snakemake.params["bedgraph_path"]
max_download_attempts = 3
backoff_factor = 2  # seconds

# Create a temporary directory for downloads (because of cluster)
temp_dir = tempfile.mkdtemp()

for bedGraph in bedGraphs:
    accession_number = bedGraph.split("_")[0]
    source_url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5649nnn/"
        + accession_number
        + "/suppl/"
        + bedGraph
        + ".bedGraph.gz"
    )
    final_output_path = os.path.join(output_dir, bedGraph + ".bedGraph.gz")
    temp_output_path = os.path.join(temp_dir, bedGraph + ".bedGraph.gz")

    os.makedirs(os.path.dirname(final_output_path), exist_ok=True)

    download_successful = False

    for attempt in range(1, max_download_attempts + 1):
        try:
            # Attempt to download with streaming
            with requests.get(source_url, stream=True, timeout=30) as response:
                response.raise_for_status()  # Raise error for bad responses
                with open(temp_output_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            # Move file from temp dir to final location
            shutil.move(temp_output_path, final_output_path)
            download_successful = True
            break  # Exit retry loop if successful

        except Exception as e:
            print(
                f"Error downloading {source_url} (attempt {attempt} of {max_download_attempts}): {e}"
            )
            # Clean up partial file if it exists
            if os.path.exists(temp_output_path):
                os.remove(temp_output_path)
            # Wait before retrying (exponential backoff)
            sleep_time = backoff_factor * (2 ** (attempt - 1))
            print(f"Waiting {sleep_time} seconds before retrying...")
            time.sleep(sleep_time)

    if not download_successful:
        raise Exception(
            f"Failed to download {source_url} after {max_download_attempts} attempts."
        )

# Clean up temporary download directory
shutil.rmtree(temp_dir)
