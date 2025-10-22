import sys
import requests
import logging

# Redirect stderr to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")

# Parameters from Snakemake
url = snakemake.params.url
output_file = snakemake.output.alignment

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def download_file(url: str, destination: str) -> None:
    """
    Download a file from a given URL and save it locally.

    Parameters:
        url (str): Remote file URL.
        destination (str): Local path where the file will be saved.
    """
    try:
        logging.info(f"Starting download: {url}")
        with requests.get(url, stream=True) as response:
            response.raise_for_status()  # Raise exception for HTTP errors
            with open(destination, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
        logging.info(f"File successfully downloaded to: {destination}")

    except requests.exceptions.RequestException as e:
        logging.error(f"Download failed: {e}")
        raise

    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        raise


if __name__ == "__main__":
    download_file(url, output_file)
