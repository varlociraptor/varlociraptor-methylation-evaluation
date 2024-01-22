import requests
import logging

url = snakemake.params.url
output_path = snakemake.output.alignment

logging.basicConfig(level=logging.INFO)

try:
    # Download the file
    response = requests.get(url, stream=True)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Write the file in chunks to handle large downloads
        with open(output_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    file.write(chunk)
        print("The file was successfully downloaded.")
    else:
        print(f"Error downloading the file. Status code: {
              response.status_code}")

except requests.exceptions.RequestException as e:
    print(f"An error occurred during the request: {e}")

except Exception as e:
    print(f"An unexpected error occurred: {e}")
