import requests

url = snakemake.params.url
output_path = snakemake.output
# Download the file
response = requests.get(url)

# Check if the request was successful (status code 200)
if response.status_code == 200:
    # Write the file in binary mode
    with open(output_path, 'wb') as file:
        file.write(response.content)
    print("The file was successfully downloaded.")
else:
    print(f"Error downloading the file. Status code: {response.status_code}")

