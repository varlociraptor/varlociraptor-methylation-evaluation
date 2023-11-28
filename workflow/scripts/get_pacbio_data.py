import requests

url = "https://downloads.pacbcloud.com/public/dataset/HG002-CpG-methylation-202202/HG002.GRCh38.haplotagged.bam"
target_file_path = "/home/adrian/Documents/Promotion/varlociraptor-methylation-evaluation/resources/pacbio/HG002.GRCh38.haplotagged.bam"

# Download the file
response = requests.get(url)

# Check if the request was successful (status code 200)
if response.status_code == 200:
    # Write the file in binary mode
    with open(target_file_path, 'wb') as file:
        file.write(response.content)
    print("The file was successfully downloaded.")
else:
    print(f"Error downloading the file. Status code: {response.status_code}")

