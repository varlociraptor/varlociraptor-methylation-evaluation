import sys
import pysam

# Redirect stderr to the Snakemake log file for proper workflow logging
sys.stderr = open(snakemake.log[0], "w")


def rename_chromosomes_bam(input_bam: str, output_bam: str) -> None:
    """
    Reads a BAM file and removes the 'chr' prefix from all chromosome names
    in both the header and alignment records.

    Parameters:
        input_bam (str): Path to the input BAM file.
        output_bam (str): Path where the modified BAM file will be written.
    """
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        # Extract and copy original BAM header
        original_header = infile.header.to_dict()

        # Initialize a new header dictionary
        new_header = {
            "HD": original_header.get("HD", {}),
            "RG": original_header.get("RG", []),
            "SQ": [],
        }

        # Create mapping between old and new reference names
        name_map = {}
        for sq in original_header.get("SQ", []):
            old_name = sq["SN"]
            new_name = (
                old_name.removeprefix("chr") if old_name.startswith("chr") else old_name
            )
            new_header["SQ"].append({**sq, "SN": new_name})
            name_map[old_name] = new_name

        # Map reference names to new IDs
        ref_name_to_id = {sq["SN"]: i for i, sq in enumerate(new_header["SQ"])}

        # Build a new header object
        new_header_obj = pysam.AlignmentHeader.from_dict(new_header)

        # Create the output BAM file with the updated header
        with pysam.AlignmentFile(output_bam, "wb", header=new_header_obj) as outfile:
            for read in infile:
                read_dict = read.to_dict()

                # Update reference names and IDs
                if read_dict["ref_name"] in name_map:
                    new_ref = name_map[read_dict["ref_name"]]
                    read_dict["ref_name"] = new_ref
                    read_dict["ref_id"] = ref_name_to_id[new_ref]

                if read_dict["next_ref_name"] in name_map:
                    new_next_ref = name_map[read_dict["next_ref_name"]]
                    read_dict["next_ref_name"] = new_next_ref
                    read_dict["next_ref_id"] = ref_name_to_id[new_next_ref]

                # Reconstruct and write the modified read
                new_read = pysam.AlignedSegment.from_dict(
                    read_dict, header=new_header_obj
                )
                outfile.write(new_read)


if __name__ == "__main__":
    # Get input/output paths from Snakemake
    input_bam = snakemake.input[0]
    output_bam = snakemake.output[0]

    # Execute main renaming function
    rename_chromosomes_bam(input_bam, output_bam)
