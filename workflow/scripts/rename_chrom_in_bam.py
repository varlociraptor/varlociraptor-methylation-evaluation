import pysam

# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


def rename_chromosomes_bam(input_bam, output_bam):
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        # Original-Header extrahieren und kopieren
        original_header = infile.header.to_dict()

        # Neuen Header bauen (ohne 'chr')
        new_header = {"HD": original_header.get("HD", {})}
        new_sq = []
        name_map = {}  # alt_name -> neuer Name

        for sq in original_header.get("SQ", []):
            old_name = sq["SN"]
            new_name = (
                old_name.removeprefix("chr") if old_name.startswith("chr") else old_name
            )
            new_sq.append({**sq, "SN": new_name})
            name_map[old_name] = new_name

        new_header["SQ"] = new_sq

        # Mapping von Referenznamen auf neue Referenz-IDs
        ref_name_to_id = {sq["SN"]: i for i, sq in enumerate(new_sq)}

        # Konvertiere neuen Header zurück in pysam-kompatibles Format
        new_header_obj = pysam.AlignmentHeader.from_dict(new_header)

        with pysam.AlignmentFile(output_bam, "wb", header=new_header_obj) as outfile:
            for read in infile:
                read_dict = read.to_dict()

                # ref_name + ref_id ändern
                if read_dict["ref_name"] in name_map:
                    new_ref_name = name_map[read_dict["ref_name"]]
                    read_dict["ref_name"] = new_ref_name
                    read_dict["ref_id"] = ref_name_to_id[new_ref_name]

                # next_ref_name + next_ref_id ändern
                if read_dict["next_ref_name"] in name_map:
                    new_next_ref_name = name_map[read_dict["next_ref_name"]]
                    read_dict["next_ref_name"] = new_next_ref_name
                    read_dict["next_ref_id"] = ref_name_to_id[new_next_ref_name]

                # Neues Read-Objekt schreiben
                new_read = pysam.AlignedSegment.from_dict(
                    read_dict, header=new_header_obj
                )
                outfile.write(new_read)


if __name__ == "__main__":
    input_bam = snakemake.input[0]
    output_bam = snakemake.output[0]
    rename_chromosomes_bam(input_bam, output_bam)
