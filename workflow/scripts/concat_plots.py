import fitz  # pymupdf

input_pdfs = snakemake.input.plots
output_pdf = snakemake.output[0]

titles = [
    "PacBio and MethylSeq",
    "Nanopore and PacBio",
    "Nanopore and MethylSeq",
]

# Größe der ersten PDF als Referenz
first_doc = fitz.open(input_pdfs[0])
first_page = first_doc[0]
w = first_page.rect.width
h = first_page.rect.height
first_doc.close()

title_height = 40

# Neue PDF mit einer Seite
out = fitz.open()
page = out.new_page(width=3 * w, height=h + title_height)

for i, pdf in enumerate(input_pdfs):
    x0 = i * w
    src = fitz.open(pdf)

    page.show_pdf_page(
        fitz.Rect(x0, title_height, x0 + w, title_height + h),
        src,
        0,
    )

    src.close()

out.save(output_pdf)
out.close()
