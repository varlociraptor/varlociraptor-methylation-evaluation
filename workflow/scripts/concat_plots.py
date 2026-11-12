import sys

import fitz  # pymupdf

sys.stderr = open(snakemake.log[0], "w")

input_pdfs = snakemake.input.plots
output_pdf = snakemake.output[0]
# Titles written above each tile (empty list = no titles)
titles = snakemake.params.get("titles", [])
# "horizontal": tiles side by side, "vertical": tiles stacked
direction = snakemake.params.get("direction", "horizontal")

if titles and len(titles) != len(input_pdfs):
    raise ValueError(f"Got {len(titles)} titles for {len(input_pdfs)} plots")
if direction not in ["horizontal", "vertical"]:
    raise ValueError(f"Unknown direction: {direction}")

title_height = 40 if titles else 0
font_size = 20

# Size of each input page
sizes = []
for pdf in input_pdfs:
    with fitz.open(pdf) as src:
        sizes.append((src[0].rect.width, src[0].rect.height))

if direction == "horizontal":
    page_width = sum(w for w, _ in sizes)
    page_height = max(h for _, h in sizes) + title_height
else:
    page_width = max(w for w, _ in sizes)
    page_height = sum(h + title_height for _, h in sizes)

out = fitz.open()
page = out.new_page(width=page_width, height=page_height)

x0, y0 = 0, 0
for i, (pdf, (w, h)) in enumerate(zip(input_pdfs, sizes)):
    if titles:
        title = titles[i]
        text_width = fitz.get_text_length(title, fontname="hebo", fontsize=font_size)
        # Centered above the tile, baseline slightly above the tile
        page.insert_text(
            (x0 + (w - text_width) / 2, y0 + title_height - 10),
            title,
            fontname="hebo",
            fontsize=font_size,
        )
    with fitz.open(pdf) as src:
        page.show_pdf_page(
            fitz.Rect(x0, y0 + title_height, x0 + w, y0 + title_height + h), src, 0
        )
    if direction == "horizontal":
        x0 += w
    else:
        y0 += h + title_height

out.save(output_pdf)
out.close()
