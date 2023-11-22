import pdfminer
import pdfminer.high_level
from pdfminer.high_level import extract_pages
for page_layout in extract_pages("totA2.pdf"):
    for element in page_layout:
        print(element)