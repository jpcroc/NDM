import pdfminer
from pdfminer.high_level import extract_pages
for page_layout in extract_pages("cat3A1.pdf"):
    for element in page_layout:
        print(element)