## Create jupyter-book_howto

- Create `.ipynb` files as documentation **and** execution (... experimentation).
- Example on is246206, change SANDBOX repository *at yours needs*.

```bash
export MYROOT=/export/home/catA/wambeke/NDM_ROODIR
alias cd_doc="cd ${MYROOT}/NDM/docs/doc_2023"

f_init_conda
conda activate py3qt5
conda list | grep jupyter-book
which jupyter-notebook
which jupyter-book

cd_doc
jupyter-notebook
```

- **... and now you make .ipynb files (as NDM_howto.ipynb etc...) at yours needs**
- **... and after you copy yours .ipynb files etc... to NDM/docs/..., and modify docs/.../_toc.yml**



## Create documentation html plus pdf with jupyter-book

- On is246206

### Prerequisites jupyter-book

```bash
pip install jupyter-book

sudo dnf remove texlive
sudo dnf install texlive-scheme-full

# https://unix.stackexchange.com/questions/199638/how-to-fully-install-latex-in-fedora
# sudo dnf install latexmk
# sudo dnf install texlive-pgfplots
# sudo dnf install texlive-bbm

```

### Create doc html plus pdf

```bash
cd_doc

# to get first example _config.yml and _toc.yml
# jupyter-book create example_tmp
# cp example_tmp._* .
# rm -rf example_tmp

export BOOK=$(pwd)             # as first empty example

# jupyter-book create ${BOOK}
tree ${BOOK}
rm -rf ${BOOK}/_build ${BOOK}/_tmp ${BOOK}/__pycache__

# html
jupyter-book build ${BOOK}
firefox ${BOOK}/_build/html/index.html

# pdf OK
# using CFI.colorize_file_to_markdown("PACKAGESPY_howto.py") etc.
jupyter-book build ${BOOK} --builder pdflatex
evince ${BOOK}/_build/latex/PACKAGESPY.pdf
cp ${BOOK}/_build/latex/PACKAGESPY.pdf ${BOOK}/.

```
