# Physikalischen Fortgeschrittenenpraktikum 2
Notes for the [advanced physics lab course 2 at the University of Konstanz](https://fp.physik.uni-konstanz.de).

Required software stack to interact with the files:
- [Git](https://git-scm.com/)
- [LaTeX](https://www.latex-project.org/)\
  With diverse extensions. On debian you can install them with `sudo apt install texlive texlive-fonts-extra`.
- [Python](https://www.python.org/)\
`sudo apt install python python-pip`

To run the python scripts create a virtual environment and install the requirements:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
Then to install the local `src` folder:
```bash
pip install -e src
```