# Physikalischen Fortgeschrittenenpraktikum 2
Notes for the [advanced physics lab course 2 at the University of Konstanz](https://fp.physik.uni-konstanz.de).

## Installing the necessary requirements

Required software stack to interact with the files:
- [Git](https://git-scm.com/)
  also `git-lfs`
- [LaTeX](https://www.latex-project.org/)\
  With diverse extensions. On debian you can install them with `sudo apt install texlive texlive-fonts-extra`.
- [Python](https://www.python.org/)\
`sudo apt install python python-pip`
- Roboto font is used for plots. You can install it with `sudo apt install fonts-roboto`.
- [Julia](https://julialang.org/)\
  For the Julia scripts.

## Activating environments

To activate the **julia** environment run in the julia shell:
```julia
]
activate .
```

To run the **python** scripts create a virtual environment and install the requirements:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
Then to install the local `src` package for shared code:
```bash
pip install -e .
```