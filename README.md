[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Conda](https://anaconda.org/russel88/cctyper/badges/installer/conda.svg)](https://anaconda.org/russel88/cctyper) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# PyCCTyper

Fork of [CCTyper](https://github.com/Russel88/CRISPRCasTyper), a tool
to predict CRISPR arrays and _cas_ genes in bacterial genome sequences.
This fork aims to replace the core depencencies of CCTyper with more
efficient reimplementations to reduce waiting times 🕰️⬇️ and energy use 🔌🌳.

The main tools that are replaced are:

| **Original** | **Replacement**                                       | **Implemented** |
| ------------ | ----------------------------------------------------- | --------------- |
| Prodigal     | [Pyrodigal](https://pyrodigal.readthedocs.io/)        | ✅              |
| HMMer        | [PyHMMER](https://pyhmmer.readthedocs.io/en/stable/)  | ✅              |
| MinCED       | [Diced](https://diced.readthedocs.io/en/stable/)      | ❌              |

All of these are reimplementations of the original code to run more effeciently.
They should produce identical output while using less resources.
(Big thanks to [@althonos](https://github.com/althonos)!)

## Technical details and benchmarks

The program CCTyper starts by identifying _cas_ genes in your input sequences.
It does this using Prodigal to predict Open Reading Frames (ORFs), and then
screen these using HMMer against a database of known _cas_ genes. These are
actually the most compute intensive steps in the whole process and take
about 20 and 75% of the total runtime, respectively.

Next up is the identification of CRISPR arrays, by looking for direct repeats.
This is done with [MinCED](https://github.com/ctSkennerton/minced),
which is based on the trusty old [CRT](http://www.room220.com/crt/).
This takes only about 1% of the total runtime, so a potential speed up
here will only have a minor effect.

By replacing each of these tools with an optimised replacement,
I hope to improve the general effeciency of the tool and improve
the feasibility of large-scale CRISPR-Cas identification in datasets
of hundreds of thousands of bacterial genomes.

Despite changing the core modules of CCTyper, I want to ensure that the
**output remains identical**. Therefore, benchmarks will evaluate output
files as well as runtime and use of computational resources.

## Acknowledgements

This tool can only exist thanks to the work of:

- Jakob Russel: [CCTyper](https://github.com/Russel88/CRISPRCasTyper)
([paper](https://doi.org/10.1089/crispr.2020.0059),
[freely available preprint](https://doi.org/10.1101/2020.05.15.097824))
- Martin Larralde: [Pyrodigal](https://pyrodigal.readthedocs.io/en/stable/)
([paper](https://doi.org/10.21105/joss.04296)),
[PyHMMER](https://pyhmmer.readthedocs.io/en/stable/)
([paper](https://doi.org/10.1093/bioinformatics/btad214))
- Sean Eddy: [HMMER](https://hmmer.org)
- Hyatt, D., Chen, GL., LoCascio, P.F., _et al._: [Prodigal](https://doi.org/10.1186/1471-2105-11-119)
- Connor T. Skennerton: [MinCED](https://github.com/ctSkennerton/minced)
- Camacho C., Coulouris G., Avagyan V., _et al._: [BLAST](https://doi.org/10.1186/1471-2105-10-421)
