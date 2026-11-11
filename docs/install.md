---
layout: default
title: "Installation"
nav_order: 2
---

# Installation

```bash
pip install git+https://github.com/yucj-520/altharm.git
```

or, from a copy of the repository (needed to run the tests and the stored examples),

```bash
git clone https://github.com/yucj-520/altharm.git
cd altharm
pip install ".[dev]"
pytest -q
```

Requirements: Python 3.10 or newer. NumPy, SciPy, spglib and ASE are installed
automatically. The `dev` extra adds pytest and the `plot` extra adds Matplotlib.

Check the installation:

```bash
altharm --help
```

## Quantum ESPRESSO

`altharm` reads the `data-file-schema.xml` file of a finished `pw.x` run. It was developed
and tested with Quantum ESPRESSO 7. No change to Quantum ESPRESSO is needed.

The optional executable `harmwave.x`, which does the same analysis inside Quantum
ESPRESSO, is described in
[`qe_native/README.md`](https://github.com/yucj-520/altharm/blob/master/qe_native/README.md).
