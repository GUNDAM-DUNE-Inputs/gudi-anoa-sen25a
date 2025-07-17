# Setup with python venv

Create a local virtual environment with

```bash
python3 -m venv --system-site-packages .
```

Start the virtual environment using 

```bash
. bin/activate
```
Leave the virtual environment using

```bash
deactivate
```

Install science plots

```bash
pip install scienceplots
```

Create your scripts using

```python
#!/usr/bin/env python3
# Make sure we get the "local" version of python3 so our paths are right

print("hello")
```

# Generate the production height histogram.

Make sure that ROOT is in your path and that you are in the virtual environment (run `source bin/activate`), and then

```bash
./makeProductionHeights.py
```

This can take long time!  Running it without arguments will create our
default binning (Note to us: be sure and update this if the default
changes).  Check the source if you want to change the binning (there are
command line arguments for that).

This will make a file named `ProdHeightDist_E<eBINS>_cZ<zBINS>_N1_hbin2.5.root` that will need to be copied into `<top-level>/inputs/parameters/oscProb/`.

# Generate the grid for neutrino energy and zenith cosine

Make sure that ROOT is in your path, and that you are in the virtual
environment (run `source bin/activate`).  You only want to use this if you
are NOT using production height averaging (if you are averaging, then the
production height histsogram will provide the grid).  Generate the grid
using

```bash
./makeEnergyZenith.py
```

This is rather quick.  Running it without arguments will create our default
binning (Note to us: be sure and update this if the default changes).  
