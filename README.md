# rnmrfit


## Overview

rnmrfit is a command-line software package for fitting NMR data using customizable YAML configuration files. It processes input data and outputs results via a combination of CSV data tables and interactive HTML plots.

The original R package version of this code is no longer supported but can still be accessed as release `v1.0.0`.

The underlying approach is described in this article: https://doi.org/10.1016/j.jmr.2018.11.004, with a new manuscript under preparation detailing more recent changes.

Please note that both the software and this documentation is under active development and may change with limited notice.

## Usage

Example configuration YAML files are provided in the `examples` directory with details below.
  
**On Windows:**

Run using Command Prompt:

```cmd
rnmrfit.exe <path_to_yaml_file>
```

Alternatively, create a batch script in the same location as the `rnmrfit.exe` and configuration files with the following:

```bat
start powershell -command "& '.\rnmrfit.exe' *.yaml -o fit"
```

Save as something like `fit_all.bat` and double-click to run. See `Examples` section below for more details.

**On Linux:**
 
Run in terminal:

```bash
./rnmrfit <path_to_yaml_file>
```


## Examples 

Two examples are provided in the `examples` directory, one for `2H` and one for `2H` data. The following is based on the `2H` example, but `2H` is largely equivalent.


### Files
  
- `simple.yaml`

    An example set of a barebones configuration file that results in reasonable fit precision for the provided data.

- `optimized.yaml`

    An example set of "optimized" configuration file that results in best fit precision for the provided data.

- `fit_all.bat`
  
    An optional batch script for Windows users to process all `.yaml` files in the same directory as the script with a single double-click. The script uses PowerShell to execute the NMRfit software (outputting all results in a new directory called `fit`):
  
    ```bat
    start powershell -command "& '.\rnmrfit.exe' *.yaml -o fit"
    ```

    Just double-click the script to run (or run `rnmrfit.exe` directly via the PowerShell).

- `2H-VAN-1` - `2H-VAN-5`

    Five sets of Bruker NMR data, each containing 3 experiments.

### Simple YAML configuration

Configuration files are composed of multiple instructions separated by YAML page breaks "---". The following instructions are mandatory:

**Loading data**

The following loads all experiments ("*") in "./2H-VAN-1" as dataset "VAN-1" using Bruker processed directory "1" under "procs".

```yaml
---
dataset: "VAN-1"
path: "./2H-VAN-1"

experiments: "*"
processing: "1"
```

**Defining peaks**

Peaks are grouped under species and resonances. The following lists eight singlets, but other multiplets are possible, e.g. "2 dd 2 4" to specify a double doublet at 2 ppm with coupling constants of 2 and 4 Hz.

```yaml
---
species: "vanillin"
resonances:
  "9.82 s": "H1"
  "8.29 s": "OH2"
  "7.45 s": "H3"
  "7.05 s": "H4"
  "3.91 s": "H5"
```

**Running fit**

The minimum set of options necessary to run a fit is a name for debugging purposes, a ppm range, and a list of species to fit:

```yaml
fit: "2H"
ppm: [50, 200]
species: ["vanillin"]
```

**Output**

Although not strictly mandatory to run `rnmrfit`, no output is generated unless at least one of the following is added.

`plot` generates an interactive plot saved to an HTML file:

```yaml
---
plot: "2H"
ppm: [50, 200]
experiments: "*"
species: ["vanillin"]
```

`calculate_areas` calculates areas and saves them to a CSV file:

```yaml
---
calculate_areas: "areas"
sum_by: resonance
```

`list_parameters` lists fit parameters like position, height, and width:

```yaml
---
list_parameters: "parameters"
```

## FAQs

### Why does the software not run after modifying the YAML settings?
The YAML file structure is sensitive to formatting. Ensure there are no extra spaces, incorrect indentation, or unintended characters. Pay special attention to separators (`---`) between sections; having multiple or misplaced separators can cause parsing errors. Note that most errors are logged to `log.txt` with sufficient detail for correction.


### Why am I getting an error when running the batch script on Windows?
This can happen if the path to your configuration file is incorrect or contains special characters. Ensure the configuration file is in the same directory as `rnmrfit.exe`, and use relative paths without spaces or special characters.


### Can I process multiple datasets simultaneously?
Yes, the `fit_all.bat` script provided under the `examples` directory processes all `.yaml` files in the same directory. Ensure each `.yaml` file is correctly configured.


### Why are my output files empty or missing key data?
This is generally result of an error with the most common culprits:
- incorrect or non-existent dataset paths
- fitting parameters or bounds that are too restrictive, causing the fit process to fail
- general syntax errors in the configuration including extra/blank separators (`---`)

Double-check your input file paths, parameters, and resonance settings in your configuration.


### Why are my baseline corrections or phase adjustments not working as expected?
The baseline and phase parameters in the YAML file must be fine-tuned for your specific data. Adjust the `span`, `polynomial degree`, and `phase order` settings to suit the noise level and signal profile of your dataset. 


### Where is the rest of the documentation
The software is still under active development and documentation will remain limited until features are fully tested and no immediate changes are expected.
