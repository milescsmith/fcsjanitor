# fcsjanitor

Automate cleanup of CyTOF fcs files using the algorithm as outined in
Fluidigm techincal note PN 400248 B1
*"Approach to Bivariate Analysis of Data Acquired Using the Maxpar Direct
Immune Profiling Assay"*

## Install

Requires Python 3.10+

```
pip install git+https://github.com/milescsmith/fcsjanitor
```
## Usage

```
Usage: janitor [OPTIONS]                                                                                                                                                                           
                                                                                                                                                                                                    
 Clean one or more CyTOF FCS files.                                                                                                                                                                 
                                                                                                                                                                                                    
 *  --input    -i      PATH           Either a list of input files or directories in which to look for FCS files
                                      [default: None]
                                      [required]
    --output   -o      PATH           Location to write cleaned FCS files
                                      [default: /Users/milessmith/workspace/fcsjanitor]
    --method   -m      TEXT           Method to use when filtering events
                                      [default: FilterMethod.sd]
    --suffix   -s      TEXT           Suffix to append to the new files
                                      [default: None]
    --format   -f      [anndata|fcs]  Save the cleaned data as an FCS or AnnData object?
                                      [default: OutputFormat.fcs]
    --verbose  -v
    --version
    --help                            Show this message and exit.
```