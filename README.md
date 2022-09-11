# Introduction
Analyze the double conjunctions between the SAMPEX satellite and the THEMIS all-sky imager (ASI) array, and also triple conjunctions between SAMPEX,THEMIS ASIs, and THEMIS satellites.

# Dependencies
Before you can run many of these scripts, you will need to install [IRBEM](https://github.com/PRBEM/IRBEM). See their installation instructions for more information.

# Installation
- To install as a developer run:
  ```bash
    git clone git@github.com:mshumko/sampex_themis_survey.git
    cd sampex_themis_survey

    # Then one of these (see comment in requirement.txt):
    python3 -m pip install -e .
    ```
    or 
    ```bash
    python3 -m pip install -r requirements.txt 
    ```

# Configuration
- If you already saved the sampex HILT and attitude files on your computer, and in a unique location, run `python3 -m sampex init` and paste the data directory there.
- Similarly, if you already have some of the THEMIS ASI data on your computer in a unique location, run `python3 -m asilib init` and paste the data directory there.