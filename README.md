# Introduction
This is the package containing the sampex_themis_survey code.


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
A lot of data science projects load data external to the source code (a good practice) so `project/__main__.py` creates a `project/config.ini` file that is loaded on import by `project/__init__.py`. 

To execute `project/__main__.py`, first install `project` with the steps above and then run `python3 -m project config` and answer the prompt. You will now see a `project/config.ini` with two paths: one to this project and the other to the specified data directory. One loaded by `project/__init__.py`, this dictionary is accessed via `import project.config`.