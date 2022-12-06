import warnings
import pathlib
import configparser

__version__ = '0.0.1'

# Load the configuration settings.
here = pathlib.Path(__file__).parent.resolve()
settings = configparser.ConfigParser()
settings.read(here / 'config.ini')

config = {
    'code_dir': here, 
    'data_dir': here.parent / 'data',
    'plots_dir': here.parent / 'plots'
    }