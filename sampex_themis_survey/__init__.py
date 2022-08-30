import warnings
import pathlib
import configparser

__version__ = '0.0.1'

# Load the configuration settings.
here = pathlib.Path(__file__).parent.resolve()
settings = configparser.ConfigParser()
settings.read(here / 'config.ini')

config = {'code_dir': here}

# # Go here if config.ini exists (don't crash if the project is not yet configured.)
# if 'Paths' in settings:  
#     try:
#         data_dir = settings['Paths']['data_dir']
#         config = {'code_dir': here, 'data_dir': data_dir}
#     except KeyError as err:
#         warnings.warn('sampex_themis_survey did not find the config.ini file. '
#             'Did you run "python3 -m sampex_themis_survey config"?')
            
# warnings.warn('sampex_themis_survey did not find the config.ini file. '
#             'Did you run "python3 -m sampex_themis_survey config"?')
    