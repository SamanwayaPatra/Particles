import logging 
import os

LOG_FILE_NAME = "main.log"
if os.path.exists(LOG_FILE_NAME):
  os.remove(LOG_FILE_NAME)
logger = logging.getLogger('MainLogger')
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(LOG_FILE_NAME)
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(fmt)
logger.addHandler(fh)

