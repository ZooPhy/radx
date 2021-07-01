import logging
import os
import subprocess
import time


def run_cmd(cmd_list, timeout=None):
    ret = None
    # TODO: Subprocess doesn't support the pipe command by default
    if "|" not in cmd_list and ">" not in cmd_list:
        logging.info("--- Running SP '%s'", " ".join(cmd_list))
        ret = subprocess.run(cmd_list, check=True)
        if ret.returncode != 0:
            logging.info("--- Return code '%s'", ret.returncode)
    else:
        logging.info("--- Running SYS '%s'", " ".join(cmd_list))
        ret = os.system(" ".join(cmd_list))
        if ret != 0:
            logging.info("--- Return code '%s'", ret)
    return ret
