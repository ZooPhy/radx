"""This file contains methods for downloading resources 
   needed for the RADx pipeline automatically.
"""
import time
import os
import tarfile

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

from radx.settings import (GISAID_PASSWORD, GISAID_USERNAME,
                           PATH_TO_GISAID, PATH_TO_DOWNLOADS)


class GISAIDDownloader(object):
    def __init__(self, driver="chrome"):
        self.uname = GISAID_USERNAME
        self.passwd = GISAID_PASSWORD
        self.driver = webdriver.Chrome() if driver=="chrome" else webdriver.Firefox()

    def dump_gisaid_data(self):
        self.login()
        self.download_fasta()
        self.download_metadata()
        self.wait_for_download()
        self.shutdown_driver()
        self.move_metadata_file()
        self.move_sequence_file()

    def login(self):
        self.driver.get("https://www.epicov.org/epi3/start")
        assert "GISAID" in self.driver.title

        unamefld = self.driver.find_element_by_id("elogin")
        unamefld.clear()
        unamefld.send_keys(self.uname)

        passwdfld = self.driver.find_element_by_id("epassword")
        passwdfld.clear()
        passwdfld.send_keys(self.passwd)

        passwdfld.send_keys(Keys.RETURN)
        time.sleep(3)

    def download_fasta(self):
        self.driver.find_element(By.XPATH, '//div[text()="Downloads"]').click()
        time.sleep(3)
        iframe = self.driver.find_element_by_xpath("//iframe[@src='about:blank']")
        self.driver.switch_to.frame(iframe)
        self.driver.find_element(By.XPATH, '//div[text()="FASTA"]').click()
        time.sleep(2)
        iframe = self.driver.find_element_by_xpath("//iframe[@src='about:blank']")
        self.driver.switch_to.frame(iframe)
        self.driver.find_element(By.XPATH, '//input[@type="checkbox"]').click()
        time.sleep(2)
        self.driver.find_element(By.XPATH, '//button[text()="Download"]').click()
        time.sleep(5)

    def download_metadata(self):
        iframe = self.driver.find_element_by_xpath("//iframe[@src='about:blank']")
        self.driver.switch_to.frame(iframe)
        self.driver.find_element(By.XPATH, '//div[text()="metadata"]').click()
        time.sleep(2)
        iframe = self.driver.find_element_by_xpath("//iframe[@src='about:blank']")
        self.driver.switch_to.frame(iframe)
        self.driver.find_element(By.XPATH, '//input[@type="checkbox"]').click()
        time.sleep(2)
        self.driver.find_element(By.XPATH, '//button[text()="Download"]').click()

    def wait_for_download(self):
        time.sleep(20)

    def shutdown_driver(self):
        self.driver.close()

    def move_metadata_file(self):
        if PATH_TO_DOWNLOADS:
            down_dir = PATH_TO_DOWNLOADS
        else:
            from pathlib import Path
            down_dir = str(Path.home() / "Downloads") 
        work_dir = PATH_TO_GISAID
        metadata_files = [x for x in os.listdir(down_dir) if x[:13]=="metadata_tsv_" and x[-6:]=="tar.xz"]
        metadata_arxiv = sorted(metadata_files, reverse=True)[0]
        src_file = os.path.join(down_dir, metadata_arxiv)
        dst_file = os.path.join(work_dir, metadata_arxiv)
        # delete if the archive file and directory exists
        if os.path.exists(dst_file):
            os.rmdir(dst_file)
        if os.path.exists(dst_file[:-6]):
            os.rmdir(dst_file[:-6])
        # move the file and extract
        os.rename(src_file, dst_file)
        if os.path.exists(dst_file):
            tar = tarfile.open(dst_file, "r:xz")
            tar.extractall(work_dir)
            tar.close()

    def move_sequence_file(self):
        if PATH_TO_DOWNLOADS:
            down_dir = PATH_TO_DOWNLOADS
        else:
            from pathlib import Path
            down_dir = str(Path.home() / "Downloads") 
        work_dir = PATH_TO_GISAID
        metadata_files = [x for x in os.listdir(down_dir) if x[:15]=="sequence_fasta_" and x[-6:]=="tar.xz"]
        metadata_arxiv = sorted(metadata_files, reverse=True)[0]
        src_file = os.path.join(down_dir, metadata_arxiv)
        dst_file = os.path.join(work_dir, metadata_arxiv)
        # delete if the archive file and directory exists
        if os.path.exists(dst_file):
            os.rmdir(dst_file)
        if os.path.exists(dst_file[:-6]):
            os.rmdir(dst_file[:-6])
        # move the file and extract
        os.rename(src_file, dst_file)
        if os.path.exists(dst_file):
            tar = tarfile.open(dst_file, "r:xz")
            tar.extractall(work_dir)
            tar.close()
