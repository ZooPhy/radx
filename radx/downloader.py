"""This file contains methods for downloading resources 
   needed for the RADx pipeline automatically.
"""
import time

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

from radx.settings import GISAID_PASSWORD, GISAID_USERNAME


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
