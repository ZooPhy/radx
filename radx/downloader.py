"""This file contains methods for downloading resources 
   needed for the RADx pipeline automatically.
"""
import json
import logging
import os
import re
import tarfile
import time

import numpy as np
import pandas as pd
import requests
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

from radx.settings import (GISAID_PASSWORD, GISAID_USERNAME, PATH_TO_DOWNLOADS,
                           PATH_TO_GISAID, PATH_TO_MUTATIONS, PATH_TO_ALCOV)


class GISAIDDownloader(object):
    def __init__(self, driver="chrome"):
        self.uname = GISAID_USERNAME
        self.passwd = GISAID_PASSWORD
        self.driver = webdriver.Chrome() if driver=="chrome" else webdriver.Firefox()

    def dump_gisaid_data(self):
        self.login()
        self.download_fasta()
        self.wait_for_download()
        self.shutdown_driver()
        self.move_sequence_file()

    def dump_gisaid_metadata(self):
        self.login()
        self.download_metadata()
        self.wait_for_download()
        self.shutdown_driver()
        self.move_metadata_file()

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
        self.driver.find_element(By.XPATH, '//div[text()="Downloads"]').click()
        time.sleep(3)
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
        time.sleep(90)

    def shutdown_driver(self):
        self.driver.close()

    def move_metadata_file(self):
        if PATH_TO_DOWNLOADS:
            down_dir = PATH_TO_DOWNLOADS
        else:
            from pathlib import Path
            down_dir = str(Path.home() / "Downloads") 
        work_dir = PATH_TO_GISAID
        try:
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
        except PermissionError as _:
            logging.info("Permission error in path: %s. Please move files manually.", down_dir)

    def move_sequence_file(self):
        if PATH_TO_DOWNLOADS:
            down_dir = PATH_TO_DOWNLOADS
        else:
            from pathlib import Path
            down_dir = str(Path.home() / "Downloads") 
        work_dir = PATH_TO_GISAID
        try:
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
        except PermissionError as _:
            logging.info("Permission error in path: %s. Please move files manually.", down_dir)

class VariantDownloader(object):
    def __init__(self):
        self.descendents = {}
        self.voc = []
        self.voi = []
        self.afm = []
        self.mutations = {}

    def download(self):
        self.get_descendent_lineages()
        self.get_voc_voi()
        self.get_mutations(self.voc+self.voi)
        self.write_mutations()

    def get_voc_voi(self):
        url = 'https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/'
        html = requests.get(url).content
        df_list = pd.read_html(html)
        voc_table, voi_table, afm_table = df_list[0], df_list[1], df_list[2]
        self.voc = [re.sub('[#§ ]', '', x) for x in list(voc_table["Pango lineage*"])]
        self.voi = [re.sub('[#§ ]', '', x) for x in list(voi_table["Pango lineage*"])]
        self.afm = [re.sub('[#§ ]', '', x) for x in list(afm_table["Pango lineage*"])]
        # In case we need to process the WHO labels
        # print("VOC", list(voc["Pango lineage*"]), "\t", list(voc["WHO\xa0label"]))
        # print("VOI", list(voi["Pango lineage*"]), "\t", list(voi["WHO\xa0label"]))

    def get_descendent_lineages(self):
        self.descendents = {}
        url = 'https://cov-lineages.org/lineage_list.html'
        html = requests.get(url).content
        df_list = pd.read_html(html)
        df = df_list[0]
        print("Rows:", len(df.index), "Columns:", list(df.columns))
        df["AliasOf"] = df["Description"].str.extract(r'(?:(?:Alias of )|(?:Previously ))([A-Z0-9.]+),', expand=True)
        # Get parents based on Aliases. ASL = AliasSubLevel , AP = AliasParent
        df["ASL"] = df["AliasOf"].str.count(r"\.").fillna(0.0).astype(int)
        df["AP"] = df["AliasOf"].str.split(r"\.")
        df["AP"] = df.apply(lambda x: '.'.join(x["AP"][0:-1]) if x["ASL"]>1 else np.nan, axis=1)
        # Get parents based on NomenClature. DSL = DirectSubLevel , DP = DirectParent
        df["DSL"] = df["Lineage"].str.count(r"\.").fillna(0.0).astype(int)
        df["DPS"] = df["Lineage"].str.split(r"\.")
        df["DP"] = df.apply(lambda x: '.'.join(x["DPS"][0:-1]) if x["DSL"]>1 else np.nan, axis=1)
        # Get dictionary of lineages
        lineages = set(list(df["Lineage"]))
        print("Number of valid lineages:", len(lineages), "\n")
        df['Parent'] = df.apply(lambda x: x.DP if x.DP in lineages else (x.AP if x.AP in lineages else np.nan) , axis=1)
        parents = dict(zip(df.Lineage, df.Parent))
        # Establish lineage hierarchy, GP : GrandParent, GGP: GreatGrandParent and so on
        df["GP"] = df.apply(lambda x: parents[x["Parent"]] if x["Parent"] in parents else np.nan, axis=1)
        df["GGP"] = df.apply(lambda x: parents[x["GP"]] if x["GP"] in parents else np.nan, axis=1)
        df["GGGP"] = df.apply(lambda x: parents[x["GGP"]] if x["GGP"] in parents else np.nan, axis=1)
        df["GGGGP"] = df.apply(lambda x: parents[x["GGGP"]] if x["GGGP"] in parents else np.nan, axis=1)
        df["GGGGGP"] = df.apply(lambda x: parents[x["GGGGP"]] if x["GGGGP"] in parents else np.nan, axis=1)
        # Check if withdrawn
        df["Withdrawn"] = df["Description"].str.contains('Withdrawn:', regex=False)
        # Create dictionary of children
        for _, row in df.iterrows():
            lin = row["Lineage"]
            c = (list(df[df["Parent"].isin([lin])]["Lineage"]) + list(df[df["GP"].isin([lin])]["Lineage"])
                + list(df[df["GGP"].isin([lin])]["Lineage"]) + list(df[df["GGGP"].isin([lin])]["Lineage"])
                + list(df[df["GGGGP"].isin([lin])]["Lineage"]) + list(df[df["GGGGGP"].isin([lin])]["Lineage"]))
            # get unique children and sort
            children = sorted(list(set(c)))
            if len(children)>0:
                self.descendents[lin] = children
        df["Children"] = df.apply(lambda x: ",".join(self.descendents[x["Lineage"]]) if x["Lineage"] in self.descendents else "", axis=1)

    def get_mutations(self, pango_lineages, batch_size=10, frequency=0.8):
        lineage_mutations = {}
        self.mutations = []
        for _, pango_lineage in enumerate(pango_lineages):
            print("Processing:", pango_lineage)
            if pango_lineage in self.descendents:
                lin_desc = self.descendents[pango_lineage]
            else:
                print("ERROR: No descendents found for", pango_lineage)
                lin_desc = [pango_lineage]
            print("Processing", pango_lineage, "with", len(lin_desc), "descendents")
            # Process in batches
            for i in range(0, len(lin_desc), batch_size):
                batch = lin_desc[i:i + batch_size]
                print("Getting mutations for", batch)
                url = "https://api.outbreak.info/genomics/lineage-mutations?pangolin_lineage="+",".join(batch)+"&frequency="+str(frequency)
                response = requests.get(url)
                muts_res = response.json()["results"]
                for lin, muts in muts_res.items():
                    lineage_mutations[lin] = muts
                    print(lin, "has", len(muts), "mutations with", frequency, "frequency")
                    with open(os.path.join(PATH_TO_MUTATIONS, lin+".json"), "w") as ofile:
                        print(json.dumps(muts, indent=4, sort_keys=True), file=ofile)
                    for mut in muts:
                        self.mutations.append(mut)

    # CREDITS for method : ALCOV - https://github.com/Ellmen/alcov
    def fix_mut_name(self, old_mut_name):
        mut_name = old_mut_name.upper()
        if '/' in mut_name:
            mut_name = mut_name[:mut_name.find('/')]
        if not mut_name.startswith('ORF'):
            return mut_name
        else:
            col_idx = mut_name.find(':')
            return "ORF" + mut_name[3:col_idx].lower() + mut_name[col_idx:]

    # CREDITS for method : ALCOV - https://github.com/Ellmen/alcov
    def write_mutations(self, filename="mutations.json"):
        muts = list(set([self.fix_mut_name(m['mutation']) for m in self.mutations]))
        lins = list(set([m['lineage'].upper() for m in self.mutations]))
        mut_lins = {mut: {lin: 0 for lin in lins} for mut in muts}
        for raw_m in self.mutations:
            mut = self.fix_mut_name(raw_m['mutation'])
            lin = raw_m['lineage'].upper()
            prev = raw_m['prevalence']
            mut_lins[mut][lin] = prev
        with open(os.path.join(PATH_TO_MUTATIONS, filename), "w") as ofile:
            print(json.dumps(mut_lins, indent=4, sort_keys=True), file=ofile)
        # write the python file to the alcov directory for loading
        with open(os.path.join(PATH_TO_ALCOV, 'alcov', 'mutations.py'), "w") as ofile:
            ofile.write('mutations = {}'.format(mut_lins))
