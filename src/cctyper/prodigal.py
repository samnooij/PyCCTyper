import os
import subprocess
import logging
import sys
import re

import pandas as pd


class Prodigal(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run_prodigal(self):

        if not self.redo:
            logging.info("Predicting ORFs with Pyrodigal")

            # Run prodigal
            with open(self.out + "prodigal.log", "w") as prodigal_log:
                subprocess.run(
                    [
                        "pyrodigal",
                        "-i",
                        self.fasta,
                        "-a",
                        self.out + "proteins.faa",
                        "-p",
                        self.prod,
                        "-j",
                        str(self.threads),
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=prodigal_log,
                )

            # Check if succesful
            self.check_rerun()

            # Make gene table
            self.get_genes()

    def check_rerun(self):
        # Check prodigal output
        if os.stat(self.prot_path).st_size == 0:
            if self.prodigal == "single":
                logging.warning("Prodigal failed. Trying in metagenome mode")
                self.prodigal = "meta"
                self.run_prodigal()
            else:
                logging.critical("Prodigal failed! Check the log")
                sys.exit()

    def get_genes(self):

        with open(self.out + "genes.tab", "w") as gene_tab:
            subprocess.run(
                ["grep", "^>", self.out + "proteins.faa"], stdout=gene_tab
            )

        genes = pd.read_csv(
            self.out + "genes.tab",
            sep=r"\s+",
            header=None,
            usecols=(0, 2, 4, 6),
            names=("Contig", "Start", "End", "Strand"),
        )

        genes["Contig"] = [re.sub("^>", "", x) for x in genes["Contig"]]
        genes["Pos"] = [int(re.sub(".*_", "", x)) for x in genes["Contig"]]
        genes["Contig"] = [re.sub("_[0-9]*$", "", x) for x in genes["Contig"]]

        self.genes = genes

        genes.to_csv(self.out + "genes.tab", index=False, sep="\t")
