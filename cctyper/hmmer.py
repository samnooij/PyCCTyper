import os
import sys
import logging
import pyhmmer
import re
import tqdm
import collections

import pandas as pd
from pathlib import Path


class HMMER(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def main_hmm(self):
        # If redo just get the table
        if self.redo:
            self.read_hmm()
        # Else run HMMER load and write data
        else:
            hmm_df = self.run_hmm()
            self.load_hmm(dataframe=hmm_df)
            self.write_hmm()

        # Check if any cas genes
        self.check_hmm()

        # Parse
        self.parse_hmm()

    # Run pyHMMER and parse required information
    def hmmsearch(self, progress=bool):

        hmms = []
        hmm_files = list(Path(self.pdir).glob("*.hmm"))
        for hmm_file in hmm_files:
            with pyhmmer.plan7.HMMFile(hmm_file) as hmmfile:
                hmm = hmmfile.read()
            hmms.append(hmm)

        Result = collections.namedtuple(
            "Result",
            [
                "target_name",
                "target_accession",
                "tlen",
                "query_name",
                "query_accession",
                "qlen",
                "evalue",
                "bitscore",
                "bias",
                "domain_number",
                "out_of_domains",
                "c_evalue",
                "i_evalue",
                "domain_bitscore",
                "domain_bias",
                "hmm_from",
                "hmm_to",
                "ali_from",
                "ali_to",
                "env_from",
                "env_to",
                "posterior_probabilities",
                "description",
            ],
        )

        result_list = []

        with pyhmmer.easel.SequenceFile(
            self.prot_path, digital=True
        ) as seqs_file:
            sequences = seqs_file.read_block()

        def collect_results(hits):
            hmm_name = hits.query.name
            for hit in hits:
                if hit.included:
                    orf_name = hit.name
                    evalue = hit.evalue
                    bitscore = hit.score
                    domains = len(hit.domains)
                    for i in range(domains):
                        domain_number = i + 1  # Because it counts 0-based
                        target_length = hit.domains[i].alignment.target_length
                        hmm_length = hit.domains[i].alignment.hmm_length

                        hmm_start = hit.domains[i].alignment.hmm_from
                        hmm_end = hit.domains[i].alignment.hmm_to

                        ali_start = hit.domains[i].alignment.target_from
                        ali_end = hit.domains[i].alignment.target_to

                        env_start = hit.domains[i].env_from
                        env_end = hit.domains[i].env_to

                        # PyHMMER does not record strand by default; only in the 'long targets pipeline', which does
                        # not support the cpus= option!
                        # However, Pyrodigal saves this in the FASTA ID, which is stored as hit description!
                        strand = hit.description.split("#")[3].strip(" ")
                        if strand is None:
                            strand = 0

                        result_list.append(
                            Result(
                                orf_name,
                                hit.accession,
                                target_length,
                                hmm_name,
                                hit.domains[i].alignment.hmm_accession,
                                hmm_length,
                                evalue,
                                bitscore,
                                hit.bias,
                                domain_number,
                                domains,
                                hit.domains[i].c_evalue,
                                hit.domains[i].i_evalue,
                                hit.domains[i].score,
                                hit.domains[i].bias,
                                hmm_start,
                                hmm_end,
                                ali_start,
                                ali_end,
                                env_start,
                                env_end,
                                hit.domains[
                                    i
                                ].alignment.posterior_probabilities,
                                hit.description,
                            )
                        )

        if progress:
            for hits in tqdm.tqdm(
                pyhmmer.hmmer.hmmsearch(hmms, sequences, cpus=self.threads),
                total=len(hmms),
            ):
                collect_results(hits=hits)
        else:
            for hits in pyhmmer.hmmer.hmmsearch(
                hmms, sequences, cpus=self.threads
            ):
                collect_results(hits=hits)

        result_df = pd.DataFrame(result_list, columns=Result._fields)

        return result_df

    # Parallel search of all HMMs
    def run_hmm(self):

        logging.info("Running pyHMMER against Cas profiles")

        # Make dir
        os.mkdir(self.out + "hmmer")
        # Each HMM
        if self.lvl == "DEBUG" or self.simplelog:
            hmm_df = self.hmmsearch(progress=False)
        else:
            hmm_df = self.hmmsearch(progress=True)

        logging.info("Write pyHMMER output to file")
        hmm_df.to_csv(
            os.path.join(self.out + "hmmer", "Cas_HMMer-like.tab"),
            sep="\t",
            index=False,
        )

        return hmm_df

    # Load data
    def load_hmm(self, dataframe):

        logging.debug("Loading HMMER output")

        # Load relevant columns from pyHMMER output
        hmm_df = dataframe.loc[
            :,
            [
                "query_name",
                "target_name",
                "tlen",
                "qlen",
                "evalue",
                "bitscore",
                "hmm_from",
                "hmm_to",
                "ali_from",
                "ali_to",
                "env_from",
                "env_to",
                "description",
            ],
        ]
        # Rename some columns to match CCTyper's default output
        hmm_df = hmm_df.rename(
            columns={
                "query_name": "Hmm",
                "target_name": "ORF",
                "evalue": "Eval",
                "bitscore": "score",
            }
        )

        # Split the 'description' field to find the start and end positions, and strand
        hmm_df[["nothing", "start", "end", "strand", "info"]] = (
            hmm_df.description.str.split("#", expand=True)
        )
        # Note that these are parsed as strings, with whitespace surrounding them...
        hmm_df["start"] = hmm_df["start"].str.strip().astype(int)
        hmm_df["end"] = hmm_df["end"].str.strip().astype(int)
        hmm_df["strand"] = hmm_df["strand"].str.strip().astype(int)
        # And remove the remaining, unused columns
        hmm_df.drop(["nothing", "info"], axis=1)

        # Add columns
        hmm_df["Acc"] = [re.sub("_[0-9]*$", "", x) for x in hmm_df["ORF"]]
        hmm_df["Pos"] = [int(re.sub(".*_", "", x)) for x in hmm_df["ORF"]]

        # Coverages of aligments
        def covs(df_sub):
            df_sub["Cov_seq"] = (
                len(
                    set(
                        [
                            x
                            for sublst in [
                                list(range(i, j))
                                for i, j in zip(
                                    df_sub["ali_from"], df_sub["ali_to"] + 1
                                )
                            ]
                            for x in sublst
                        ]
                    )
                )
                / df_sub["tlen"]
            )
            df_sub["Cov_hmm"] = (
                len(
                    set(
                        [
                            x
                            for sublst in [
                                list(range(i, j))
                                for i, j in zip(
                                    df_sub["hmm_from"], df_sub["hmm_to"] + 1
                                )
                            ]
                            for x in sublst
                        ]
                    )
                )
                / df_sub["qlen"]
            )
            df_sub = df_sub[
                [
                    "Hmm",
                    "ORF",
                    "tlen",
                    "qlen",
                    "Eval",
                    "score",
                    "start",
                    "end",
                    "Acc",
                    "Pos",
                    "Cov_seq",
                    "Cov_hmm",
                    "strand",
                ]
            ]
            df_sub = df_sub.drop_duplicates()
            return df_sub

        hmm_df = hmm_df.groupby(["Hmm", "ORF"]).apply(covs)
        hmm_df.reset_index(drop=True, inplace=True)
        self.hmm_df = hmm_df.drop_duplicates()

    # Write to file
    def write_hmm(self):
        self.hmm_df.to_csv(self.out + "hmmer.tab", sep="\t", index=False)

    # Read from file
    def read_hmm(self):
        try:
            self.hmm_df = pd.read_csv(self.out + "hmmer.tab", sep="\t")
        except Exception:
            logging.error("No matches to Cas HMMs")
            sys.exit()

    # Check if any cas genes
    def check_hmm(self):
        if len(self.hmm_df) == 0:
            logging.info("No Cas proteins found.")
        else:
            self.any_cas = True

    # Parse
    def parse_hmm(self):

        if self.any_cas:

            logging.debug("Parsing HMMER output")

            # Pick best hit
            self.hmm_df.sort_values("score", ascending=False, inplace=True)
            self.hmm_df.drop_duplicates("ORF", inplace=True)
