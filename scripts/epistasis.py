"""
Python dependencies: pandas
    - sudo pip install pandas

The script take as input a matrix whit n sampples and p features(PCAs, SNPs and phenotype) and perform
epistatic anaylis.

The epistatis analysis use two SNPs. With the two SNPs we create two logistic regression
models. The first one takes PCAs and the two SNPs, the second one takes PCAs, SNPs and the
interaction term of the two SNPs. Afterwards we compare the two model using likelihood ratio
test.

samples          PCAs      PHENO   21:14603577:C:T:b37   21:14604174:C:T:b37   21:14604361:C:T:b37
WTCCC125636      1.029     case           2                    2                      0.3
WTCCC125637      1.05      control        1                    2                       2
WTCCC125638      0.997     control        0                    1                       1
WTCCC125639      1.5       case           1                    1.9                     5
...                 ..                     ..                   ..
..                  .                      .                    .


It is possible to apply filters on MAF and Rsq. The outputs are going

How to run on the command line:
$python create_pheno+dose.py -h

"""

import pandas as pd
import numpy as np
import os
import argparse
from likelihood_ratio_test import LRT
from scipy.stats import chi2
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import log_loss


class epistasis():
    """
    This class contain functions to run SNP-SNP epistasis interactions giving
    a phenotype+dosage matrix processed byt create_pheno+dose.py script.
    """

    def __init__(self):
        """
        Initialise the class
        """

    def get_snps_loops(self, chr, loop):
        """
        This function takes the chromosome and loop number and
        store SNPs from the two loop anchors in two lists, A and B

        Parameters
        ----------
        chr: Chromosome number
        loop: Loop number

        Return
        ------
        Both lists
        """
        loopA = []
        loopB = []
        with open("../tests/ld/chr"+chr+"/WTCCC.chr"+str(chr)+".hrc.imputed.plink.dose."
                  "info.0.6.hardcalls.loop"+str(loop)+".A.prune.in", "r") as file_A:
            for line in file_A:
                line_hndl = line.rstrip()
                loopA.append(line_hndl)

        with open("../tests/ld/chr"+chr+"/WTCCC.chr"+str(chr)+".hrc.imputed.plink.dose."
                  "info.0.6.hardcalls.loop"+str(loop)+".B.prune.in", "r") as file_B:
            for line in file_B:
                line_hndl = line.rstrip()
                loopB.append(line_hndl)

        return loopA, loopB

    def logit(self, matrix, loopA, loopB):
        """
        This functions takes a matrix containing PCAs, SNPs and phenotype and
        perform the logistic regression for the two models. The first one "aditive model"
        contain just the values for the two SNPs, the secong model the "interaction model"
        contin the two SNPs plus the interaction term. for every of the model

        Aguments
        --------
        matrix: str
            matrix processed byt create_pheno+dose.py script
        snp1: str
            name of the SNP from anchor A
        snp2: str
            name of the SNP from anchor B
        """
        df = pd.read_csv(MATRIX, sep="\t")

        for SNPa in loopA:
            for SNPb in loopB:

                model = SGDClassifier(loss="log",
                                      penalty="l2")
                                      #n_iter=needed_sgd_iter(N_SAMPLES))


                features_null = np.array(df[["PC1", "PC2","PC3",
                                             "PC4", "PC5", "PC6",
                                              "PC7", "PC8", "PC9",
                                              "PC10", SNPa, SNPb]])

                print(features_null)
                break
            break
            """
                features_alternate = fd[["PC1", "PC2","PC3",
                                             "PC4", "PC5", "PC6",
                                              "PC7", "PC8", "PC9",
                                              "PC10", SNPa, SNPb]]
            """

            #features_alternate, labels, lr_model, features_null=None):








if __name__ == "__main__":

    #set up the command line parameters
    PARSER = argparse.ArgumentParser(
        description="""
                    The script take as input a matrix whit n sampples
                    and p features(PCAs, SNPs and phenotype) and perform
                    epistatic anaylis.

                    The epistatis analysis use two SNPs. With the two SNPs we create two logistic regression
                    models. The first one takes PCAs and the two SNPs, the second one takes PCAs, SNPs and the
                    interaction term of the two SNPs. Afterwards we compare the two model using likelihood ratio
                    test. """)

    PARSER.add_argument(
        "--matrix",
        help="Matrix generated by create_pheno+dose.py",
        type=str)

    PARSER.add_argument(
        "--chr", help="chromosome",
        type=str)

    PARSER.add_argument(
        "--loop", help="chromatine loop",
        type=str)


    #Get matching parameters from the command line
    ARGS = PARSER.parse_args()

    MATRIX = ARGS.matrix
    CHR = ARGS.chr
    LOOP = ARGS.loop

    #Read phenotype dosages
    loopA, loopB = epistasis().get_snps_loops(CHR, LOOP)
    print(loopA, loopB)

    epistasis().logit(MATRIX, loopA, loopB)

