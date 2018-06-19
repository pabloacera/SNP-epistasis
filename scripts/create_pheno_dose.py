"""
Python dependencies: pandas
    - sudo pip install pandas

The script creates a matrix with phenotypes and dosages for all samples and SNPs.

samples     21:14603577:C:T:b37   21:14604174:C:T:b37   21:14604361:C:T:b37   PCAs  PHENO
WTCCC125636        1.029                   2                    2              0.3  Case
WTCCC125637        1.05                    1                    2              2    Control
WTCCC125638        0.997                   0                    1              1    Case
WTCCC125639        1.5                     1                    1.9            5    Case
...                 ..                     ..                   ..
..                  .                      .                    .


It is possible to apply filters on MAF and Rsq. The outputs are going

How to run on the command line:
$python create_pheno+dose.py -h

"""

import pandas as pd
import sys
import os
import argparse


class create_pheno_dosage():

    def __init__(self):

        """
        Initiate the class
        """

    def preprocess_matrix(self, MAF, Rsq, PHENO, DOSE,
                          output="../pheno+dose/"+output):
        """
        This function takes phenotype and dosage files and
        merge them into a single matrix to do the epistasis analysis

        Arguments:
        ----------
        MAF:
            MAF Frecuency to filter SNPs.
            The number introduced will be the minimal SNPs
            frecuency allowed in the files, including that number
        Rsq:
            Imputation quality to filter SNPs.
            The number introduced will be the minimal SNPs
            imputation quality allowed in the files, including that number
        PHENO:
            complete path to phenotype file

        DOSE:
            complete path to dosage file

        Return
        ------
        Print the output Matrix dosage+phenotype into the
        dosage+phenotype folder if no output path is specifies
        """

        if "/" in PHENO:
            output_pheno = PHENO.split("/")[-1]
        else:
            output_pheno = PHENO

        if "/" in DOSE:
            output_dose = DOSE.split("/")[-1]
        else:
            output_dose = DOSE

        output = output_dose+"+"+output_pheno

        if not os.path.isdir("../pheno+dose/"):
            os.mkdir("../pheno+dose/")

        pheno = pd.read_csv(PHENO, sep=" ")
        dose = pd.read_csv(DOSE, sep=" ")

        #Filter the SNPs from dose according with MAF and Rsq paramters
        if MAF:
            dose = dose[dose.MAF >= MAF]
        if Rsq:
            dose = dose[dose.Rsq >= Rsq]

        #Drop columns, uniqueSNPs should have just dosage information

        dose.drop(
            ["SNP", "POS", "REF.0.",
             "ALT.1.", "ALT_Frq", "MAF",
             "Rsq", "SNP.1", "A1", "A2"],
            axis=1, inplace=True)


        #Transpose matrix
        dose.set_index("uniqueSNP", inplace=True)
        dose = dose.transpose()


        pheno = pheno.drop(["FID"], axis=1)
        pheno = pheno.set_index("IID")

        pheno_dose = pd.concat([pheno,dose], join="inner", axis=1)

        pheno_dose.to_csv("../pheno+dose/"+output, sep="\t")

if __name__ == "__name__":

    #set up the command line parameters
    PARSER = argparse.ArgumentParser(
        description="This script builds a matrix containin \
                     information about phenotype and dosages \
                     for all samples and SNPs")

    PARSER.add_argument(
        "--MAF", help="MAF Frecuency to filter SNPs.\
        The number introduced will be the minimal SNPs\
        frecuency allowed in the files, including that number",
        type=float)

    PARSER.add_argument(
        "--Rsq", help="Imputation quality to filter SNPs.\
        The number introduced will be the minimal SNPs\
        imputation quality allowed in the files, including that number",
        type=float)

    PARSER.add_argument(
        "--pheno", help="Location of phenotype file",
        type=str)
    PARSER.add_argument(
        "--dose", help="Location of phenotype file",
        type=str)


    #Get matching parameters from the command line
    ARGS = PARSER.parse_args()

    MAF = ARGS.MAF
    Rsq = ARGS.Rsq
    PHENO = ARGS.pheno
    DOSE = ARGS.dose

    create_pheno_dosage_hndl = create_pheno_dosage()
    create_pheno_dosage_hndl.preprocess_matrix(MAF, Rsq, PHENO, DOSE)