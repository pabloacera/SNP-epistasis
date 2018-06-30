"""
Testing the logistic regression analysis and likelihood ratio test
try the log and lkt in the titanic dataset
"""
from __future__ import print_function #pylint: disable=syntax-error
from pydataset import data
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np

from scipy import stats


def LRT(results_res, results):
	llf_full = results.llf
	llf_restr = results_res.llf
	df_full = results.df_resid
	df_restr = results_res.df_resid
	lrdf = (df_restr - df_full)
	lrstat = -2*(llf_restr - llf_full)
	lr_pvalue = stats.chi2.sf(lrstat, df=lrdf)
	return lr_pvalue

bank = pd.read_csv("banking.csv")

print(bank.head())


#TITANIC = data("titanic")

#[class, age,sex, survived]


aditive_model = "y ~ age + nr_employed"
interaction_model = "y ~ age + age * nr_employed"


logitfit_adi = smf.logit(formula=str(aditive_model), data=bank).fit()
logitfit_inte = smf.logit(formula= str(interaction_model), data=bank).fit()
print(logitfit_adi.summary())
print(logitfit_inte.summary())

print(LRT(logitfit_adi, logitfit_inte))
#print(LRT(logitfit_adi.llr, logitfit_inte.llr, logitfit_inte))
#Other formula
#import statsmodels.discrete.discrete_model

#discrete_model.Logit(bank.y, np.array())