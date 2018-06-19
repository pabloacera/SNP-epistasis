##### Command-line arguments: chromosome, loop
##### The script runs as follows: Rscript tadtagger.R chromosome loop_number
chr=commandArgs(trailingOnly=T)[1]
loop=commandArgs(trailingOnly=T)[2]

##### Read in the genotype dosages
dose <- read.table(paste("dosage/WTCCC.chr",chr,".hrc.imputed.plink.dose.info.0.6.maf.0.01.txt.gz",sep=""),header=T)

##### Read in the covariates and phenotype
pheno <- read.table("./WTCCC.hrc.imputed.plink.dose.info.0.9.hardcalls.pca.txt",header=T)

##### Create a matrix to store all of the results
##### The matrix will print the loop number, the chromosome/position, the SNPs tested (snp A, snp B), the GWAS p-value of the SNPs
##### and 4 LRT tests: the interaction of snpA and snpB vs (1) the null model, (2) snp A only, (3) snp B only, (4) snp A + snp B
result <- matrix(ncol=13,nrow=0)
colnames(result) <- c("loop","chrA","posA","snpA","chrB","posB","snpB","snpA.gwas","snpB.gwas","lrt.null","lrt.A","lrt.B","lrt.AB")

##### Read in the looping SNPs
##### "A" and "B" are the anchor regions of the topologically associating domain (TAD)
loop.name <- paste("loop",loop,sep="")
loopA <- read.table(paste("./ld/chr",chr,"/WTCCC.chr",chr,".hrc.imputed.plink.dose.info.0.6.hardcalls.loop",loop,".A.prune.in",sep=""),header=F)
loopB <- read.table(paste("./ld/chr",chr,"/WTCCC.chr",chr,".hrc.imputed.plink.dose.info.0.6.hardcalls.loop",loop,".B.prune.in",sep=""),header=F)

##### Just to be sure, make sure that the looping SNPs overlap 100% with the SNPs in the dataframe
##### If not, drop the SNPs that don't overlap the dosage dataframe
snps <- as.matrix(dose$uniqueSNP)
loopA <- merge(loopA,snps)
loopB <- merge(loopB,snps)

##### Use each SNP in loopA (i.e., anchor region A) as the 'index' and then tick through all possible SNPs in loopB (i.e., anchor region B)

for(j in 1:dim(loopA)[1]) {
	
		snpA <- loopA[j,1]
		
		### Temporarily bind together the dosage for the SNP
		tmp.doseA <- as.matrix(t(subset(dose,dose$uniqueSNP==as.character(snpA))[12:dim(dose)[2]]))
		tmp.idA <- rownames(tmp.doseA)
		rownames(tmp.doseA) <- NULL
		tmp.doseA <- cbind(tmp.idA,tmp.doseA); colnames(tmp.doseA) <- c("IID","snpA")
		
		### Merge together the dosages with the phenotype
	 	tmp.data <- merge(pheno,tmp.doseA,by="IID")
 	
 	### Now for every SNP in loopB, run a series of tests
 	for(k in 1:dim(loopB)[1]) {
 		
 		snpB <- loopB[k,1]
 		
 		### Add on these dosages (for snp B)
 		tmp.doseB <- as.matrix(t(subset(dose,dose$uniqueSNP==as.character(snpB))[12:dim(dose)[2]]))
		tmp.idB <- rownames(tmp.doseB)
		rownames(tmp.doseB) <- NULL
		tmp.doseB <- cbind(tmp.idB,tmp.doseB); colnames(tmp.doseB) <- c("IID","snpB")
	 	tmp.data <- merge(tmp.data,tmp.doseB,by="IID")
 	
 		### Now test snpA and snpB individually in the logistic framework (i.e., this is the GWAS)
		### GWAS result, snp A
 		modelA <- formula(paste("as.factor(PHENO) ~ as.numeric(snpA) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10",sep=""))
 		logit.A <- glm(modelA,data=tmp.data,family="binomial")
 		
		### GWAS result, snp B
 		modelB <- formula(paste("as.factor(PHENO) ~ as.numeric(snpB) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10",sep=""))
 		logit.B <- glm(modelB,data=tmp.data,family="binomial")
 		
 		### Logistic for snpA (extract the p-value)
 		snpA.logitp <- summary(logit.A)$coeff[2,4]
 		
 		### Logistic for snpB (extract the p-value)
 		snpB.logitp <- summary(logit.B)$coeff[2,4]
 	
 		### Now, move on to the LRT tests
		### Should we only do this if snp A or snp B p < 0.05? Or all snps?
 		### Compare the interaction test to the null model (PCs only) as well as to the models that contain only snpA or only snpB
 
 		### The null model, covariates only
  		null <- formula(paste("as.factor(PHENO) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10",sep=""))
 		logit.null <- glm(null, data=tmp.data, family="binomial")

 		### The additive model
 		modelAB <- formula(paste("as.factor(PHENO) ~ as.numeric(snpA) + as.numeric(snpB) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10",sep=""))
 		logit.AB <- glm(modelAB, data=tmp.data, family="binomial")

 		### The interaction model
 		modelABinter <- formula(paste("as.factor(PHENO) ~ as.numeric(snpA) + as.numeric(snpB) + as.numeric(snpA)*as.numeric(snpB) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10",sep=""))
 		logit.ABinter <- glm(modelABinter, data=tmp.data, family="binomial")

		### Now run all of the LRT tests that compare the full (interaction) model to all other models
		### Extract the p-value from each LRT test
 		lrt.null <- anova(logit.null, logit.ABinter, test="Chisq")$P[2]
 		lrt.A <- anova(logit.A, logit.ABinter, test="Chisq")$P[2]
 		lrt.B <- anova(logit.B, logit.ABinter, test="Chisq")$P[2]
 		lrt.AB <- anova(logit.AB, logit.ABinter, test="Chisq")$P[2]
 		
 		### Now assemble all the information on the SNP pair, for plotting/further analysis
 		chr.snpA <- strsplit(as.character(snpA),":")[[1]][1]; pos.snpA <- strsplit(as.character(snpA),":")[[1]][2]
 		chr.snpB <- strsplit(as.character(snpB),":")[[1]][1]; pos.snpB <- strsplit(as.character(snpB),":")[[1]][2]
 		
 		result <- rbind(result, c(as.character(loop.name),chr.snpA,pos.snpA,as.character(snpA),chr.snpB,pos.snpB,as.character(snpB),snpA.logitp,snpB.logitp,lrt.null,lrt.A,lrt.B,lrt.AB))
		
		### Slice the second SNP off of tmp.data so we can start again
		tmp.data <- tmp.data[,1:14]
 	
 	}
	
}

### When the analysis is all done, write out the results to this file
write.table(result,file=paste("./result/chr",chr,"/WTCCC.chr",chr,".hrc.imputed.plink.dose.info.0.6.maf.0.01.loop",loop,".lrt.results.txt",sep=""),quote=F,row.names=F,col.names=T)
