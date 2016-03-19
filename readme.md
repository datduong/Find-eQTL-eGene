# Egene 
Egene is a gene that has at least one genetic variant associated with its expression (measured by amount of RNA produced). There is also Pgene which is the analogous, where the gene expression is measured by the amount of protein produced. 

# Finding Egenes
This executable finds genes whose expressions are associated with at least one genetic variant. The code uses simple permutation test where the expression measurements are swapped among the individuals, while the genetic matrix stays unchanged so that linkage disequilibrium is kept throughout. 

To run, use the following 

./computeOneGenePval_permutation --expr [gene expression] --geno [variants] --final [permutation test p-value] --permutate [# permutations] --pval [minimum p-value among all SNPs] --numSNP [# SNPs] --seed [random seed] -- numCov [# components such as age, gender etc] --w [prior] --alt [alternative mean] --gradient [maximum likelihood ratio] 
            
