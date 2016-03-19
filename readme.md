Work in progress. Instruction is incomplete. 

# Egene 
Egene is a gene that has at least one genetic variant associated with its expression (measured by amount of RNA produced). There is also Pgene which is the analogous, where the gene expression is measured by the amount of protein produced. 

# Finding Egenes
This executable finds genes whose expressions are associated with at least one genetic variant. The code uses simple permutation test where the expression measurements are swapped among the individuals, while the genetic matrix stays unchanged so that linkage disequilibrium is kept throughout. 

To run, use the following 

./computeOneGenePval_permutation 
            --var [covariates age, gender etc]
            --numCov [# components such as age, gender etc] 
            --expr [gene expression] 
            --geno [variants] 
            --numSNP [# variants] 
            --w [prior]
            --final [permutation test p-value output] 
            --pval [minimum p-value among all variants output] 
            --lrt [maximum likelihood ratio output] 
            --permutation [# permutations] 
            --seed [random seed] 
            --alt [alternative mean] 
            --reuse [reuse estimated mean as alternative mean] 
            
