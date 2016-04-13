Work in progress.

Please begin at this parent website http://genetics.cs.ucla.edu/egene-mvn/index.html. 

# Egene 
Egene is a gene that has at least one genetic variant associated with its expression (measured by amount of RNA produced). There is also Pgene which is the analogous, where the gene expression is measured by the amount of protein produced. 

# Finding Egenes
This executable finds genes whose expressions are associated with at least one genetic variant. The code uses simple permutation test where the expression measurements are swapped among the individuals, while the genetic matrix stays unchanged so that linkage disequilibrium is kept throughout. 

A second faster method relies on the MVN distribution. One would generate the null density by sampling from the MVN (with mean 0, and covariance as the SNP-covariance matrix). It is a bit too involved to fully explain the detail here. Please refer to the parent website at UCLA Genetics Department. http://genetics.cs.ucla.edu/egene-mvn/index.html 

To run, use the executable, and read the instruction.txt. 


            
