library(bio3d)

# Attempt to automate the read.pdb function
easypdb.plot <- function(x, b = 1) {
  
  y <- read.pdb(x)
  
  A.chainA <- trim.pdb(y, chain="A", elety="CA")

  Z <- A.chainA$atom$b
   
  plotb3(Z, sse=A.chainA, typ="l", ylab="B Factor", col = "coral", lwd=2, main = x)
  
}

easypdb.plot("4AKE")
 



  
  
