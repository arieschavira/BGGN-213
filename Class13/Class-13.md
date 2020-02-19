Class13: Drug Discovery
================
Aries Chavira
2/19/2020

Retrieve and process starting HIV protease structure
----------------------------------------------------

Here we download and clean the HIV-Pr structure (pdb codeL 1HSG) from PDB database. We will separate the "protein-only" and "ligand only" into two PDB files

``` r
library(bio3d)
```

``` r
file.name <- get.pdb("1HSG") 
```

    ## Warning in get.pdb("1HSG"): ./1HSG.pdb exists. Skipping download

We will use `read.pdb`, `atom.select`, and `write.pdb` functions to make our seperate "Protein-only" and "Ligand only" PDB files

``` r
HIV <- read.pdb("1HSG")
```

    ##   Note: Accessing on-line PDB file

``` r
HIV
```

    ## 
    ##  Call:  read.pdb(file = "1HSG")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
Protein <- atom.select(HIV, "protein", value = TRUE)

Ligand <- atom.select(HIV, "ligand", value = TRUE)
```

``` r
prot.pdb <- write.pdb(Protein, file = "1hsg_protein.pdb")

lig.pdb <- write.pdb(Ligand, file = "1hsg_ligand.pdb")
```

Reading docking results
-----------------------

converting results into a .pdb instead of .pdbqt

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```
