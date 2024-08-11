

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# library("devtools");
# install_github("anRichment", username = "plangfelder");


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# library("githubinstall");
# githubinstall("anRichment");


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


# install.packages("BiocManager");
# BiocManager::install(c("AnnotationDBI", "GO.db", "org.Hs.eg.db", "org.Mm.eg.db", "XML", "WGCNA", 
#            "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Mmusculus.UCSC.mm10.knownGene"));


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# BiocManager::install(c("org.Rn.eg.db"));


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# install.packages("path/to/anRichment", repos = NULL, type = "source")


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# setwd("path/to/analysis");


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


options(stringsAsFactors = FALSE);
library("anRichment");
library("Hmisc");
library("WGCNA");


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Read in the module assignment data for Mike Oldham's modules
data = read.csv(file = bzfile("Data/MO_FxOrg_ColorVectors.csv.bz2"), sep = ",", header = TRUE)
# We will only keep the CTX modules
symbol.0 = data$CTX_Gene;
moduleColor = data$CTX_Module;
table(moduleColor)
# Some gene symbols have the form "XYZ /// ABC". Keep only the first symbol of all such multi-symbols.
split = strsplit(symbol.0, split = " /// ", fixed = TRUE);
symbol = sapply(split, function(x) x[1]);
# Convert symbols to Entrez IDs
entrez = convert2entrez(organism = "human", symbol = symbol);
# How many conversions were successful?
table(is.na(entrez))
# Build a Entrez ID-symbol translation table
entrez2symbol = data.frame(Entrez = entrez, Symbol = symbol);


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# GOcollection = buildGOcollection(organism = "human")


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


GOcollection = buildGOcollection(organism = "human", verbose = 0)


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


GOenrichment = enrichmentAnalysis.Entrez(
     classLabels = moduleColor, identifiers = entrez,
     refCollection = GOcollection,
     useBackground = "given",
     threshold = 1e-4,
     thresholdType = "Bonferroni",
     getOverlapEntrez = TRUE,
     getOverlapSymbols = TRUE,
     ID2symbol = entrez2symbol,
     ignoreLabels = "grey");

gc();


#=====================================================================================
#
#  Code chunk 12
#
#=====================================================================================


names(GOenrichment)


#=====================================================================================
#
#  Code chunk 13
#
#=====================================================================================


names(GOenrichment$enrichmentTable);


#=====================================================================================
#
#  Code chunk 14
#
#=====================================================================================


table.display = GOenrichment$enrichmentTable;
table.display$overlapGenes = shortenStrings(table.display$overlapGenes, maxLength = 70, 
                                            split = "|");
head(table.display);


#=====================================================================================
#
#  Code chunk 15
#
#=====================================================================================


write.csv(GOenrichment$enrichmentTable, file = "Results/GOenrichment-enrichmentTable.csv",
          row.names = FALSE);


#=====================================================================================
#
#  Code chunk 16
#
#=====================================================================================


names(GOenrichment$dataSetDetails)


#=====================================================================================
#
#  Code chunk 17
#
#=====================================================================================


names(GOenrichment$dataSetDetails[[1]][[1]])


#=====================================================================================
#
#  Code chunk 18
#
#=====================================================================================


GOenrichment$dataSetDetails$black[[3]]$commonGeneIDs


#=====================================================================================
#
#  Code chunk 19
#
#=====================================================================================


# GOenrichment$dataSetDetails[[1]][[3]]$commonGeneIDs


#=====================================================================================
#
#  Code chunk 20
#
#=====================================================================================


knownGroups(GOcollection)


#=====================================================================================
#
#  Code chunk 21
#
#=====================================================================================


GO.BPcollection = subsetCollection(GOcollection, tags = "GO.BP")


#=====================================================================================
#
#  Code chunk 22
#
#=====================================================================================


internalColl = internalCollection(organism = "human");


#=====================================================================================
#
#  Code chunk 23
#
#=====================================================================================


knownGroups(internalColl, sortBy = "size")


#=====================================================================================
#
#  Code chunk 24
#
#=====================================================================================


dataSetNames(internalColl, groups = "BrainLists.MO")
dataSetNames(internalColl, groups = "Nonexistent group")


#=====================================================================================
#
#  Code chunk 25
#
#=====================================================================================


ids = dataSetIDs(internalColl, groups = "BrainLists.MO")
dataSetNames(internalColl, dataSets = ids)


#=====================================================================================
#
#  Code chunk 26
#
#=====================================================================================


biosysCollection = BioSystemsCollection("human")


#=====================================================================================
#
#  Code chunk 27
#
#=====================================================================================


genomicPosCollection = genomicPositionCollection(
  organism = "human", 
  spacings = 5e6, 
  overlapFactor = 2)


#=====================================================================================
#
#  Code chunk 28
#
#=====================================================================================


HDSigCollection = HDSigDBCollection(organism = "human")


#=====================================================================================
#
#  Code chunk 29
#
#=====================================================================================


# msdbColl = MSigDBCollection(file = "path/to/msigdb.xlm", organism = "human")


#=====================================================================================
#
#  Code chunk 30
#
#=====================================================================================


# phenopediaColl = PhenopediaCollection(organism = "human")


#=====================================================================================
#
#  Code chunk 31
#
#=====================================================================================


WGCNA.HD.coll = HuntingtonsDiseaseWGCNACollection("human")


#=====================================================================================
#
#  Code chunk 32
#
#=====================================================================================


knownGroups(WGCNA.HD.coll, sortBy = "size")


#=====================================================================================
#
#  Code chunk 33
#
#=====================================================================================


scbtCollection = SCSBrainCellTypeCollection("human");
bdCollection = BrainDiseaseCollection("human");


#=====================================================================================
#
#  Code chunk 34
#
#=====================================================================================


combinedCollection = mergeCollections(
   GOcollection, 
   internalColl, 
   biosysCollection, 
   HDSigCollection);
knownGroups(combinedCollection, sortBy = "size")
@ 

One can then calculate the enrichment using the combined collection:

<<>>=
combinedEnrichment =  enrichmentAnalysis(classLabels = moduleColor, identifiers = entrez,
                                   refCollection = combinedCollection,
                                   useBackground = "given",
                                   threshold = 1e-4,
                                   thresholdType = "Bonferroni");


#=====================================================================================
#
#  Code chunk 35
#
#=====================================================================================


head(combinedEnrichment$enrichmentTable[, -ncol(combinedEnrichment$enrichmentTable)])


#=====================================================================================
#
#  Code chunk 36
#
#=====================================================================================


eLabels = enrichmentLabels(
   combinedEnrichment$enrichmentTable,
   focusOnGroups = c("all", "GO", "Cell type markers", "Brain region markers", "HDSigDB"),
   groupShortNames = c("all", "GO", "CT", "BR", "HD"),
   minSize = 0.05,
   numericClassLabels = FALSE);


#=====================================================================================
#
#  Code chunk 37
#
#=====================================================================================


bloodAtlasEnrichment =  enrichmentAnalysis(classLabels = moduleColor, identifiers = entrez,
                                   refCollection = combinedCollection,
                                   useGroups = c("BloodAtlases", "ImmunePathways"),
                                   useBackground = "given",
                                   threshold = 5e-2,
                                   nBestDataSets = 3,
                                   thresholdType = "Bonferroni");
head(bloodAtlasEnrichment$enrichmentTable[, -16])


#=====================================================================================
#
#  Code chunk 38
#
#=====================================================================================


active = entrez[moduleColor=="blue"];
all = entrez;
GOenrichment.blue = enrichmentAnalysis(active = active, inactive = all,
                                     refCollection = GOcollection,
                                     useBackground = "intersection",
                                     threshold = 1e-4,
                                     thresholdType = "Bonferroni");
head(GOenrichment.blue$enrichmentTable[, -16])


#=====================================================================================
#
#  Code chunk 39
#
#=====================================================================================


brainListColl = subsetCollection(internalColl, tags = "BrainLists");
nDataSets(brainListColl)


#=====================================================================================
#
#  Code chunk 40
#
#=====================================================================================


noBrainListColl = subsetCollection(internalColl, tags = "BrainLists", invertSearch = TRUE);
nDataSets(noBrainListColl)
@ 
We now have 166 gene sets left.

\subsection{Adding user-defined gene sets and collections programmatically}

An important capability of \anRichment\ is the ability for the user to create custom gene sets, groups,
and collections. We illustrate this procedure on a simple example of genes that are either in the blue or
black module. Our approach here is a bit naive since in a real analysis one should be more careful about
the probe to gene mapping. We start by selecting the blue and black genes and dropping missing Entrez
identifiers, 

<<>>=
moduleColorX = moduleColor;
moduleColorX[is.na(moduleColor)] = "grey";
bbGeneEntrez.0 = entrez[ moduleColorX %in% c("blue", "black") ];
# Some of the entrez codes are missing; we will drop them
bbGeneEntrez.1 = bbGeneEntrez.0[ !is.na(bbGeneEntrez.0) ];
# Multiple entrez codes are represented by several probes: keep only one copy of each.
bbGeneEntrez = unique(bbGeneEntrez.1);


#=====================================================================================
#
#  Code chunk 41
#
#=====================================================================================


knownEvidenceCodes()[, c(1:3)]


#=====================================================================================
#
#  Code chunk 42
#
#=====================================================================================


bbGeneSet = newGeneSet(
    geneEntrez = bbGeneEntrez, 
    geneEvidence = "IEP", 
    geneSource = paste0("Oldham MC et al, ", 
                        "Functional organization of the transcriptome in human brain",
                        "Nature Neuroscience 11, 1271 - 1282 (2008)"),
    ID = "dummy000001",
    name = "bb_M00_CTX", 
    description = "Blue or black genes from CTX network", 
    source = paste0("Oldham MC et al, ", 
                    "Functional organization of the transcriptome in human brain",
                    "Nature Neuroscience 11, 1271 - 1282 (2008)"), 
    organism = "human", 
    internalClassification = c("PL", "dummy"), 
    groups = "PL",
    lastModified = "2011-11-01");
@ 

In addition to the gene information, a gene set also contains several pieces of meta-information: a
unique identifier, a name (should also be unique but need not be), a short description, source (article
reference etc), organism for which the gene set is defined, internal classification (a vector of keywords
that are in principle arbitrary but should hopefully be organized in a hierarchical structure, from most
general to most specific), names of groups the gene set belongs to, and date of last modification. 
The meta-information helps the user
identify the meaning and source of gene sets and it is very important that as much information as possible
be included. 

We next create a group \code{"PL"} that is referenced in the gene set we just created.
<<>>=
PLgroup =  newGroup(name = "PL", description = "PL's experimental group of gene sets",
                                 source = "Personal imagination");


#=====================================================================================
#
#  Code chunk 43
#
#=====================================================================================


PLcollection = newCollection(dataSets = list(bbGeneSet), groups = list(PLgroup));


#=====================================================================================
#
#  Code chunk 44
#
#=====================================================================================


PLcollection = newCollection()
PLcollection = addToCollection(PLcollection, bbGeneSet, PLgroup)


#=====================================================================================
#
#  Code chunk 45
#
#=====================================================================================


PLenrichment = enrichmentAnalysis.Entrez(
    classLabels = moduleColor, identifiers = entrez,
    refCollection = PLcollection,
    useBackground = "given",
    threshold = 5e-2,
    ID2symbol = entrez2symbol,
    nBestDataSets = 3,
    thresholdType = "Bonferroni");

head(PLenrichment$enrichmentTable[, -16])


#=====================================================================================
#
#  Code chunk 46
#
#=====================================================================================


cellTypeColl = subsetCollection(internalColl, tags = "Cell type markers");
cellTypeDF = collection2dataFrames(cellTypeColl);


#=====================================================================================
#
#  Code chunk 47
#
#=====================================================================================


names(cellTypeDF)


#=====================================================================================
#
#  Code chunk 48
#
#=====================================================================================


head(shortenStrings(cellTypeDF$geneSetInfo))


#=====================================================================================
#
#  Code chunk 49
#
#=====================================================================================


cellTypeColl.2 = collectionFromDataFrames(
  geneSetInfoDF = cellTypeDF$geneSetInfo,
  geneSetContentDF = cellTypeDF$geneSetContent,
  groupDF = cellTypeDF$groupInfo);


#=====================================================================================
#
#  Code chunk 50
#
#=====================================================================================


groupA = newGroup(name = "A", description = "A", parents = "B");
groupB = newGroup(name = "B", description = "B", parents = "C");
groupC = newGroup(name = "C", description = "C");

groupLst = list(groupA, groupB, groupC);


#=====================================================================================
#
#  Code chunk 51
#
#=====================================================================================


impliedGroups(groupLst, get = "parents");


#=====================================================================================
#
#  Code chunk 52
#
#=====================================================================================


impliedGroups(groupLst, get = "children");


#=====================================================================================
#
#  Code chunk 53
#
#=====================================================================================


organismLabels()


#=====================================================================================
#
#  Code chunk 54
#
#=====================================================================================


names(bbGeneSet)


#=====================================================================================
#
#  Code chunk 55
#
#=====================================================================================


class(bbGeneSet)


#=====================================================================================
#
#  Code chunk 56
#
#=====================================================================================


names(PLgroup)


#=====================================================================================
#
#  Code chunk 57
#
#=====================================================================================


names(internalColl)


#=====================================================================================
#
#  Code chunk 58
#
#=====================================================================================


sessionInfo()


