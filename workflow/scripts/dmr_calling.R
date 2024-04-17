library(methylKit)
# file.list=list( system.file("extdata", 
#                             "test1.myCpG.txt", package = "methylKit"),
#                 system.file("extdata",
#                             "test2.myCpG.txt", package = "methylKit"),
#                 system.file("extdata", 
#                             "control1.myCpG.txt", package = "methylKit"),
#                 system.file("extdata", 
#                             "control2.myCpG.txt", package = "methylKit") )

file.list <- c(
    snakemake@input[['_1']],
    snakemake@input[['_2']],
    snakemake@input[['_3']],
    snakemake@input[['_4']],
)



# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("test1","test2","ctrl1","ctrl2"),
           assembly="hg18",
           treatment=c(1,1,0,0),
           context="CpG",
           mincov = 10
           )