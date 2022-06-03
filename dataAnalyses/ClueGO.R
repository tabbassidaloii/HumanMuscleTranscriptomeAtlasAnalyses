
# required libraries that need to be installed before starting the script
require (xml2)
require (RJSONIO)
require (httr) # needs the libssl-dev and libcurl4-openssl-dev packages in linux (ubuntu: sudo apt-get install libssl-dev libcurl4-openssl-dev)

rm (list = ls ())

# Helper function to transform output into a data frame
text.to.data.frame <- function (table.text) {
  table <- NULL
  rows <- unlist (strsplit (result.table.text, split = "\n"))
  header <- t (unlist (strsplit (rows[1], split = "\t")))
  for (i in 2:length (rows)) {
    if (is.null (table)) {
      table <- t (unlist (strsplit (rows[i], split = "\t")))
    } else {
      table <- rbind (table, t (unlist (strsplit (rows[i], split = "\t"))))
    }
  }
  table <- as.data.frame (table)
  names (table) <- header
  return (table)
}


#### choose the example to run with one or two gene lists ####
Enrichment.Analysis = "ClueGO"

#### Basic settings for cyREST ####
home.folder = "OurDataset/MuscleDataset/Outputs"
output.folder = paste (home.folder, Enrichment.Analysis, sep = "/")
dir.create (output.folder, showWarnings = F)

cluego.home.folder = paste ("C://Users/", "ClueGOConfiguration", "v2.5.7", sep = "/")
port.number = 1234
host.address <- "localhost"

# define base urls
cytoscape.base.url = paste ("http://", host.address, ":", toString (port.number), "/v1", sep = "")
cluego.base.url = paste (cytoscape.base.url, "apps", "cluego", "cluego-manager", sep = "/")

print (paste ("Gene Modules ", Enrichment.Analysis, sep = ""))
print (paste ("User Home Folder: ", home.folder, sep = ""))
print (paste ("User Output Folder: ", output.folder, sep = ""))
print (paste ("Cytoscape Base URL: ", cytoscape.base.url, sep = ""))
print (paste ("ClueGO Base URL: ", cluego.base.url, sep = ""))

# Modules <- setdiff (paste0 ("M.", 1:34), paste0 ("M.", c (2, 20, 28, 22, 19, 9, 26)))
# Modules <- paste0 ("M.", c (14, 16, 23, 4, 33, 31, 18, 13))
# Modules <- paste0 ("M.", c (24, 17, 7, 3, 32, 30, 15, 12, 6))
# Modules <- paste0 ("M.", c (10, 11, 25, 5, 8, 27))
Modules <- paste0 ("M.", c (34, 29, 21, 1))

for (M in Modules) {
response <- POST ("http://localhost:1234/v1/apps/cluego/start-up-cluego", encode = "json")
  
#### 0.0 Start up ClueGO in case it is not running yet
response <- POST (url = paste (cytoscape.base.url, "apps", "cluego", "start-up-cluego", sep = "/"),
                  encode = "json")
# wait 2 seconds to make sure ClueGO is started
if (grepl ("HTTP ERROR 500", response)) {
  print ("wait 2 secs")
  Sys.sleep (2)
}

#### 1.0 Select the ClueGO Organism to analyze ####
organism.name = "Homo Sapiens" # (run "1.1 Get all ClueGO organisms" to get all options)
print (paste ("1.0 Select the ClueGO Organism to analyze: ", organism.name, sep = ""))
response <- PUT (url = paste (cluego.base.url, "organisms", "set-organism", URLencode (organism.name), 
                              sep = "/"), encode = "json")

#### 2.0 Upload IDs for a specific module ####
print (paste ("2.0 Upload IDs for a specific module", sep = ""))

## [optional functions and settings, un comment and modify if needed]
# 2.1 Select the ClueGO ID type
id.type.name = "EnsemblGeneID" # (run "2.3 Get ClueGO ID types" to get all options)
response <- PUT (url = paste (cluego.base.url, "ids", "set-id-type", id.type.name, sep = "/"), encode = "json")

# 2.4 Set the number of modules
max.input.panel.number = 1
response <- PUT (url = paste (cluego.base.url, "cluster", "max-input-panel", max.input.panel.number, sep = "/"), encode = "json")


## [optional functions and settings, un comment and modify if needed]
# 2.5 Set analysis properties for a module, this is to repeat for each input module
cluster1 = "1"

input.panel.index = cluster1 # set here the module number e.g. "1"
node.shape = "Ellipse" # ("Ellipse", "Diamond", "Hexagon", "Octagon", "Parallelogram", "Rectangle", "Round Rectangle", "Triangle", "V")
cluster.color = "#ff0000" # The color in hex, e.g. #F3A455
min.number.of.genes.per.term = 4
min.percentage.of.genes.mapped = 4
no.restrictions = FALSE # TRUE for no restricions in number and percentage per term
response <- PUT (url = paste (cluego.base.url, "cluster", "set-analysis-properties", input.panel.index, node.shape, URLencode (cluster.color, reserved = TRUE), min.number.of.genes.per.term, min.percentage.of.genes.mapped, no.restrictions, sep = "/"), encode = "json")

# set gene list for cluster 1
file.location = paste0 (home.folder, "/", "ClueGO/GenesModules/", M, ".txt")
gene.list1 <- toJSON (read.table (file.location, as.is = TRUE)[[1]])
response <- PUT (url = paste (cluego.base.url, "cluster", "upload-ids-list", URLencode (cluster1), sep = "/"), body = gene.list1, encode = "json", content_type_json ())

# 2.6 Select visual style
visual.style = "ShowGroupDifference"  # (ShowGroupDifference, ShowSignificanceDifference, ShowClusterDifference)
response <- PUT (url = paste (cluego.base.url, "cluster", "select-visual-style", visual.style, sep = "/"), encode = "json")

####  3.0 Select Ontologies
print (paste ("3.0 Select Ontologies", sep = ""))

selected.ontologies <- toJSON (c ("3;Ellipse", "4;Hexagon", "6;Octagon", 
                                  "8;Triangle", "9;Rectangle", "11;Round Rectangle")) # (run "3.1 Get all available Ontologies" to get all options)
response <- PUT (url = paste (cluego.base.url, "ontologies", "set-ontologies", sep = "/"), body = selected.ontologies, encode = "json", content_type_json ())

## [optional functions and settings, un comment and modify if needed]
# 3.1 Get all available Ontologies
# response <- GET (paste (cluego.base.url, "ontologies", "get-ontology-info", sep = "/"))
# print (content (response))
# 
# # 3.2 Select Evidence Codes
# evidence.codes <- toJSON (c ("All")) # (run "3.3 Get all available Evidence Codes" to get all options)
# response <- PUT (url = paste (cluego.base.url, "ontologies", "set-evidence-codes", sep = "/"), body = evidence.codes, encode = "json", content_type_json ())
# 
# # 3.3 Get all available Evidence Codes
# response <- GET (paste (cluego.base.url, "ontologies", "get-evidence-code-info", sep = "/"))
# print (content (response))
# 
# 3.4 Set min, max GO tree level
min.level = 5
max.level = 6
all.levels = FALSE
response <- PUT (url = paste (cluego.base.url, "ontologies", "set-min-max-levels", min.level, max.level, all.levels, sep = "/"), encode = "json")

# # 3.5 Use GO term significance cutoff
p.value.cutoff = 0.05
use.significance.cutoff = TRUE
response <- PUT (url = paste (cluego.base.url, "ontologies", use.significance.cutoff, p.value.cutoff, sep = "/"), encode = "json")
# 
# # 3.6 Use GO term fusion
# use.go.term.fusion = TRUE
# response <- PUT (url = paste (cluego.base.url, "ontologies", use.go.term.fusion, sep = "/"), encode = "json")

# # 3.7 Set statistical parameters
# ("Enrichment (Right-sided hypergeometric test)", "Depletion (Left-sided hypergeometric test)", "Enrichment/Depletion (Two-sided hypergeometric test)")
enrichment.type = "Enrichment (Right-sided hypergeometric test)" 
multiple.testing.correction = "Benjamini-Hochberg"
use.mid.pvalues = FALSE
use.doubling = FALSE
use.custom.reference.set = TRUE
custom.reference.set = paste (output.folder, "GenesModules/All.txt", sep = "/")

response <- PUT (url = paste (cluego.base.url, "stats", enrichment.type, multiple.testing.correction, 
                              use.mid.pvalues, use.doubling, use.custom.reference.set, custom.reference.set,
                              sep = "/"), encode = "json")
# 
# # 3.8 Set the Kappa Score level
# kappa.score = 0.4
# response <- PUT (url = paste (cluego.base.url, "ontologies", "set-kappa-score-level", kappa.score, sep = "/"), encode = "json")
# 
# # 3.9 Set grouping parameters
# do.grouping = TRUE
# coloring.type = "Random" # ("Random", "Fix")
# group.leading.term = "Highest Significance" # ("Highest Significance", "#Genes / Term", "%Genes / Term", "%Genes / Term vs Cluster")
# grouping.type = "Kappa Score" # ("Kappa Score", "Tree")
# init.group.size = 1
# perc.groups.for.merge = 50
# perc.terms.for.merge = 50
# response <- PUT (url = paste (cluego.base.url, "grouping", do.grouping, coloring.type, group.leading.term, grouping.type, init.group.size, perc.groups.for.merge, perc.terms.for.merge, sep = "/"), encode = "json")


#### 4.0 Run ClueGO Analysis ####
print (paste ("4.0 Run ClueGO Analysis", sep = ""))
# Run the analysis an save log file
analysis.name <- M
analysis.option <- "Cancel and refine selection" # ("Continue analysis", "Skip the grouping", "Cancel and refine selection")  -> Analysis option in case there are more than 1000 terms found!
response <- GET (paste (cluego.base.url, URLencode (analysis.name), URLencode (analysis.option), sep = "/"))
if (!grepl ("error", http_status ( response)["category"])) {
output.folder.M <- paste (output.folder, M, sep = "/")
dir.create (output.folder.M)
log.file.name = paste0 (output.folder.M, "/", "ClueGOlog", M, ".txt")
writeLines (content (response, encoding = "UTF-8"), log.file.name)
# print (content (response, encode = "text"))

# 4.1 Get network id (SUID) (CyRest function from Cytoscape)
response <- GET (paste (cytoscape.base.url, "networks", "currentNetwork", sep = "/"))
current.network.suid <- content (response, encode = "json")$data$networkSUID
# print (current.network.suid)

print (paste ("Save results", sep = ""))

# Get network graphics (CyRest function from Cytoscape)
image.type = "svg" # png, pdf
response <- GET (paste (cytoscape.base.url, "networks", current.network.suid, "views", paste ("first.", image.type, sep = ""), sep = "/"))
image.file.name = paste (output.folder.M, paste ("ClueGONetwork.", image.type, sep = ""), sep = "/")
writeBin (content (response, encode = "raw"), image.file.name)

# 4.2 Get ClueGO result table
response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-table", current.network.suid, sep = "/"))
result.table.text <- content (response, encode = "text", encoding = "UTF-8")
table.file.name = paste (output.folder.M, "ClueGO-Result-Table.txt", sep = "/")
write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")

# 4.3 Get ClueGO genes and main functions
# number.of.functions.to.add = 3
# response <- GET (paste (cluego.base.url, "analysis-results", "get-main-functions", current.network.suid, number.of.functions.to.add, sep = "/"))
# result.table.text <- content (response, encode = "text", encoding = "UTF-8")
# table.file.name = paste (output.folder.M, "ClueGO-Genes-With-Main-Functions.txt", sep = "/")
# write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")

# 4.4 Get genes result table
# response <- GET (paste (cluego.base.url, "analysis-results", "get-gene-table", current.network.suid, sep = "/"))
# result.table.text <- content (response, encode = "text", encoding = "UTF-8")
# table.file.name = paste (output.folder.M, "ClueGO-Gene-Table.txt", sep = "/")
# write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")

# 4.5 Get Kappascore Matrix
response <- GET (paste (cluego.base.url, "analysis-results", "get-kappascore-matrix", current.network.suid, sep = "/"))
result.table.text <- content (response, encode = "text", encoding = "UTF-8")
table.file.name = paste (output.folder.M, "ClueGO-Kappascore-Matrix.txt", sep = "/")
write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")

# 4.6 Get binary Gene-Term Matrix
response <- GET (paste (cluego.base.url, "analysis-results", "get-binary-gene-term-matrix", current.network.suid, sep = "/"))
result.table.text <- content (response, encode = "text", encoding = "UTF-8")
table.file.name = paste (output.folder.M, "ClueGO-Binary-Gene-Term-Matrix.txt", sep = "/")
write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")

# 4.7 ClueGO Result Chart
# Get result charts for both cluster as pie chart
chart.type = "BarChart" # ("PieChart", "BarChart")
image.type = "svg" # ("svg", "png", "pdf")
response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-result-chart", current.network.suid, cluster1, chart.type, image.type, sep = "/"))
image.file.name = paste (output.folder.M, paste ("ClueGO-", chart.type, "-For-Cluster", cluster1, ".", image.type, sep = ""), sep = "/")
writeBin (content (response, encode = "raw"), image.file.name)
if (Enrichment.Analysis ==  "ClueGO Rmultiple gene lists") {
  response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-result-chart", current.network.suid, cluster2, chart.type, image.type, sep = "/"))
  image.file.name = paste (output.folder.M, paste ("ClueGO-", chart.type, "-For-Cluster", cluster2, ".", image.type, sep = ""), sep = "/")
  writeBin (content (response, encode = "raw"), image.file.name)
}

## [optional functions and settings, un comment and modify if needed]
# # 4.8 Remove ClueGO analysis result
# print (paste ("Remove ClueGO Network", sep = ""))
# # # Remove analysis to reduce memory usage. This is important when using patch modes that create lots of analyses.
# response <- DELETE (paste (cluego.base.url, "remove-cluego-analysis-result", current.network.suid, sep = "/"))
# 
# # 4.9 Remove all ClueGO analysis results
# print (paste ("Remove all ClueGO Networks", sep = ""))
# response <- DELETE (paste (cluego.base.url, "remove-all-cluego-analysis-results", sep = "/"))
# 
# print (paste (Enrichment.Analysis, M, "done", sep = " "))

rm (list = setdiff(ls(), c ("cluego.base.url", "cluego.home.folder", "cytoscape.base.url",
                            "Enrichment.Analysis", "home.folder", "host.address", "output.folder",
                            "port.number", "text.to.data.frame", "Modules" )))
}
}

filenames <- list.files (path = "OurDataset/MuscleDataset/Outputs/ClueGO",
                         recursive = TRUE, pattern = "Result-Table.txt", full.names = T)

filenames <- filenames [grep ("\\/M\\.", filenames)]
ClueGO_enrichment <- c ()
for (file in filenames){
  ClueGO <- read.csv (file, sep = "\t")
  ClueGO$Module <- gsub (".*ClueGO/|/ClueGO.*", "", file)
  ClueGO_enrichment <- rbind (ClueGO_enrichment, ClueGO)
}

ClueGO_enrichment$GroupModule <- ifelse (ClueGO_enrichment$Module %in%
                                           paste0 ("M.", c (14, 16, 23, 4, 33, 31, 18, 13)), "HigerInGRVV",
                                         ifelse (ClueGO_enrichment$Module %in%
                                                   paste0 ("M.", c (24, 17, 7, 3, 32, 30, 15, 12, 6)), "HigerInSSG", 
                                                 ifelse (ClueGO_enrichment$Module %in% 
                                                           paste0 ("M.", c (10, 11, 25, 5, 8, 27)), "HigerInG", "Other")))
ClueGO_enrichment <- ClueGO_enrichment [, c (14, 15, 1:13)]

write.csv (ClueGO_enrichment, "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableInterstingModules.csv", row.names = F)
save (ClueGO_enrichment, file = "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableInterstingModules.RData")


# DC
dir.create (paste (output.folder, "DifferentialCoexpression", sep = "/"))
"HigherConnectivityInSEMthanGRA"
load ("Outputs/NetworkPerMuscleWGCNA/GeneConnectivity/DifferentiallyCoexpressedGenes.RData")
for (DiffCo in names (DiffCoExpGenes)){
  response <- POST ("http://localhost:1234/v1/apps/cluego/start-up-cluego", encode = "json")
  
  #### 0.0 Start up ClueGO in case it is not running yet
  response <- POST (url = paste (cytoscape.base.url, "apps", "cluego", "start-up-cluego", sep = "/"),
                    encode = "json")
  # wait 2 seconds to make sure ClueGO is started
  if (grepl ("HTTP ERROR 500", response)) {
    print ("wait 2 secs")
    Sys.sleep (2)
  }
  
  #### 1.0 Select the ClueGO Organism to analyze ####
  organism.name = "Homo Sapiens" # (run "1.1 Get all ClueGO organisms" to get all options)
  print (paste ("1.0 Select the ClueGO Organism to analyze: ", organism.name, sep = ""))
  response <- PUT (url = paste (cluego.base.url, "organisms", "set-organism", URLencode (organism.name), 
                                sep = "/"), encode = "json")
  
  #### 2.0 Upload IDs for a specific list ####
  print (paste ("2.0 Upload IDs for", DiffCo, sep = " "))
  
  ## [optional functions and settings, un comment and modify if needed]
  # 2.1 Select the ClueGO ID type
  id.type.name = "EnsemblGeneID" # (run "2.3 Get ClueGO ID types" to get all options)
  response <- PUT (url = paste (cluego.base.url, "ids", "set-id-type", id.type.name, sep = "/"), encode = "json")
  
  # 2.4 Set the number of modules
  max.input.panel.number = 1
  response <- PUT (url = paste (cluego.base.url, "cluster", "max-input-panel", max.input.panel.number, sep = "/"), encode = "json")
  
  
  ## [optional functions and settings, un comment and modify if needed]
  # 2.5 Set analysis properties for a module, this is to repeat for each input module
  cluster1 = "1"
  
  input.panel.index = cluster1 # set here the module number e.g. "1"
  node.shape = "Ellipse" # ("Ellipse", "Diamond", "Hexagon", "Octagon", "Parallelogram", "Rectangle", "Round Rectangle", "Triangle", "V")
  cluster.color = "#ff0000" # The color in hex, e.g. #F3A455
  min.number.of.genes.per.term = 4
  min.percentage.of.genes.mapped = 4
  no.restrictions = FALSE # TRUE for no restricions in number and percentage per term
  response <- PUT (url = paste (cluego.base.url, "cluster", "set-analysis-properties", input.panel.index, node.shape, URLencode (cluster.color, reserved = TRUE), min.number.of.genes.per.term, min.percentage.of.genes.mapped, no.restrictions, sep = "/"), encode = "json")
  
  # set gene list for cluster 1
  gene.list1 <- toJSON (DiffCoExpGenes [[DiffCo]])
  response <- PUT (url = paste (cluego.base.url, "cluster", "upload-ids-list", URLencode (cluster1), sep = "/"), body = gene.list1, encode = "json", content_type_json ())
  # 2.6 Select visual style
  visual.style = "ShowGroupDifference"  # (ShowGroupDifference, ShowSignificanceDifference, ShowClusterDifference)
  response <- PUT (url = paste (cluego.base.url, "cluster", "select-visual-style", visual.style, sep = "/"), encode = "json")
  ####  3.0 Select Ontologies
  print (paste ("3.0 Select Ontologies", sep = ""))
  selected.ontologies <- toJSON (c ("3;Ellipse", "4;Hexagon", "6;Octagon", 
                                    "8;Triangle", "9;Rectangle", "11;Round Rectangle")) # (run "3.1 Get all available Ontologies" to get all options)
  response <- PUT (url = paste (cluego.base.url, "ontologies", "set-ontologies", sep = "/"), body = selected.ontologies, encode = "json", content_type_json ())
  
  # 3.4 Set min, max GO tree level
  min.level = 5
  max.level = 6
  all.levels = FALSE
  response <- PUT (url = paste (cluego.base.url, "ontologies", "set-min-max-levels", min.level, max.level, all.levels, sep = "/"), encode = "json")
  # # 3.5 Use GO term significance cutoff
  p.value.cutoff = 0.05
  use.significance.cutoff = TRUE
  response <- PUT (url = paste (cluego.base.url, "ontologies", use.significance.cutoff, p.value.cutoff, sep = "/"), encode = "json")
   
  # # 3.7 Set statistical parameters
  # ("Enrichment (Right-sided hypergeometric test)", "Depletion (Left-sided hypergeometric test)", "Enrichment/Depletion (Two-sided hypergeometric test)")
  enrichment.type = "Enrichment (Right-sided hypergeometric test)" 
  multiple.testing.correction = "Benjamini-Hochberg"
  use.mid.pvalues = FALSE
  use.doubling = FALSE
  use.custom.reference.set = TRUE
  custom.reference.set = paste (output.folder, "GenesModules/All.txt", sep = "/")
  response <- PUT (url = paste (cluego.base.url, "stats", enrichment.type, multiple.testing.correction, 
                                use.mid.pvalues, use.doubling, use.custom.reference.set, custom.reference.set,
                                sep = "/"), encode = "json")
  #### 4.0 Run ClueGO Analysis ####
  print (paste ("4.0 Run ClueGO Analysis", sep = ""))
  # Run the analysis an save log file
  analysis.name <- DiffCo
  analysis.option <- "Cancel and refine selection" # ("Continue analysis", "Skip the grouping", "Cancel and refine selection")  -> Analysis option in case there are more than 1000 terms found!
  response <- GET (paste (cluego.base.url, URLencode (analysis.name), URLencode (analysis.option), sep = "/"))
  if (!grepl ("error", http_status ( response)["category"])) {
    output.folder.Diff <- paste (output.folder, "DifferentialCoexpression", DiffCo, sep = "/")
    dir.create (output.folder.Diff)
    log.file.name = paste0 (output.folder.Diff, "/", "ClueGOlog", DiffCo, ".txt")
    writeLines (content (response, encoding = "UTF-8"), log.file.name)
    # 4.1 Get network id (SUID) (CyRest function from Cytoscape)
    response <- GET (paste (cytoscape.base.url, "networks", "currentNetwork", sep = "/"))
    current.network.suid <- content (response, encode = "json")$data$networkSUID
    # print (current.network.suid)
    
    print (paste ("Save results", sep = ""))
    
    # Get network graphics (CyRest function from Cytoscape)
    image.type = "svg" # png, pdf
    response <- GET (paste (cytoscape.base.url, "networks", current.network.suid, "views", paste ("first.", image.type, sep = ""), sep = "/"))
    image.file.name = paste (output.folder.Diff, paste ("ClueGONetwork.", image.type, sep = ""), sep = "/")
    writeBin (content (response, encode = "raw"), image.file.name)
    
    # 4.2 Get ClueGO result table
    response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-table", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Result-Table.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name,
                 row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.3 Get ClueGO genes and main functions
    # number.of.functions.to.add = 3
    # response <- GET (paste (cluego.base.url, "analysis-results", "get-main-functions", 
    # current.network.suid, number.of.functions.to.add, sep = "/"))
    # result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    # table.file.name = paste (output.folder.Diff, "ClueGO-Genes-With-Main-Functions.txt", sep = "/")
    # write.table (text.to.data.frame (result.table.text), file = table.file.name,
    # row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.4 Get genes result table
    # response <- GET (paste (cluego.base.url, "analysis-results", "get-gene-table", current.network.suid, sep = "/"))
    # result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    # table.file.name = paste (output.folder.M, "ClueGO-Gene-Table.txt", sep = "/")
    # write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.5 Get Kappascore Matrix
    response <- GET (paste (cluego.base.url, "analysis-results", "get-kappascore-matrix", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Kappascore-Matrix.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, 
                 quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.6 Get binary Gene-Term Matrix
    response <- GET (paste (cluego.base.url, "analysis-results", "get-binary-gene-term-matrix", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Binary-Gene-Term-Matrix.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, 
                 quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.7 ClueGO Result Chart
    # Get result charts for both cluster as pie chart
    chart.type = "BarChart" # ("PieChart", "BarChart")
    image.type = "svg" # ("svg", "png", "pdf")
    response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-result-chart", 
                            current.network.suid, cluster1, chart.type, image.type, sep = "/"))
    image.file.name = paste (output.folder.Diff, paste ("ClueGO-", chart.type, "-For-Cluster", cluster1, ".", 
                                                     image.type, sep = ""), sep = "/")
    writeBin (content (response, encode = "raw"), image.file.name)
    ## [optional functions and settings, un comment and modify if needed]
    # 4.8 Remove ClueGO analysis result
    # print (paste ("Remove ClueGO Network", sep = ""))
    # # Remove analysis to reduce memory usage. This is important when using patch modes that create lots of analyses.
    # response <- DELETE (paste (cluego.base.url, "remove-cluego-analysis-result", current.network.suid, sep = "/"))

    # 4.9 Remove all ClueGO analysis results
    # print (paste ("Remove all ClueGO Networks", sep = ""))
    # response <- DELETE (paste (cluego.base.url, "remove-all-cluego-analysis-results", sep = "/"))

    print (paste (Enrichment.Analysis, DiffCo, "done", sep = " "))

    rm (list = setdiff (ls(), c ("cluego.base.url", "cluego.home.folder", "cytoscape.base.url",
                                "Enrichment.Analysis", "home.folder", "host.address", "output.folder",
                                "port.number", "text.to.data.frame", "DiffCoExpGenes" )))
  }
  
  
}

filenames <- list.files (path = "OurDataset/MuscleDataset/Outputs/ClueGO/DifferentialCoexpression",
                         recursive = TRUE, pattern = "Result-Table.txt", full.names = T)

ClueGO_enrichment <- c ()
for (file in filenames){
  ClueGO <- read.csv (file, sep = "\t")
  ClueGO$DC <- gsub (".*DifferentialCoexpression/|/ClueGO.*", "", file)
  ClueGO_enrichment <- rbind (ClueGO_enrichment, ClueGO)
}

ClueGO_enrichment <- ClueGO_enrichment [, c (14, 1:13)]

write.csv (ClueGO_enrichment, "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableDiffCoExpresses.csv", row.names = F)
save (ClueGO_enrichment, file = "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableDiffCoExpresses.RData")

# DC sub-networks
dir.create (paste (output.folder, "DifferentialCoexpressionSubNetworks", sep = "/"))
load ("Outputs/NetworkPerMuscleWGCNA/SubNetworks/CytoscapeInput-TOM-nodes_TO.075_kTotal.5.RData")
head(SubNetworkNodes)
for (Mus in c ("GAL", "GRA", "SEM", "SED", "REF", "VAM", "VAL")) {
  response <- POST ("http://localhost:1234/v1/apps/cluego/start-up-cluego", encode = "json")
  #### 0.0 Start up ClueGO in case it is not running yet
  response <- POST (url = paste (cytoscape.base.url, "apps", "cluego", "start-up-cluego", sep = "/"),
                    encode = "json")
  # wait 2 seconds to make sure ClueGO is started
  if (grepl ("HTTP ERROR 500", response)) {
    print ("wait 2 secs")
    Sys.sleep (2)
  }
  
  #### 1.0 Select the ClueGO Organism to analyze ####
  organism.name = "Homo Sapiens" # (run "1.1 Get all ClueGO organisms" to get all options)
  print (paste ("1.0 Select the ClueGO Organism to analyze: ", organism.name, sep = ""))
  response <- PUT (url = paste (cluego.base.url, "organisms", "set-organism", URLencode (organism.name), 
                                sep = "/"), encode = "json")
  #### 2.0 Upload IDs for a specific list ####
  print (paste ("2.0 Upload IDs for", Mus, sep = " "))
  
  ## [optional functions and settings, un comment and modify if needed]
  # 2.1 Select the ClueGO ID type
  id.type.name = "EnsemblGeneID" # (run "2.3 Get ClueGO ID types" to get all options)
  response <- PUT (url = paste (cluego.base.url, "ids", "set-id-type", id.type.name, sep = "/"), encode = "json")
  
  # 2.4 Set the number of modules
  max.input.panel.number = 1
  response <- PUT (url = paste (cluego.base.url, "cluster", "max-input-panel", max.input.panel.number, sep = "/"), encode = "json")
  
  
  ## [optional functions and settings, un comment and modify if needed]
  # 2.5 Set analysis properties for a module, this is to repeat for each input module
  cluster1 = "1"
  
  input.panel.index = cluster1 # set here the module number e.g. "1"
  node.shape = "Ellipse" # ("Ellipse", "Diamond", "Hexagon", "Octagon", "Parallelogram", "Rectangle", "Round Rectangle", "Triangle", "V")
  cluster.color = "#ff0000" # The color in hex, e.g. #F3A455
  min.number.of.genes.per.term = 4
  min.percentage.of.genes.mapped = 4
  no.restrictions = FALSE # TRUE for no restricions in number and percentage per term
  response <- PUT (url = paste (cluego.base.url, "cluster", "set-analysis-properties", input.panel.index, node.shape, URLencode (cluster.color, reserved = TRUE), min.number.of.genes.per.term, min.percentage.of.genes.mapped, no.restrictions, sep = "/"), encode = "json")
  
  # set gene list for cluster 1
  gene.list1 <- toJSON (SubNetworkNodes$nodeName [SubNetworkNodes$Networks == Mus])
  response <- PUT (url = paste (cluego.base.url, "cluster", "upload-ids-list", URLencode (cluster1), sep = "/"), body = gene.list1, encode = "json", content_type_json ())
  # 2.6 Select visual style
  visual.style = "ShowGroupDifference"  # (ShowGroupDifference, ShowSignificanceDifference, ShowClusterDifference)
  response <- PUT (url = paste (cluego.base.url, "cluster", "select-visual-style", visual.style, sep = "/"), encode = "json")
  ####  3.0 Select Ontologies
  print (paste ("3.0 Select Ontologies", sep = ""))
  selected.ontologies <- toJSON (c ("3;Ellipse", "4;Hexagon", "6;Octagon", 
                                    "8;Triangle", "9;Rectangle", "11;Round Rectangle")) # (run "3.1 Get all available Ontologies" to get all options)
  response <- PUT (url = paste (cluego.base.url, "ontologies", "set-ontologies", sep = "/"), body = selected.ontologies, encode = "json", content_type_json ())

  # 3.4 Set min, max GO tree level
  min.level = 5
  max.level = 6
  all.levels = FALSE
  response <- PUT (url = paste (cluego.base.url, "ontologies", "set-min-max-levels", min.level, max.level, all.levels, sep = "/"), encode = "json")
  # # 3.5 Use GO term significance cutoff
  p.value.cutoff = 0.05
  use.significance.cutoff = TRUE
  response <- PUT (url = paste (cluego.base.url, "ontologies", use.significance.cutoff, p.value.cutoff, sep = "/"), encode = "json")
  
  # # 3.7 Set statistical parameters
  # ("Enrichment (Right-sided hypergeometric test)", "Depletion (Left-sided hypergeometric test)", "Enrichment/Depletion (Two-sided hypergeometric test)")
  enrichment.type = "Enrichment (Right-sided hypergeometric test)" 
  multiple.testing.correction = "Benjamini-Hochberg"
  use.mid.pvalues = FALSE
  use.doubling = FALSE
  use.custom.reference.set = TRUE
  custom.reference.set = paste (output.folder, "GenesModules/All.txt", sep = "/")
  response <- PUT (url = paste (cluego.base.url, "stats", enrichment.type, multiple.testing.correction, 
                                use.mid.pvalues, use.doubling, use.custom.reference.set, custom.reference.set,
                                sep = "/"), encode = "json")
  #### 4.0 Run ClueGO Analysis ####
  print (paste ("4.0 Run ClueGO Analysis", sep = ""))
  # Run the analysis an save log file
  analysis.name <- Mus
  analysis.option <- "Cancel and refine selection" # ("Continue analysis", "Skip the grouping", "Cancel and refine selection")  -> Analysis option in case there are more than 1000 terms found!
  response <- GET (paste (cluego.base.url, URLencode (analysis.name), URLencode (analysis.option), sep = "/"))
  if (!grepl ("error", http_status (response)["category"])) {
    output.folder.Mus <- paste (output.folder, "DifferentialCoexpressionSubNetworks", Mus, sep = "/")
    dir.create (output.folder.Mus)
    log.file.name = paste0 (output.folder.Mus, "/", "ClueGOlog", Mus, ".txt")
    writeLines (content (response, encoding = "UTF-8"), log.file.name)
    # 4.1 Get network id (SUID) (CyRest function from Cytoscape)
    response <- GET (paste (cytoscape.base.url, "networks", "currentNetwork", sep = "/"))
    current.network.suid <- content (response, encode = "json")$data$networkSUID
    # print (current.network.suid)
    
    print (paste ("Save results", sep = ""))
    
    # Get network graphics (CyRest function from Cytoscape)
    image.type = "svg" # png, pdf
    response <- GET (paste (cytoscape.base.url, "networks", current.network.suid, "views", paste ("first.", image.type, sep = ""), sep = "/"))
    image.file.name = paste (output.folder.Mus, paste ("ClueGONetwork.", image.type, sep = ""), sep = "/")
    writeBin (content (response, encode = "raw"), image.file.name)
    
    # 4.2 Get ClueGO result table
    response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-table", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Mus, "ClueGO-Result-Table.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name,
                 row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.3 Get ClueGO genes and main functions
    # number.of.functions.to.add = 3
    # response <- GET (paste (cluego.base.url, "analysis-results", "get-main-functions", 
    # current.network.suid, number.of.functions.to.add, sep = "/"))
    # result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    # table.file.name = paste (output.folder.Diff, "ClueGO-Genes-With-Main-Functions.txt", sep = "/")
    # write.table (text.to.data.frame (result.table.text), file = table.file.name,
    # row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.4 Get genes result table
    # response <- GET (paste (cluego.base.url, "analysis-results", "get-gene-table", current.network.suid, sep = "/"))
    # result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    # table.file.name = paste (output.folder.M, "ClueGO-Gene-Table.txt", sep = "/")
    # write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.5 Get Kappascore Matrix
    response <- GET (paste (cluego.base.url, "analysis-results", "get-kappascore-matrix", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Mus, "ClueGO-Kappascore-Matrix.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, 
                 quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.6 Get binary Gene-Term Matrix
    response <- GET (paste (cluego.base.url, "analysis-results", "get-binary-gene-term-matrix", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Mus, "ClueGO-Binary-Gene-Term-Matrix.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, 
                 quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.7 ClueGO Result Chart
    # Get result charts for both cluster as pie chart
    chart.type = "BarChart" # ("PieChart", "BarChart")
    image.type = "svg" # ("svg", "png", "pdf")
    response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-result-chart", 
                            current.network.suid, cluster1, chart.type, image.type, sep = "/"))
    image.file.name = paste (output.folder.Mus, paste ("ClueGO-", chart.type, "-For-Cluster", cluster1, ".", 
                                                        image.type, sep = ""), sep = "/")
    writeBin (content (response, encode = "raw"), image.file.name)
    ## [optional functions and settings, un comment and modify if needed]
    # 4.8 Remove ClueGO analysis result
    # print (paste ("Remove ClueGO Network", sep = ""))
    # # Remove analysis to reduce memory usage. This is important when using patch modes that create lots of analyses.
    # response <- DELETE (paste (cluego.base.url, "remove-cluego-analysis-result", current.network.suid, sep = "/"))
    
    # 4.9 Remove all ClueGO analysis results
    # print (paste ("Remove all ClueGO Networks", sep = ""))
    # response <- DELETE (paste (cluego.base.url, "remove-all-cluego-analysis-results", sep = "/"))
    
    print (paste (Enrichment.Analysis, Mus, "done", sep = " "))
    
    rm (list = setdiff (ls(), c ("cluego.base.url", "cluego.home.folder", "cytoscape.base.url",
                                 "Enrichment.Analysis", "home.folder", "host.address", "output.folder",
                                 "port.number", "text.to.data.frame", "SubNetworkNodes" )))
  }
  
}

filenames <- list.files (path = "OurDataset/MuscleDataset/Outputs/ClueGO/DifferentialCoexpressionSubNetworks",
                         recursive = TRUE, pattern = "Result-Table.txt", full.names = T)

ClueGO_enrichment <- c ()
for (file in filenames){
  ClueGO <- read.csv (file, sep = "\t")
  ClueGO$DC <- gsub (".*DifferentialCoexpressionSubNetworks/|/ClueGO.*", "", file)
  ClueGO_enrichment <- rbind (ClueGO_enrichment, ClueGO)
}

ClueGO_enrichment <- ClueGO_enrichment [, c (14, 1:13)]

write.csv (ClueGO_enrichment, "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableDiffCoExpressesSubNetworks.csv", row.names = F)
save (ClueGO_enrichment, file = "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableDiffCoExpressesSubNetworks.RData")


# DE 
# rm (list = ls())
dir.create (paste (output.folder, "DifferentialExpression", sep = "/"))

load ("Outputs/CellTypeComposition/DEGeneNotRelatedToCellTypes.RData")
for (DiffDE in names (DEnotCellTypeRelated)[42]){
  response <- POST ("http://localhost:1234/v1/apps/cluego/start-up-cluego", encode = "json")
  
  #### 0.0 Start up ClueGO in case it is not running yet
  response <- POST (url = paste (cytoscape.base.url, "apps", "cluego", "start-up-cluego", sep = "/"),
                    encode = "json")
  # wait 2 seconds to make sure ClueGO is started
  if (grepl ("HTTP ERROR 500", response)) {
    print ("wait 2 secs")
    Sys.sleep (2)
  }
  
  #### 1.0 Select the ClueGO Organism to analyze ####
  organism.name = "Homo Sapiens" # (run "1.1 Get all ClueGO organisms" to get all options)
  print (paste ("1.0 Select the ClueGO Organism to analyze: ", organism.name, sep = ""))
  response <- PUT (url = paste (cluego.base.url, "organisms", "set-organism", URLencode (organism.name), 
                                sep = "/"), encode = "json")
  
  #### 2.0 Upload IDs for a specific list ####
  print (paste ("2.0 Upload IDs for", DiffDE, sep = " "))
  
  ## [optional functions and settings, un comment and modify if needed]
  # 2.1 Select the ClueGO ID type
  id.type.name = "EnsemblGeneID" # (run "2.3 Get ClueGO ID types" to get all options)
  response <- PUT (url = paste (cluego.base.url, "ids", "set-id-type", id.type.name, sep = "/"), encode = "json")
  
  # 2.4 Set the number of modules
  max.input.panel.number = 1
  response <- PUT (url = paste (cluego.base.url, "cluster", "max-input-panel", max.input.panel.number, sep = "/"), encode = "json")
  
  
  ## [optional functions and settings, un comment and modify if needed]
  # 2.5 Set analysis properties for a module, this is to repeat for each input module
  cluster1 = "1"
  
  input.panel.index = cluster1 # set here the module number e.g. "1"
  node.shape = "Ellipse" # ("Ellipse", "Diamond", "Hexagon", "Octagon", "Parallelogram", "Rectangle", "Round Rectangle", "Triangle", "V")
  cluster.color = "#ff0000" # The color in hex, e.g. #F3A455
  min.number.of.genes.per.term = 4
  min.percentage.of.genes.mapped = 4
  no.restrictions = FALSE # TRUE for no restricions in number and percentage per term
  response <- PUT (url = paste (cluego.base.url, "cluster", "set-analysis-properties", input.panel.index, node.shape, URLencode (cluster.color, reserved = TRUE), min.number.of.genes.per.term, min.percentage.of.genes.mapped, no.restrictions, sep = "/"), encode = "json")
  
  # set gene list for cluster 1
  gene.list1 <- toJSON (DEnotCellTypeRelated [[DiffDE]])
  response <- PUT (url = paste (cluego.base.url, "cluster", "upload-ids-list", URLencode (cluster1), sep = "/"), body = gene.list1, encode = "json", content_type_json ())
  # 2.6 Select visual style
  visual.style = "ShowGroupDifference"  # (ShowGroupDifference, ShowSignificanceDifference, ShowClusterDifference)
  response <- PUT (url = paste (cluego.base.url, "cluster", "select-visual-style", visual.style, sep = "/"), encode = "json")
  ####  3.0 Select Ontologies
  print (paste ("3.0 Select Ontologies", sep = ""))
  selected.ontologies <- toJSON (c ("3;Ellipse", "4;Hexagon", "6;Octagon", 
                                    "8;Triangle", "9;Rectangle", "11;Round Rectangle")) # (run "3.1 Get all available Ontologies" to get all options)
  response <- PUT (url = paste (cluego.base.url, "ontologies", "set-ontologies", sep = "/"), body = selected.ontologies, encode = "json", content_type_json ())
  
  # 3.4 Set min, max GO tree level
  min.level = 5
  max.level = 6
  all.levels = FALSE
  response <- PUT (url = paste (cluego.base.url, "ontologies", "set-min-max-levels", min.level, max.level, all.levels, sep = "/"), encode = "json")
  # # 3.5 Use GO term significance cutoff
  p.value.cutoff = 0.05
  use.significance.cutoff = TRUE
  response <- PUT (url = paste (cluego.base.url, "ontologies", use.significance.cutoff, p.value.cutoff, sep = "/"), encode = "json")
  
  # # 3.7 Set statistical parameters
  # ("Enrichment (Right-sided hypergeometric test)", "Depletion (Left-sided hypergeometric test)", "Enrichment/Depletion (Two-sided hypergeometric test)")
  enrichment.type = "Enrichment (Right-sided hypergeometric test)" 
  multiple.testing.correction = "Benjamini-Hochberg"
  use.mid.pvalues = FALSE
  use.doubling = FALSE
  use.custom.reference.set = TRUE
  custom.reference.set = paste (output.folder, "GenesModules/All.txt", sep = "/")
  response <- PUT (url = paste (cluego.base.url, "stats", enrichment.type, multiple.testing.correction, 
                                use.mid.pvalues, use.doubling, use.custom.reference.set, custom.reference.set,
                                sep = "/"), encode = "json")
  #### 4.0 Run ClueGO Analysis ####
  print (paste ("4.0 Run ClueGO Analysis", sep = ""))
  # Run the analysis an save log file
  analysis.name <- DiffDE
  analysis.option <- "Cancel and refine selection" # ("Continue analysis", "Skip the grouping", "Cancel and refine selection")  -> Analysis option in case there are more than 1000 terms found!
  response <- GET (paste (cluego.base.url, URLencode (analysis.name), URLencode (analysis.option), sep = "/"))
  if (!grepl ("error", http_status ( response)["category"])) {
    output.folder.Diff <- paste (output.folder, "DifferentialExpression", DiffDE, sep = "/")
    dir.create (output.folder.Diff)
    log.file.name = paste0 (output.folder.Diff, "/", "ClueGOlog", DiffDE, ".txt")
    writeLines (content (response, encoding = "UTF-8"), log.file.name)
    # 4.1 Get network id (SUID) (CyRest function from Cytoscape)
    response <- GET (paste (cytoscape.base.url, "networks", "currentNetwork", sep = "/"))
    current.network.suid <- content (response, encode = "json")$data$networkSUID
    # print (current.network.suid)
    
    print (paste ("Save results", sep = ""))
    
    # Get network graphics (CyRest function from Cytoscape)
    image.type = "svg" # png, pdf
    response <- GET (paste (cytoscape.base.url, "networks", current.network.suid, "views", paste ("first.", image.type, sep = ""), sep = "/"))
    image.file.name = paste (output.folder.Diff, paste ("ClueGONetwork.", image.type, sep = ""), sep = "/")
    writeBin (content (response, encode = "raw"), image.file.name)
    
    # 4.2 Get ClueGO result table
    response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-table", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Result-Table.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name,
                 row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.3 Get ClueGO genes and main functions
    # number.of.functions.to.add = 3
    # response <- GET (paste (cluego.base.url, "analysis-results", "get-main-functions", 
    # current.network.suid, number.of.functions.to.add, sep = "/"))
    # result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    # table.file.name = paste (output.folder.Diff, "ClueGO-Genes-With-Main-Functions.txt", sep = "/")
    # write.table (text.to.data.frame (result.table.text), file = table.file.name,
    # row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.4 Get genes result table
    # response <- GET (paste (cluego.base.url, "analysis-results", "get-gene-table", current.network.suid, sep = "/"))
    # result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    # table.file.name = paste (output.folder.M, "ClueGO-Gene-Table.txt", sep = "/")
    # write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.5 Get Kappascore Matrix
    response <- GET (paste (cluego.base.url, "analysis-results", "get-kappascore-matrix", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Kappascore-Matrix.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, 
                 quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.6 Get binary Gene-Term Matrix
    response <- GET (paste (cluego.base.url, "analysis-results", "get-binary-gene-term-matrix", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Binary-Gene-Term-Matrix.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, 
                 quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.7 ClueGO Result Chart
    # Get result charts for both cluster as pie chart
    chart.type = "BarChart" # ("PieChart", "BarChart")
    image.type = "svg" # ("svg", "png", "pdf")
    response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-result-chart", 
                            current.network.suid, cluster1, chart.type, image.type, sep = "/"))
    image.file.name = paste (output.folder.Diff, paste ("ClueGO-", chart.type, "-For-Cluster", cluster1, ".", 
                                                        image.type, sep = ""), sep = "/")
    writeBin (content (response, encode = "raw"), image.file.name)
    ## [optional functions and settings, un comment and modify if needed]
    # 4.8 Remove ClueGO analysis result
    # print (paste ("Remove ClueGO Network", sep = ""))
    # # Remove analysis to reduce memory usage. This is important when using patch modes that create lots of analyses.
    # response <- DELETE (paste (cluego.base.url, "remove-cluego-analysis-result", current.network.suid, sep = "/"))
    
    # 4.9 Remove all ClueGO analysis results
    print (paste ("Remove all ClueGO Networks", sep = ""))
    # response <- DELETE (paste (cluego.base.url, "remove-all-cluego-analysis-results", sep = "/"))
    
    print (paste (Enrichment.Analysis, DiffDE, "done", sep = " "))
    
    rm (list = setdiff (ls(), c ("cluego.base.url", "cluego.home.folder", "cytoscape.base.url",
                                 "Enrichment.Analysis", "home.folder", "host.address", "output.folder",
                                 "port.number", "text.to.data.frame", "DEnotCellTypeRelated" )))
  }
  
  
}

filenames <- list.files (path = "OurDataset/MuscleDataset/Outputs/ClueGO/DifferentialExpression",
                         recursive = TRUE, pattern = "Result-Table.txt", full.names = T)

ClueGO_enrichment <- c ()
for (file in filenames){
  ClueGO <- read.csv (file, sep = "\t")
  ClueGO$DC <- gsub (".*DifferentialExpression/|/ClueGO.*", "", file)
  ClueGO_enrichment <- rbind (ClueGO_enrichment, ClueGO)
}

ClueGO_enrichment <- ClueGO_enrichment [, c (14, 1:13)]

write.csv (ClueGO_enrichment, "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableDifferentialExpression.csv", row.names = F)
save (ClueGO_enrichment, file = "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableDifferentialExpression.RData")


## LFC1
dir.create (paste (output.folder, "DifferentialExpressionLFC1", sep = "/"))

load ("Outputs/CellTypeComposition/DEGeneNotRelatedToCellTypes.RData")
for (DiffDE in names (DELFC1_notCellTypeRelated)){
  response <- POST ("http://localhost:1234/v1/apps/cluego/start-up-cluego", encode = "json")
  
  #### 0.0 Start up ClueGO in case it is not running yet
  response <- POST (url = paste (cytoscape.base.url, "apps", "cluego", "start-up-cluego", sep = "/"),
                    encode = "json")
  # wait 2 seconds to make sure ClueGO is started
  if (grepl ("HTTP ERROR 500", response)) {
    print ("wait 2 secs")
    Sys.sleep (2)
  }
  
  #### 1.0 Select the ClueGO Organism to analyze ####
  organism.name = "Homo Sapiens" # (run "1.1 Get all ClueGO organisms" to get all options)
  print (paste ("1.0 Select the ClueGO Organism to analyze: ", organism.name, sep = ""))
  response <- PUT (url = paste (cluego.base.url, "organisms", "set-organism", URLencode (organism.name), 
                                sep = "/"), encode = "json")
  
  #### 2.0 Upload IDs for a specific list ####
  print (paste ("2.0 Upload IDs for", DiffDE, sep = " "))
  
  ## [optional functions and settings, un comment and modify if needed]
  # 2.1 Select the ClueGO ID type
  id.type.name = "EnsemblGeneID" # (run "2.3 Get ClueGO ID types" to get all options)
  response <- PUT (url = paste (cluego.base.url, "ids", "set-id-type", id.type.name, sep = "/"), encode = "json")
  
  # 2.4 Set the number of modules
  max.input.panel.number = 1
  response <- PUT (url = paste (cluego.base.url, "cluster", "max-input-panel", max.input.panel.number, sep = "/"), encode = "json")
  
  
  ## [optional functions and settings, un comment and modify if needed]
  # 2.5 Set analysis properties for a module, this is to repeat for each input module
  cluster1 = "1"
  
  input.panel.index = cluster1 # set here the module number e.g. "1"
  node.shape = "Ellipse" # ("Ellipse", "Diamond", "Hexagon", "Octagon", "Parallelogram", "Rectangle", "Round Rectangle", "Triangle", "V")
  cluster.color = "#ff0000" # The color in hex, e.g. #F3A455
  min.number.of.genes.per.term = 4
  min.percentage.of.genes.mapped = 4
  no.restrictions = FALSE # TRUE for no restricions in number and percentage per term
  response <- PUT (url = paste (cluego.base.url, "cluster", "set-analysis-properties", input.panel.index, node.shape, URLencode (cluster.color, reserved = TRUE), min.number.of.genes.per.term, min.percentage.of.genes.mapped, no.restrictions, sep = "/"), encode = "json")
  
  # set gene list for cluster 1
  gene.list1 <- toJSON (DELFC1_notCellTypeRelated [[DiffDE]])
  response <- PUT (url = paste (cluego.base.url, "cluster", "upload-ids-list", URLencode (cluster1), sep = "/"), body = gene.list1, encode = "json", content_type_json ())
  # 2.6 Select visual style
  visual.style = "ShowGroupDifference"  # (ShowGroupDifference, ShowSignificanceDifference, ShowClusterDifference)
  response <- PUT (url = paste (cluego.base.url, "cluster", "select-visual-style", visual.style, sep = "/"), encode = "json")
  ####  3.0 Select Ontologies
  print (paste ("3.0 Select Ontologies", sep = ""))
  selected.ontologies <- toJSON (c ("3;Ellipse", "4;Hexagon", "6;Octagon", 
                                    "8;Triangle", "9;Rectangle", "11;Round Rectangle")) # (run "3.1 Get all available Ontologies" to get all options)
  response <- PUT (url = paste (cluego.base.url, "ontologies", "set-ontologies", sep = "/"), body = selected.ontologies, encode = "json", content_type_json ())
  
  # 3.4 Set min, max GO tree level
  min.level = 5
  max.level = 6
  all.levels = FALSE
  response <- PUT (url = paste (cluego.base.url, "ontologies", "set-min-max-levels", min.level, max.level, all.levels, sep = "/"), encode = "json")
  # # 3.5 Use GO term significance cutoff
  p.value.cutoff = 0.05
  use.significance.cutoff = TRUE
  response <- PUT (url = paste (cluego.base.url, "ontologies", use.significance.cutoff, p.value.cutoff, sep = "/"), encode = "json")
  
  # # 3.7 Set statistical parameters
  # ("Enrichment (Right-sided hypergeometric test)", "Depletion (Left-sided hypergeometric test)", "Enrichment/Depletion (Two-sided hypergeometric test)")
  enrichment.type = "Enrichment (Right-sided hypergeometric test)" 
  multiple.testing.correction = "Benjamini-Hochberg"
  use.mid.pvalues = FALSE
  use.doubling = FALSE
  use.custom.reference.set = TRUE
  custom.reference.set = paste (output.folder, "GenesModules/All.txt", sep = "/")
  response <- PUT (url = paste (cluego.base.url, "stats", enrichment.type, multiple.testing.correction, 
                                use.mid.pvalues, use.doubling, use.custom.reference.set, custom.reference.set,
                                sep = "/"), encode = "json")
  #### 4.0 Run ClueGO Analysis ####
  print (paste ("4.0 Run ClueGO Analysis", sep = ""))
  # Run the analysis an save log file
  analysis.name <- DiffDE
  analysis.option <- "Cancel and refine selection" # ("Continue analysis", "Skip the grouping", "Cancel and refine selection")  -> Analysis option in case there are more than 1000 terms found!
  response <- GET (paste (cluego.base.url, URLencode (analysis.name), URLencode (analysis.option), sep = "/"))
  if (!grepl ("error", http_status ( response)["category"])) {
    output.folder.Diff <- paste (output.folder, "DifferentialExpressionLFC1", DiffDE, sep = "/")
    dir.create (output.folder.Diff)
    log.file.name = paste0 (output.folder.Diff, "/", "ClueGOlog", DiffDE, ".txt")
    writeLines (content (response, encoding = "UTF-8"), log.file.name)
    # 4.1 Get network id (SUID) (CyRest function from Cytoscape)
    response <- GET (paste (cytoscape.base.url, "networks", "currentNetwork", sep = "/"))
    current.network.suid <- content (response, encode = "json")$data$networkSUID
    # print (current.network.suid)
    
    print (paste ("Save results", sep = ""))
    
    # Get network graphics (CyRest function from Cytoscape)
    image.type = "svg" # png, pdf
    response <- GET (paste (cytoscape.base.url, "networks", current.network.suid, "views", paste ("first.", image.type, sep = ""), sep = "/"))
    image.file.name = paste (output.folder.Diff, paste ("ClueGONetwork.", image.type, sep = ""), sep = "/")
    writeBin (content (response, encode = "raw"), image.file.name)
    
    # 4.2 Get ClueGO result table
    response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-table", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Result-Table.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name,
                 row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.3 Get ClueGO genes and main functions
    # number.of.functions.to.add = 3
    # response <- GET (paste (cluego.base.url, "analysis-results", "get-main-functions", 
    # current.network.suid, number.of.functions.to.add, sep = "/"))
    # result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    # table.file.name = paste (output.folder.Diff, "ClueGO-Genes-With-Main-Functions.txt", sep = "/")
    # write.table (text.to.data.frame (result.table.text), file = table.file.name,
    # row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.4 Get genes result table
    # response <- GET (paste (cluego.base.url, "analysis-results", "get-gene-table", current.network.suid, sep = "/"))
    # result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    # table.file.name = paste (output.folder.M, "ClueGO-Gene-Table.txt", sep = "/")
    # write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.5 Get Kappascore Matrix
    response <- GET (paste (cluego.base.url, "analysis-results", "get-kappascore-matrix", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Kappascore-Matrix.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, 
                 quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.6 Get binary Gene-Term Matrix
    response <- GET (paste (cluego.base.url, "analysis-results", "get-binary-gene-term-matrix", 
                            current.network.suid, sep = "/"))
    result.table.text <- content (response, encode = "text", encoding = "UTF-8")
    table.file.name = paste (output.folder.Diff, "ClueGO-Binary-Gene-Term-Matrix.txt", sep = "/")
    write.table (text.to.data.frame (result.table.text), file = table.file.name, row.names = FALSE, 
                 quote = FALSE, na = "", col.names = TRUE, sep = "\t")
    
    # 4.7 ClueGO Result Chart
    # Get result charts for both cluster as pie chart
    chart.type = "BarChart" # ("PieChart", "BarChart")
    image.type = "svg" # ("svg", "png", "pdf")
    response <- GET (paste (cluego.base.url, "analysis-results", "get-cluego-result-chart", 
                            current.network.suid, cluster1, chart.type, image.type, sep = "/"))
    image.file.name = paste (output.folder.Diff, paste ("ClueGO-", chart.type, "-For-Cluster", cluster1, ".", 
                                                        image.type, sep = ""), sep = "/")
    writeBin (content (response, encode = "raw"), image.file.name)
    ## [optional functions and settings, un comment and modify if needed]
    # 4.8 Remove ClueGO analysis result
    # print (paste ("Remove ClueGO Network", sep = ""))
    # # Remove analysis to reduce memory usage. This is important when using patch modes that create lots of analyses.
    # response <- DELETE (paste (cluego.base.url, "remove-cluego-analysis-result", current.network.suid, sep = "/"))
    
    # 4.9 Remove all ClueGO analysis results
    print (paste ("Remove all ClueGO Networks", sep = ""))
    # response <- DELETE (paste (cluego.base.url, "remove-all-cluego-analysis-results", sep = "/"))
    
    print (paste (Enrichment.Analysis, DiffDE, "done", sep = " "))
    
    rm (list = setdiff (ls(), c ("cluego.base.url", "cluego.home.folder", "cytoscape.base.url",
                                 "Enrichment.Analysis", "home.folder", "host.address", "output.folder",
                                 "port.number", "text.to.data.frame", "DELFC1_notCellTypeRelated" )))
  }
  
  
}

filenames <- list.files (path = "OurDataset/MuscleDataset/Outputs/ClueGO/DifferentialExpressionLFC1",
                         recursive = TRUE, pattern = "Result-Table.txt", full.names = T)

ClueGO_enrichment <- c ()
for (file in filenames){
  ClueGO <- read.csv (file, sep = "\t")
  ClueGO$DC <- gsub (".*DifferentialExpressionLFC1/|/ClueGO.*", "", file)
  ClueGO_enrichment <- rbind (ClueGO_enrichment, ClueGO)
}

ClueGO_enrichment <- ClueGO_enrichment [, c (14, 1:13)]

write.csv (ClueGO_enrichment, "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableDifferentialExpressionLFC1.csv", row.names = F)
save (ClueGO_enrichment, file = "OurDataset/MuscleDataset/Outputs/ClueGO/ClueGO-ResultTableDifferentialExpressionLFC1.RData")

