# devtools::install_github("felixfan/PubMedWordcloud")

thisyear <- as.numeric(substr(Sys.Date(),1,4))
IDs1 <- PubMedWordcloud::getPMIDs(author="Geoffrey Zahn",dFrom = 2011,dTo = thisyear)
IDs2 <- PubMedWordcloud::getPMIDs(author="Geoffrey L Zahn",dFrom = 2011,dTo = thisyear)
IDs <- PubMedWordcloud::editPMIDs(IDs1,IDs2,method = 'exclude')
abstracts <- PubMedWordcloud::getAbstracts(IDs)
cleanABS <- PubMedWordcloud::cleanAbstracts(abstracts)
cleanABS[cleanABS$word == "plants",]$freq <- 
  cleanABS[cleanABS$word == "plants",]$freq + 
  cleanABS[cleanABS$word == "plant",]$freq

cleanABS[cleanABS$word == "plant",]$freq <- 0
PubMedWordcloud::plotWordCloud(cleanABS)
