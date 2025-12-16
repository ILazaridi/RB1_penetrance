## Adaptation of read_selfreport_cancer() function from "https://raw.githubusercontent.com/hdg204/UKBB/main/UKBB_Health_Records_New_Project.R"
# - Adapts the read_selfreport_cancer() function to generate read_selfreport_cancer_p20001_i1() and read_selfreport_cancer_p20001_i2() which look through instances 1 and 2 of the self reported data

## This is run within another script "phenotype_script_ustringent_permissive.R"

read_selfreport_cancer_p20001_i1 <- function(codes,file='selfreport_participant.csv'){
  data=read.csv(file)
  coding3=read.csv('coding3.tsv',sep='\t')%>%filter(coding>1)
  outlines=NULL
  for (i in 1:length(codes)){
    if (length(coding3[coding3$coding==codes[1],'meaning'])>0){
      outlines=c(outlines,grep(coding3[coding3$coding==codes[i],'meaning'],data$p20001_i1))
    }
  }
  data_frame=data.frame(eid=data[outlines,1])
  return(data_frame)
}

read_selfreport_cancer_p20001_i2 <- function(codes,file='selfreport_participant.csv'){
  data=read.csv(file)
  coding3=read.csv('coding3.tsv',sep='\t')%>%filter(coding>1)
  outlines=NULL
  for (i in 1:length(codes)){
    if (length(coding3[coding3$coding==codes[1],'meaning'])>0){
      outlines=c(outlines,grep(coding3[coding3$coding==codes[i],'meaning'],data$p20001_i2))
    }
  }
  data_frame=data.frame(eid=data[outlines,1])
  return(data_frame)
}