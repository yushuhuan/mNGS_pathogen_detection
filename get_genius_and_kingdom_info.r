# install.packages("dplyr")
# install.packages("stringr")
library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)


## 从终端获取参数
args <- commandArgs(trailingOnly = TRUE)
classified_nt_filename <- args[1]
output_dir <- args[2]
sample_id <- args[3]

# print(classified_nt_filename)

classified_nt <- read.table(classified_nt_filename, sep = "\t", header = F, stringsAsFactors = F, quote = "")
colnames(classified_nt) <- c("percent", "all_reads", "unque_reads", "tax_level", "tax_id", "name")
classified_nt$line_num <- 1:dim(classified_nt)[1]

### 提取species水平的病原信息
species_reads_info <- classified_nt %>% filter(grepl("^S", tax_level)) %>% select(name, tax_id, all_reads, percent)
species_reads_info$name <- str_trim(species_reads_info$name)

patho_genius_and_domain <- data.frame()
for(i in 1:dim(species_reads_info)[1]){
    id = as.character(species_reads_info[i,2])
    pathogen_name = str_trim(as.character(species_reads_info[i,1]))

    patho.line_num = which(classified_nt$tax_id == id)

    genius = classified_nt %>% filter(classified_nt$line_num < patho.line_num) %>% filter(tax_level == "G") %>% arrange(desc(line_num))
    genius_name = str_trim(as.character(genius$name[1]))

    domain = classified_nt %>% filter(classified_nt$line_num < patho.line_num) %>% filter(tax_level == "D") %>% arrange(desc(line_num))
    domain_name = str_trim(as.character(domain$name[1]))

    result = as.data.frame(t(c(domain_name, genius_name, pathogen_name)))
    colnames(result) = c("Domain", "Genius", "Species")

    patho_genius_and_domain = rbind(patho_genius_and_domain, result)
}

output <- merge(patho_genius_and_domain, species_reads_info, by.x = "Species", by.y = "name")

# 过滤掉人的分类信息
output <- output %>% filter(tax_id != "9606") %>% filter(tax_id != "9598")  %>% select(Domain, Genius, Species, tax_id, all_reads, percent)

# 计算基于domain的相对丰度
# domain_id <- c("2", "2759", "10239", "708403")
# names(domain_id) <- c("Bacteria", "Eukaryota", "Viruses", "Parasitus")

output$aboundance_based_domain <- unlist(apply(output, 1, function(x){
    if(x[1]=="Bacteria"){
        result = round((as.integer(x[5])*100)/classified_nt[which(classified_nt$tax_id == "2"), 2], 2)
    }else {
       if(x[1]=="Eukaryota"){
            result = round((as.integer(x[5])*100)/classified_nt[which(classified_nt$tax_id == "2759"), 2], 2)
       }else {
          if(x[1]=="Viruses"){
            result = round((as.integer(x[5])*100)/classified_nt[which(classified_nt$tax_id == "10239"), 2], 2)
          }else {
            if(x[1]=="Parasitus"){
             result = round((as.integer(x[5])*100)/classified_nt[which(classified_nt$tax_id == "708403"), 2], 2)
            }else {
               if(x[1]=="Archaea"){
                result = round((as.integer(x[5])*100)/classified_nt[which(classified_nt$tax_id == "2157"), 2], 2)
               }
            }
          }
       }
        
    }
}))

# 计算基于全部species序列的相对丰度
output$aboundance_based_species <- apply(output, 1, function(x){
    result = round((as.integer(x[5])*100)/classified_nt[which(classified_nt$tax_id == "1"), 2], 2)
})

# 计算species reads标准化
output$CPM <- apply(output, 1, function(x){
  all_reads = classified_nt[which(classified_nt$tax_id == "0"), 2] + classified_nt[which(classified_nt$tax_id == "1"), 2]
  result = round((as.integer(x[5])*1000000)/all_reads,2)
})

# 给出rank
output <- output %>% arrange(desc(all_reads))
output$Rank <- 1:dim(output)[1]

output <- output %>% select(Rank, tax_id, Domain, Genius, Species, all_reads, CPM, percent, aboundance_based_domain, aboundance_based_species)
output$aboundance_based_domain <- as.vector(output$aboundance_based_domain)
str(output)

write.table(output, paste0(output_dir, "/", sample_id, ".pathogen_with_genius_and_domain.txt"), row.names = F, sep = "\t", quote = F, fileEncoding = "UTF-8")
