# install.packages("dplyr")
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_filename <- args[1]
NC_input_filename <- args[2]
output_dir <- args[3]
sample_id <- args[4]

input <- read.table(input_filename, header = T,sep = "\t")
NC_input <- read.table(NC_input_filename, header = T,sep = "\t")

colnames(NC_input)[3] <- "Reads_num_in_NC"
colnames(NC_input)[2] <- "taxid_NC"
colnames(NC_input)[1] <- "name"
NC_input <- NC_input %>% select(name, Reads_num_in_NC)

classified_with_NC <- merge(input, NC_input, by.x = "Species", by.y = "name", all.x = TRUE)
classified_with_NC[is.na(classified_with_NC$Reads_num_in_NC), 9] <- 0

classified_with_NC <- classified_with_NC %>% arrange(Rank)

### 获取致病性及常见背景微生物label
concerned_pathogens_V1 <- read.csv("/opt/mNGS/ZhiDe-mNGS-analysis-V3/known_pathogen_database_V1.0.csv", header = T, fill = TRUE)[,1:3]

labled_result <- merge(classified_with_NC, concerned_pathogens_V1, by.x = "Species", by.y = "species", all.x = TRUE)
labled_result <- labled_result %>% arrange(Rank)

extracted_result <- labled_result %>% filter(is.na(label)=="FALSE") %>% arrange(desc(all_reads))
all_pathogen_results <- labled_result

### 分配关注度
bacteria_attention <- function(species, reads_num, rank, genius, aboundance_based_species, reads_in_NC){
    # 如果是已知大概规则的病原，则可以进行标注，对于规则未知的病原，则不标注，不报出
    known_species = c("Cutibacterium acnes", "Staphylococcus epidermidis", "Staphylococcus aureus", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Bartonella henselae", "Treponema pallidum", "Pseudomonas aeruginosa", "Acinetobacter baumannii")
    known_genius = c("Enterococcus", "Enterobacter", "Nocardia", "Streptococcus", "Bacillus", "Escherichia", "Staphylococcus", "Acinetobacter", "Cutibacterium", "Moraxella", "Kocuria", "Corynebacterium")

    attention_level = ""

    if(species %in% known_species){
        if(species == "Cutibacterium acnes"){
            if(rank == 1){
                if(aboundance_based_species > 70){
                    attention_level = "高关注"
                }else{
                    attention_level = "背景微生物"
                }
            }else{
                attention_level = "背景微生物"
            }
        }
        if(species == "Staphylococcus epidermidis"){
            if(aboundance_based_species > extracted_result[which(extracted_result$Species == "Cutibacterium acnes"), 8]){
                attention_level = "高关注"
            }else{
                attention_level = "背景微生物"
            }
        }
        if(species == "Staphylococcus aureus"){
            if(aboundance_based_species > 10){
                attention_level = "低关注"
            }else{
                attention_level = "背景微生物"
            }
        }
        if(species == "Klebsiella pneumoniae"){if(rank == 1){attention_level = "高关注"}else{attention_level = "低关注"}}
        if(species %in% c("Mycobacterium tuberculosis", "Bartonella henselae", "Treponema pallidum")){if(reads_num < 20){attention_level = "低关注"}else{attention_level = "高关注"}}
        if(species %in% c("Pseudomonas aeruginosa", "Acinetobacter baumannii")){if(rank <= 3){attention_level = "结合临床及批次考虑"}else{attention_level = "背景微生物"}}
    }else {
       if(genius %in% known_genius){
        if(genius %in% c("Enterococcus", "Enterobacter", "Nocardia", "Streptococcus", "Bacillus", "Escherichia")){if(rank <= 5){attention_level = "高关注"}else{attention_level = "低关注"}}
        if(genius %in% c("Staphylococcus", "Acinetobacter", "Cutibacterium", "Moraxella", "Kocuria", "Corynebacterium")){if(rank ==1){attention_level = "高关注"}else{attention_level = "背景微生物"}}
       }else {
          attention_level = ""
       }
    }

    return(attention_level)
}

fungi_attention <- function(species, rank, genius){
    known_species = c("Komagataella pastoris", "Candida albicans", "Candida tropicalis", "Aspergillus fumigatus", "Aspergillus flavus", "Aspergillus niger", "Aspergillus nidulans", "Fusarium solani", "Fusarium verticillioides", "Fusarium proliferatum")
    
    attention_level = ""
    
    if(species %in% c(known_species)){
        if(species == "Komagataella pastoris"){attention_level = "污染"}
        if(species %in% c("Candida albicans", "Candida tropicalis", "Aspergillus fumigatus", "Aspergillus flavus", "Aspergillus niger", "Aspergillus nidulans", "Fusarium solani", "Fusarium verticillioides", "Fusarium proliferatum")){
            if(rank <= 3){
                attention_level = "高关注"
            }
        }
    }

    return(attention_level)
}

viruses_attention <- function(species, rank, label){
    known_species = c("Human alphaherpesvirus 1", "Human alphaherpesvirus 2", "Human alphaherpesvirus 3", "Human betaherpesvirus 5", "Human gammaherpesvirus 4", "Human betaherpesvirus 6", "Human mastadenovirus D", "Suid alphaherpesvirus 1", 
                    "Torque teno virus", "Torque teno virus 7", "Torque teno virus 10", "Torque teno virus 11", "Torque teno virus 12", "Torque teno virus 13", "Torque teno virus 16", "Torque teno virus 18", "Torque teno virus 19", "Torque teno virus 20", "Torque teno virus 22", "Torque teno virus 24", "Torque teno virus 27", "Torque teno virus 29",
                    "Betapapillomavirus 1", "Betapapillomavirus 2", "Betapapillomavirus 3", "Alphapolyomavirus octihominis", "Alphapolyomavirus quintihominis")
    
    attention_level = ""

    if(species %in% known_species){
        if(species %in% c("Human alphaherpesvirus 1", "Human alphaherpesvirus 2", "Human alphaherpesvirus 3", "Human betaherpesvirus 5")){attention_level = "结合批次检出情况确认"}
        if(species %in% c("Human gammaherpesvirus 4", "Human betaherpesvirus 6")){if(rank <=3){attention_level = "高关注"}else{attention_level = "低关注"}}
        if(species %in% c("Human mastadenovirus D", "Suid alphaherpesvirus 1")){attention_level = "结合批次检出情况确认"}
        if(species %in% c("Torque teno virus", "Torque teno virus 7", "Torque teno virus 10", "Torque teno virus 11", "Torque teno virus 12", "Torque teno virus 13", "Torque teno virus 16", "Torque teno virus 18", "Torque teno virus 19", "Torque teno virus 20", "Torque teno virus 22", "Torque teno virus 24", "Torque teno virus 27", "Torque teno virus 29",
                    "Betapapillomavirus 1", "Betapapillomavirus 2", "Betapapillomavirus 3", "Alphapolyomavirus octihominis", "Alphapolyomavirus quintihominis")){attention_level = "背景微生物"}
    }else {
       attention_level = ""
    }

    return(attention_level)
}

parasitus_attention <- function(species){
    known_species = c("Acanthamoeba castellanii", "Encephalitozoon hellem", "Toxoplasma gondii")

    attention_level = ""

    if(species %in% known_species){
        if(species %in% c("Acanthamoeba castellanii", "Encephalitozoon hellem")){attention_level = "结合批次检出情况确认"}
        if(species %in% c("Toxoplasma gondii")){attention_level = "结合批次检出情况确认"}
    }
}

extracted_result$Species <- as.character(extracted_result$Species)
extracted_result$Domain <- as.character(extracted_result$Domain)
extracted_result$Genius <- as.character(extracted_result$Genius)
extracted_result$label <- as.character(extracted_result$label)
extracted_result$speciescn <- as.character(extracted_result$speciescn)
# str(extracted_result)

extracted_result$attention_level <- apply(extracted_result, 1,function(x){
    species = x[1]
    rank = x[2]
    domain = x[4]
    genius = x[5]
    reads_num = x[6]
    aboundance_based_species = x[10]
    reads_in_NC = x[11]
    label = x[12]

    attention_level = ""

    # 根据 domain 的值调用不同的函数
    if (domain == "Bacteria") {
        attention_level = bacteria_attention(species, rank, genius, aboundance_based_species, reads_in_NC, extracted_result)
    } else if (domain == "Eukaryota") {
        attention_level = fungi_attention(species, rank, genius)
    } else if (domain == "Viruses") {
        attention_level = viruses_attention(species, rank, label)
    } else if (domain == "Parasitus") {
        attention_level = parasitus_attention(species)
    } else {
        # 可以选择在没有匹配的 domain 时输出警告或其他处理
        warning(paste("Unknown domain:", domain, "for species:", species))
    }
    
    return(attention_level)
})
# View(extracted_result)
extracted_result$percent <- paste0(extracted_result$percent, "%")
extracted_result$aboundance_based_domain <- paste0(extracted_result$aboundance_based_domain, "%")
extracted_result$aboundance_based_species <- paste0(extracted_result$aboundance_based_species, "%")

extracted_result <- extracted_result %>% select(attention_level, speciescn, tax_id, all_reads, CPM, Rank, Domain, Genius, Species, percent, aboundance_based_domain, aboundance_based_species, Reads_num_in_NC, label)
colnames(extracted_result)[14] <- "label_experience"

str(extracted_result)
write.table(extracted_result, paste0(output_dir, "/", sample_id, ".pathogen_report.txt"), row.names = F, sep = "\t", quote = F, fileEncoding = "UTF-8")
write.table(all_pathogen_results, paste0(output_dir, "/", sample_id, ".all_pathogen_report.txt"), row.names = F, sep = "\t", quote = F, fileEncoding = "UTF-8")