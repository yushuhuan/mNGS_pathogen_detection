
sample_id <- ""
seqtype <- ""

if(seqtype == "RNA"){
    sample_id=paste0(sample_id, "-R")
}else {
   if(seqtype == "cfDNA"){
        sample_id <- paste0(sample_id, "-CF")
   }
}

## 设定输出文件路径
output_dir <- paste0("/opt/ossfs/mNGS_pipeline_output/", sample_id, "/", "host_expression_info")