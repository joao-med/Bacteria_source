library(tidyverse)

bacteria_id <- analysis_data$ID %>% as.tibble()

#primeiramente vou gerar uma tibble que tenha info de conter ou nao um gene e nao
#a quantidade desse gene 

big_gene_tib <- tibble() #generating empty tibble

#Loop for that, based on bacteria ID, parse genbank files and then select genes and 
#put them on a tibble. To parse the genome genbankr from bioconductor is used 

checking_true <- 0
checking_false <- 0

for ( i in 1:length(bacteria_id)) {
  
  print(paste0( i ,' accessing genome ', bacteria_id[i]))
  
  gene_accession <- 
    genbankr::GBAccession(bacteria_id[i])
  
  parsing_genome <- 
    genbankr::readGenBank(gene_accession) #slower line to be processed 
  
  getting_genes  <- 
    parsing_genome %>%
    genbankr::genes()
  
  tib_genome <- #function that produce a tibble with gene infomration and targets informations 
    function(genome) { #that was found using rentrez package
      data_2bind <-  analysis_data %>% 
        filter(
          ID == bacteria_id[i]) %>%
        select(
          sub_strain, 
          country, 
          isolation_source, 
          host)
      
      tibble( 
        data_2bind,
        bacteria_ID = bacteria_id[i], 
        type = genome$type,
        genes_composition = genome$gene,
        locus_tag =  genome$locus_tag,
        loctype = genome$loctype,
        gene_ID = genome$gene_id
      )
    }
  
  print(paste0("seeing if there's gene annotation in ",  bacteria_id[i]))
  
  #some genomes does not have annotation, therefore, there is two options, the first
  #would be save only those that have it. the other would be make the annotations based
  #on database. The last would take a really long time studying and processing the data 
  
  if (length(getting_genes) >= 1 ) { 
    
    print(paste0(bacteria_id[i], " checked"))
    
    gene_data <-  
      tib_genome(genome = getting_genes ) %>% 
      select(bacteria_ID, country, host, isolation_source , genes_composition) %>% 
      mutate(rn = data.table::rowid(genes_composition)) %>% drop_na( genes_composition )
    
    
    gene_data <- gene_data[!duplicated((gene_data$genes_composition)),]
    
    checking_true <- checking_true + 1
    
    big_gene_tib <- bind_rows(big_gene_tib, gene_data)
    
  } else { 
    print(paste0(bacteria_id[i] ," ", length(getting_genes) >= 1))
    checking_false <- checking_false + 1
  }
} 

#check if data was correctly prospected is important. Then I sould work to find out how to 
#do that 

count_TF <- tibble(checking_false, checking_true)

expected_T <-
  big_gene_tib$bacteria_ID %>%
  unique() %>%
  length()

expected_T == count_TF$checking_true

#saving brute data

write.csv(big_gene_tib, "big_gene_tib.csv", row.names = F)

#generating data set of targets and gene presence
#be aware that the rows must not have names, ie. identification

big_gene_tib <- read.csv("big_gene_tib1500.csv")

gene_n_target_ds <- 
  big_gene_tib[,3:8] %>% 
  pivot_wider(
    names_from = genes_composition, 
    values_from = rn ) %>% 
  mutate(
    across(
      everything(),
      ~replace_na(.x,0)))

gene_n_target_ds %>% write.csv("gene_and_targets.csv")



entendendo <-  big_gene_tib %>% group_by(bacteria_ID) %>% count() %>% filter(n >=2000)

big_gene_tib[,1] %>% length()


big_gene_tib_2000 <- big_gene_tib %>% group_by(bacteria_ID) %>% filter( bacteria_ID %in% entendendo$bacteria_ID)


gene_n_target_ds <- 
  big_gene_tib_2000 %>% 
  pivot_wider(
    names_from = bacteria_ID, 
    values_from = bacteria_ID ) %>% 
  mutate(
    across(
      everything(), 
      ~replace_na(.x,0)))

big_gene_tib %>% group_by(bacteria_ID) %>% count() %>% arrange(desc(n))
adskad <- big_gene_tib %>% group_by(bacteria_ID) %>% count() %>% arrange(desc(n))
adskad$n %>% as.numeric() %>% median()
mean(n)
################

Counting_genomes %>% ggplot(aes(x = n, fill= x)) +
  geom_area()
Counting_genomes %>% ggplot(aes(x = n)) +
  geom_area(stat="bin")+
  theme_classic()


















