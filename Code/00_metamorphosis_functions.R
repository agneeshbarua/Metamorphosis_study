#Coloursets
cols.grouper=colorRampPalette(c("#25151b", "white", "#D9C343"))
cols.manini = colorRampPalette(c("#25151b", "white","#9acd32"))
cols.goby = colorRampPalette(c("#25151b","white","#BA4566"))
cols.clown = colorRampPalette(c("#25151b", "white", "#e77515"))
cols.danio = colorRampPalette(c("#25151b", "white", "#343BC2"))


#Get expression data and PCA
expression_matrix_and_pca<-function(sp,colramp){
  library(netZooR)
  
  samples<-read_csv(paste0("./01.Kallisto/",sp,"/samples.txt"),col_names = c("samples","stage"))
  files<-file.path(paste0("./01.Kallisto/",sp),samples$samples,"abundance.tsv")
  txi.kallisto<- tximport(files, type = "kallisto",txOut = T)
  cts<-txi.kallisto$counts
  groups<-factor(samples$stage)
  dge <- DGEList(cts,
             group = groups,
             genes = rownames(cts))
  keep <- rowSums(edgeR::cpm(dge) > 0.5) >= 3 
  dge<-dge[keep, ,keep.lib.sizes =F]
  dge<-calcNormFactors(dge,method = "TMM")
  dge_keep<-dge #save for later
  rownames(dge_keep$counts)<-paste0(sp,"_",rownames(dge_keep$counts))
  dge_pca<-vst(round(dge$counts),blind = T,fitType = "mean") #used only for plotting PCA
  #makePCA
  d<-SummarizedExperiment(dge_pca - rowMeans(dge_pca),colData = groups)
  pcaDat<-plotPCA(DESeqTransform(d),intgroup = "X", ntop = 500, returnData=T)
  percentVar_all<-round(100 * attr(pcaDat,"percentVar"))
  colramp = colramp(length(levels(pcaDat$group)))

  plot_pca<-ggplot(pcaDat,aes(PC1,PC2,label = name))+
    geom_point(aes(fill = group), color = "grey50", shape = 21, size = 5, alpha = 0.8)+
    scale_fill_manual(values = colramp)+
    theme_light()+
    xlab(paste0("PC1 (",percentVar_all[1],"%)")) +
    ylab(paste0("PC2 (",percentVar_all[2],"%)"))
  
  dat <- edgeR::cpm(dge, log = T,prior.count = 1e-5,normalized.lib.sizes = T)
  design <- model.matrix(~1+group, dge$samples)
  cobra_estimates <- netZooR::cobra(design, dat)
  rownames(cobra_estimates$Q) <- paste0(sp,"_",rownames(dat))
  colnames(cobra_estimates$Q) <- colnames(dat)
  Sigma_D<-t(cobra_estimates$Q)#Unlike the publication, we do not need to extract any batch specific effect. Our objective is to obtain a network of the entire transition, not stage specific ones. So we use the decomposed matrix as is.
  datExp<-t(dat)
  hclust(dist(datExp), method = "average")->clust
    return(list(dge_keep=dge_keep, pcaDat=pcaDat, plot=plot_pca, clust=clust, wgcna_in=Sigma_D))
}
expression_matrix_and_pca_clown<-function(){
  clown_ensembl<-readRDS("./01.Kallisto/clown_res/clown_biomart.rds")
  tx2gene<-data.frame(transcript=clown_ensembl$ensembl_transcript_id_version,
                      gene = clown_ensembl$ensembl_gene_id)
  samples<-read_tsv("./01.Kallisto/clown_res/samples.txt",col_names = "samples")
  files<-file.path("./01.Kallisto/clown_res",samples$samples,"abundance.tsv")
  txi.kallisto<- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  cts<-txi.kallisto$counts
  
  groups<-factor(paste0(rep(c("stage1","stage2","stage3","stage4","stage5","stage6","stage7"),each = 3)))

# Creating a DGEList object for use in edgeR.
  dge <- DGEList(cts,
                 group = groups,
                 genes = rownames(cts))
  keep <- rowSums(edgeR::cpm(dge) > 0.5) >= 3
  dge<-dge[keep, ,keep.lib.sizes =F]
  dge<-calcNormFactors(dge,method = "TMM")
  dge_keep<-dge
  dge_pca<-vst(round(dge$counts),blind = T,fitType = "mean") #used only for plotting PCA
  #dge<-calcNormFactors(dge,method = "TMM")

  d<-SummarizedExperiment(dge_pca - rowMeans(dge_pca),colData = groups)
  pcaDat<-plotPCA(DESeqTransform(d),intgroup = "X", ntop = 500, returnData=T)
  percentVar_all<-round(100 * attr(pcaDat,"percentVar"))
  colramp = colorRampPalette(c("#25151b","#a9552b","#e77515","#bd846b"))(7)
  
  plot_pca<-ggplot(pcaDat,aes(PC1,PC2,label = name))+
    geom_point(aes(fill = group), color = "grey50", shape = 21, size = 5, alpha = 0.8)+
    scale_fill_manual(values = colramp)+
    theme_light()+
    xlab(paste0("PC1 (",percentVar_all[1],"%)")) +
    ylab(paste0("PC2 (",percentVar_all[2],"%)"))
  
  dat <- edgeR::cpm(dge, log = T,prior.count = 1e-5,normalized.lib.sizes = T)
  design <- model.matrix(~1+group, dge$samples)
  cobra_estimates <- netZooR::cobra(design, dat)
  rownames(cobra_estimates$Q) <- rownames(dat)
  colnames(cobra_estimates$Q) <- colnames(dat)
  Sigma_D<-t(cobra_estimates$Q)#Unlike the publication, we do not need to extract any batch specificeffect. Our objective is to obtain a network of the entire transition, not stage specific ones. So we use decomposed matrix as is.
  datExp<-t(dat)
  hclust(dist(datExp), method = "average")->clust

  return(list(dge_keep=dge_keep, pcaDat=pcaDat, plot=plot_pca, clust=clust, wgcna_in=Sigma_D))
}

expression_matrix_and_pca_danio<-function(){
  danio_ensembl<-read_csv("./01.Kallisto/Danio_biomart.csv")
  tx2gene<-data.frame(transcript=danio_ensembl$`Transcript stable ID version`,
                      gene = danio_ensembl$`Gene stable ID`)
  samples<-read_csv("./01.Kallisto/Danio/samples.txt",col_names = c("samples","stage"))
  groups<-factor(samples$stage)
  files<-file.path("./01.Kallisto/Danio",samples$samples,"abundance.tsv")
  txi.kallisto<- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  cts<-txi.kallisto$counts
  
  #groups<-factor(paste0(rep(c("stage01","stage02","stage03","stage04","stage05","stage06","stage07","stage08","stage09","stage10","stage11"),each = 3)))

# Creating a DGEList object for use in edgeR.
  dge <- DGEList(cts,
                 group = groups,
                 genes = rownames(cts))
  keep <- rowSums(edgeR::cpm(dge) > 0.5) >= 3
  dge<-dge[keep, ,keep.lib.sizes =F]
  dge<-calcNormFactors(dge,method = "TMM")
  dge_keep<-dge
  dge_pca<-vst(round(dge$counts),blind = T,fitType = "mean") #used only for plotting PCA
  #dge<-calcNormFactors(dge,method = "TMM")

  d<-SummarizedExperiment(dge_pca - rowMeans(dge_pca),colData = groups)
  pcaDat<-plotPCA(DESeqTransform(d),intgroup = "X", ntop = 500, returnData=T)
  percentVar_all<-round(100 * attr(pcaDat,"percentVar"))
  colramp = colorRampPalette(c("#25151b","#2968D6","#29BFD6","#82D1E6"))(length(levels(groups)))
  
  plot_pca<-ggplot(pcaDat,aes(PC1,PC2,label = name))+
    geom_point(aes(fill = group), color = "grey50", shape = 21, size = 5, alpha = 0.8)+
    scale_fill_manual(values = colramp)+
    theme_light()+
    xlab(paste0("PC1 (",percentVar_all[1],"%)")) +
    ylab(paste0("PC2 (",percentVar_all[2],"%)"))
  
  dat <- edgeR::cpm(dge, log = T,prior.count = 1e-5,normalized.lib.sizes = T)
  #dat <- dge_vsd
  design <- model.matrix(~1+group, dge$samples)
  cobra_estimates <- netZooR::cobra(design, dat)
  rownames(cobra_estimates$Q) <- rownames(dat)
  colnames(cobra_estimates$Q) <- colnames(dat)
  Sigma_D<-t(cobra_estimates$Q)#Unlike the publication, we do not need to extract any batch specificeffect. Our objective is to obtain a network of the entire transition, not stage specific ones. So we use decomposed matrix as is.
  datExp<-t(dat)
  hclust(dist(datExp), method = "average")->clust

  return(list(dge_keep=dge_keep, pcaDat=pcaDat, plot=plot_pca, clust=clust, wgcna_in=Sigma_D))
}

## Define soft threshold
soft_threshold <- function(data) {
  library(WGCNA)
  allowWGCNAThreads(nThreads=10)
  # Choose a set of soft-thresholding powers, given a WGCNA data object
  powers = c(seq(from = 5, to=20, by=1))
  # Call the network topology analysis function
  sft = pickSoftThreshold(data, powerVector = powers, verbose = 5,RsquaredCut = 0.80,moreNetworkConcepts = T,networkType = "signed")
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.80,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  par(mfrow = c(1,1));
  
  return(sft)
}


## define adjacency function
wgcna_adjacency <- function(datExpr, minModuleSize=50, MEDissThres = 0.25, deepSplit = 2, power) {
  library(flashClust)
  library(WGCNA)
  # takes WGNCA data expression matrix, a soft threshold
  # optionally minimum module size, module eigengene dissimilarity threshold, and splitting threthold
  # returns adjacency matrix, module eigengenes and a data frame of colors
  enableWGCNAThreads()

  #compute adjacency
  adjacency <- adjacency(datExpr, power = power,type="signed")
  TOM <- TOMsimilarity(adjacency,TOMType="signed")
  # Calculate TOM-based connectivity
  geneTree <- flashClust::flashClust (as.dist(1-TOM), method = "average")

  # Module identification using dynamic tree cut
  dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = 1-TOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 50,);
  table(dynamicMods)
  dynamicColors = labels2colors(dynamicMods)

  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  METree = flashClust(as.dist(1-cor(MEs)), method = "average");
  plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 0)

  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs

  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs;

  # Recalculate MEs with color labels
  invisible(MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes)

  print(table(moduleColors))
  moduleColors <- as.data.frame(moduleColors)
  rownames(moduleColors) <- colnames(datExpr)
  
  # Calculate TOM-based connectivity
  MEs = orderMEs(MEs0)
  TOMConnectivity <- rowSums(TOM)
  
  # Module-specific TOM connectivity
  moduleTOMConnectivity<-function(mod){
  tibble(gene=moduleColors %>% filter(moduleColors == mod) %>% rownames(), 
         degree=TOMConnectivity[which(moduleColors == mod)],
         module = mod) %>% return()
}

purrr::map(unique(moduleColors$moduleColors),moduleTOMConnectivity) %>% purrr::list_rbind()->connectivity


  return(list(adjacency=adjacency,MEs=MEs,moduleColors=moduleColors, dynamicMods=dynamicMods, connectivity=connectivity, METree=METree, geneTree=geneTree))
}

#Module heatmaps
module_heatmaps<-function(dat,net,cols){

  colnames(dat$counts)<-dat$samples$group
#heat_dat<-dat$counts %>% reshape2::melt() %>% 
heat_dat<-edgeR::cpm(dat$counts, log = T,prior.count = 1e-5,normalized.lib.sizes = T) %>% reshape2::melt() %>% 
  dplyr::rename(genes = Var1,stages = Var2, logCPM = value) %>% 
  summarise(mean_logCPM = mean(logCPM),.by = c(genes,stages)) %>% group_by(genes) %>% 
  mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
heat_dat %>% summary()->sum
#clip the data at the 95 quantile to aid in visualisation
heat_dat<-heat_dat %>% mutate(z.clipped = case_when(
  z.score > quantile(heat_dat$z.score,0.95) ~ quantile(heat_dat$z.score,0.95),
  z.score < -quantile(heat_dat$z.score,0.95) ~ -quantile(heat_dat$z.score,0.95),
  T ~ z.score
))
qnt<-quantile(heat_dat$z.score,0.95) #used for plotting
make_heat_dat<-function(mod){
  heat_dat %>% filter(genes %in% (net$connectivity %>% filter(module == mod) %>% pull(gene))) %>% 
    summarise(mod.z = mean(z.clipped),.by = stages) %>% 
    mutate(mod = mod) %>% return()
  }

heat_plot_dat<-purrr::map(unique(net$connectivity$module),make_heat_dat) %>% list_rbind()
#Order the modules
wide<-heat_plot_dat %>% pivot_wider(names_from = stages,values_from = mod.z) 
wide<-wide %>% mutate(diff=wide[[length(wide)]]-wide[[2]])
wide<-wide %>% mutate(phase = case_when(
    diff < quantile(wide$diff, 0.33) ~ "early",
    diff > quantile(wide$diff, 0.33) & diff < quantile(wide$diff, 0.66) ~ "mid",
    diff >= quantile(wide$diff, 0.33) ~ "late"
)) %>% arrange(factor(phase, levels=c("early","mid","late")))

heat_plot_dat<-full_join(wide,heat_plot_dat, by = "mod") %>% arrange(diff) %>% dplyr::select(mod,stages,mod.z,phase)
#plot
heat<-ggplot(data = heat_plot_dat,aes(x = stages,y = as.factor(mod)))+
  facet_grid(rows = vars(factor(phase, levels = c("early","mid","late"))),
                         scales = "free")+
  geom_tile(aes(fill = mod.z), color = "grey80",) +
  scale_fill_gradientn(colors = (cols(11)), limits = c(-qnt, qnt),
                       breaks = c(-qnt, 0, qnt), labels = c(paste0("<", round(-qnt)), "0", paste0(">", round(qnt))))+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  ylab("Modules")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 6),
        strip.background = element_rect(fill = "grey95", color = "black"))
  return(list(heat=heat, sum=sum, dat=heat_plot_dat))
}

#Module lineplots
module_pres_lineplots<-function(dat,mp_net,orth_exp,plot_mods){
  dat<-clown_orth_net_in$dge_keep
  orth_exp<-clown_orth_exp
  colnames(dat$counts)<-dat$samples$group
  heat_dat<-dat$counts %>% reshape2::melt() %>% 
  dplyr::rename(orths = Var1,stages = Var2, logCPM = value) %>% 
  summarise(mean_logCPM = mean(logCPM),.by = c(orths,stages)) %>% group_by(orths) %>% 
  mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup() 
  
  heat_dat<-heat_dat %>% mutate(z.clipped = case_when(
  z.score > quantile(heat_dat$z.score,0.95) ~ quantile(heat_dat$z.score,0.95),
  z.score < -quantile(heat_dat$z.score,0.95) ~ -quantile(heat_dat$z.score,0.95),
  T ~ z.score
))
  heat_dat<-left_join(heat_dat,orth_exp) %>% drop_na() %>%
    dplyr::select(-c(contains("Sample"),sp))
#heat_dat<-heat_dat %>% dplyr::select(-c(contains("Sample"),sp))
  
  make_line_plot<-function(mod){
    plot_dat<-heat_dat %>% filter(Orthogroup %in% (mp_net$connectivity %>% filter(module == mod)
                                                   %>% pull(gene)))
    qnt<-quantile(plot_dat$z.score,0.95) #used for plotting
    # Define where to place the vertical panel (on the right side)
    x_rect_start <- length(levels(plot_dat$stages))+ .1
    x_rect_end <- length(levels(plot_dat$stages))+ 1  # Adjust width of the rectangle as needed
    rec_fill = mod
    lp<-ggplot(plot_dat,aes(x = stages, y = z.clipped,group = orths))+
      geom_line(alpha = 0.3, color = "grey70")+ 
      coord_cartesian(ylim = c(-qnt, qnt)) +
      geom_rect(aes(xmin = x_rect_start, 
                    xmax = x_rect_end, 
                    ymin = -Inf, 
                    ymax = Inf), 
                fill = mod, alpha = 0.2,color = "white") +  # Add vertical rectangle on the right
      geom_line(data = (plot_dat %>% summarise(mean_exp = mean(z.score),.by = "stages")),
                aes(x = stages, y = mean_exp, group = 1),linewidth = 1, color = "grey30")+
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            axis.ticks.x = element_blank())
    return(lp)
  }
  colors<-as.list(plot_mods$mods)
  purrr::map(colors,make_line_plot)->line_plots
  
  return(line_plots %>% patchwork::wrap_plots(ncol = 5))
}


#Module preservation heatmaps
module_heatmaps_pres_analysis<-function(dat,net,cols){
  
heat_dat<-dat %>% reshape2::melt() %>% 
  dplyr::rename(genes = Orthogroup,stages = variable, logCPM = value) %>% 
  summarise(mean_logCPM = mean(logCPM),.by = c(orths,genes,module,stages)) %>% group_by(genes) %>% 
  mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
heat_dat %>% summary()->sum
#clip the data at the 95 quantile to aid in visualisation
heat_dat<-heat_dat %>% mutate(z.clipped = case_when(
  z.score > quantile(heat_dat$z.score,0.95) ~ quantile(heat_dat$z.score,0.95),
  z.score < -quantile(heat_dat$z.score,0.95) ~ -quantile(heat_dat$z.score,0.95),
  T ~ z.score
))
qnt<-quantile(heat_dat$z.score,0.95) #used for plotting
make_heat_dat<-function(mod){
  #When running on module pres data change the orths in the filter to genes
  heat_dat %>% filter(genes %in% (net$connectivity %>% filter(module == mod) %>% pull(gene))) %>% 
    summarise(mod.z = mean(z.clipped),.by = stages) %>% 
    mutate(mod = mod) %>% return()
}

heat_plot_dat<-purrr::map(unique(net$connectivity$module),make_heat_dat) %>% list_rbind()
heat_plot_dat<-heat_plot_dat %>% mutate(stages=as.character(stages), mod=as.character(mod), mod.z=as.numeric(mod.z)) %>% distinct(stages,mod.z,mod)

#Order the modules
wide<-heat_plot_dat %>% pivot_wider(names_from = stages,values_from = mod.z) 
heat_plot_dat<-full_join(heat_plot_dat %>% pivot_wider(names_from = stages,values_from = mod.z) %>% rowwise() %>% mutate(diff = list(wide[[ncol(wide)]]-wide[[2]]
)) %>% arrange(diff) %>% 
  dplyr::select(mod),heat_plot_dat,by = c("mod"))
heat_plot_dat<-heat_plot_dat %>% mutate(mod = factor(mod,levels = unique(mod)))

heat<-ggplot(data = heat_plot_dat,aes(y = stages,x = as.factor(mod)))+
  geom_tile(aes(fill = mod.z), color = "grey80",) +
  scale_fill_gradientn(colors = (cols(11)), limits = c(-qnt, qnt),
                       breaks = c(-qnt, 0, qnt), labels = c(paste0("<", round(-qnt)), "0", paste0(">", round(qnt))))+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(panel.background = element_blank())
return(list(plot=heat, dat=heat_plot_dat, sum=sum))
}

#Module preservations and expression plots
make_tiles_and_mp_plot<-function(d,cols){
  plot_dat<-d %>% mutate(pres = case_when(
    z.score <=6  ~ "no preservation",
    z.score > 6 & z.score <= 12 ~ "low to medium preservation",
    z.score > 12 ~ "high preservation")) %>% drop_na()
  plot_dat$mods<-fct_reorder(plot_dat$mods, plot_dat$rank)
  qnt<-quantile(plot_dat$mod.z,0.95)
  p<-ggplot(plot_dat, aes(x = mods, y = stages)) +
    facet_grid(.~ pres, scales = "free", space = "free") +
    geom_tile(aes(fill = mod.z),color="white") +
    scale_fill_gradientn(colors = (cols(11)), limits = c(-qnt, qnt+0.3),
                         breaks = c(-qnt, 0, qnt),
                         labels = c(paste0("<", round(-qnt)), "0", paste0(">",round(qnt))))+
    geom_tile(aes(color = pres),alpha=0.1,size=1)+
    scale_color_manual(values = c("#497BB6", "#B6497B", "#7BB649"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90))
  return(list(plot=p, dat=plot_dat))
}

#Summaries of orth network
check_orth_exp_sumarries<-function(a,b,exp_dat){
  exp_dat %>% filter(Orthogroup %in% c(a,b)) %>% distinct(orths,.keep_all = T)->raw_dat

  exp_dat %>% filter(Orthogroup %in% c(a,b)) %>% distinct(orths,.keep_all = T) %>% 
  reframe(across(starts_with("Sample"), \(x) mean(x, na.rm = TRUE)),.by = ("Orthogroup"))->og_mean
  
  exp_dat %>% filter(Orthogroup %in% c(a,b)) %>% distinct(orths,.keep_all = T) %>% 
  reframe(across(starts_with("Sample"), \(x) sum(x, na.rm = TRUE)),.by = ("Orthogroup"))->og_sum
  
  exp_dat %>% filter(Orthogroup %in% c(a,b)) %>% distinct(orths,.keep_all = T) %>% 
  reframe(across(starts_with("Sample"), \(x) median(x, na.rm = TRUE)),.by = ("Orthogroup"))->og_median
  
  exp_dat %>% filter(Orthogroup %in% c(a,b)) %>% distinct(orths,.keep_all = T) %>% 
  filter(orths %in% (exp_dat %>% filter(Orthogroup %in% c(a,b)) %>% 
           distinct(orths,.keep_all = T) %>% rowwise(orths) %>% 
           mutate(m=mean(Sample1:Sample21)) %>% ungroup() %>%  
           slice_max(m,by = "Orthogroup") %>% pull(orths)))->og_max
  
  exp_dat %>% filter(Orthogroup %in% c(a,b)) %>% distinct(orths,.keep_all = T) %>% 
  filter(orths %in% (exp_dat %>% filter(Orthogroup %in% c(a,b)) %>% 
                       distinct(orths,.keep_all = T) %>% rowwise(orths) %>% 
                       mutate(m=mean(Sample1:Sample21)) %>% ungroup() %>%  
                       slice_min(m,by = "Orthogroup") %>% pull(orths)))->og_min
  
  raw_dat$orths<-fct_reorder(raw_dat$orths, raw_dat$Orthogroup)
  
  p1<-ggplot(raw_dat %>% pivot_longer(cols = contains("Sample")) %>% 
               mutate(Orthogroup=paste0("OG-",Orthogroup)),
       aes(x=name, y=log10(value+1), group=orths))+
  facet_wrap(~orths)+
  geom_line(aes(colour=Orthogroup),size=0.7)+
  ggsci::scale_color_aaas()+
  ggtitle("Expression of individual paralogs")+
  scale_y_continuous(limits = c(0,6))+
  theme_light()+
  ylab("log10(cpm+1)")+
  theme(axis.text.x = element_text(angle = 90),size=7)
  
  p2<-ggplot(og_mean %>% pivot_longer(cols = contains("Sample")) %>% 
               mutate(Orthogroup=paste0("OG-",Orthogroup)),
       aes(x=name, y=log10(value+1), group=Orthogroup))+
  facet_wrap(~Orthogroup)+
  geom_line(aes(colour=Orthogroup),size=0.7)+
  ggsci::scale_color_aaas()+
  scale_y_continuous(limits = c(0,6))+
  ggtitle("Mean expression of paralogs")+
  theme_light()+
  ylab("log10(cpm+1)")+
  theme(axis.text.x = element_text(angle = 90,size = 7))+
  theme(legend.position = "none")
  
  p3<-ggplot(og_median %>% pivot_longer(cols = contains("Sample")) %>% 
               mutate(Orthogroup=paste0("OG-",Orthogroup)),
           aes(x=name, y=log10(value+1), group=Orthogroup))+
  facet_wrap(~Orthogroup)+
  geom_line(aes(colour=Orthogroup),size=0.7)+
  ggsci::scale_color_aaas()+
  scale_y_continuous(limits = c(0,6))+
  ggtitle("Median expression of paralogs")+
  theme_light()+
  ylab("log10(cpm+1)")+
  theme(axis.text.x = element_text(angle = 90,size = 7))+
  theme(legend.position = "none")
  
  p4<-ggplot(og_sum %>% pivot_longer(cols = contains("Sample")) %>% 
               mutate(Orthogroup=paste0("OG-",Orthogroup)),
           aes(x=name, y=log10(value+1), group=Orthogroup))+
  facet_wrap(~Orthogroup)+
  geom_line(aes(colour=Orthogroup),size=0.7)+
  ggsci::scale_color_aaas()+
  scale_y_continuous(limits = c(0,6))+
  ggtitle("Summed expression of paralogs")+
  theme_light()+
  ylab("log10(cpm+1)")+
  theme(axis.text.x = element_text(angle = 90,size = 7))+
  theme(legend.position = "none")
  
  p5<-ggplot(og_max %>% pivot_longer(cols = contains("Sample")) %>% 
               mutate(Orthogroup=paste0("OG-",Orthogroup)),
           aes(x=name, y=log10(value+1), group=Orthogroup))+
  facet_wrap(~Orthogroup)+
  geom_line(aes(colour=Orthogroup),size=0.7)+
  ggsci::scale_color_aaas()+
  scale_y_continuous(limits = c(0,6))+
  ggtitle("Highest expressed paralog")+
  theme_light()+
  ylab("log10(cpm+1)")+
  theme(axis.text.x = element_text(angle = 90,size = 7))+
  theme(legend.position = "none")
  
  p6<-ggplot(og_min %>% pivot_longer(cols = contains("Sample")) %>% 
               mutate(Orthogroup=paste0("OG-",Orthogroup)),
           aes(x=name, y=log10(value+1), group=Orthogroup))+
  facet_wrap(~Orthogroup)+
  geom_line(aes(colour=Orthogroup),size=0.7)+
  ggsci::scale_color_aaas()+
  scale_y_continuous(limits = c(0,6))+
  ggtitle("Lowest expressed paralog")+
  theme_light()+
  ylab("log10(cpm+1)")+
  theme(axis.text.x = element_text(angle = 90,size = 7))+
  theme(legend.position = "none")
  
  p1 | p2 + p3 + p4 + p5 + p6 + patchwork::plot_layout(widths = c(1,1),guides = "keep") %>% return()
}

#Orthologous networks
make_orth_net_input<-function(d){
  
  dge <- DGEList(d,group = colnames(d),genes = rownames(d))
  dge<-calcNormFactors(dge,method = "TMM")
  dge_keep<-dge
  dat <- edgeR::cpm(dge, log = T,prior.count = 1e-5,normalized.lib.sizes = T)
  design <- model.matrix(~1+group, dge$samples)
  cobra_estimates <- cobra(design, dat)
  rownames(cobra_estimates$Q) <- rownames(dat)
  colnames(cobra_estimates$Q) <- colnames(dat)
  Sigma_D<-t(cobra_estimates$Q)#Unlike the publication, we do not need to extract any batch specific effect. Our objective is to obtain a network of the entire transition, not stage specific ones. So we use decomposed matrix as is.
  datExp<-t(dat)
  hclust(dist(datExp), method = "average")->clust
  return(list(dge_keep=dge_keep, wgcna_in=Sigma_D, clust=clust))
}

#Module preservation plots
module_pres_plots<-function(mp,heatdat,cols,c1,c2){
mp<-mp
modPreservStats = mp$preservation$Z #Extract module preservation statistics
Z_summary = modPreservStats[[c1]][[c2]]$Zsummary.pres #Summarize the Z-summary for all modules
print(Z_summary)
d<-tibble(mods = rownames(modPreservStats[[c1]][[c2]]), 
                     z.score = modPreservStats[[c1]][[c2]]$Zsummary.pres,
                     rank=mp$preservation$observed[[c1]][[c2]]$medianRank.pres) %>% arrange(rank)
d$mods<-factor(d$mods,levels = d$mods)
left_join(heatdat %>% dplyr::rename(mods=mod),d,by = "mods") %>% arrange(rank) ->plot_dat
make_tiles_and_mp_plot(d = plot_dat,cols)->p

return(p)
}

#Jaccard index and Fishers exact test
get_jaccard_estimate<-function(col1,col2,net1,net2){
  jac<-function(net1,net2,col1,col2) {
    #Get the indices of the genes in each of the modules
    g1<-which((net1 %>% rownames_to_column("g") %>% arrange(g))$moduleColors==col1)
    g2<-which((net2 %>% rownames_to_column("g") %>% arrange(g))$moduleColors==col2)
    # Jaccard index
    #checks the similarity in terms of whether the same genes are present in the same module(irrespective of colour of the modules).
    #Basically checks if genes tend to co-occur together.
    jaccard_index <- length(intersect(g1, g2)) / length(union(g1, g2))
    
    #Fisher's exact test to see of there are more significant overlaps
    uni<-length(clown_orth_net$moduleColors$moduleColors) #all orthologous networks have the same length
    overlap_matrix <- matrix(c(length(intersect(g1, g2)), length(setdiff(g1, g2)),
                               length(setdiff(g2, g1)), uni - length(union(g1, g2))), 
                             nrow = 2)
    fisher.test(overlap_matrix,alternative = "greater")->res
    data.frame(jaccard_index=jaccard_index, c1=col1, c2=col2, f.test.pval=res$p.value) %>% return()
  }
  jac(net1 = net1, net2 = net2, col1 = col1, col2 = col2) -> d
  return(d)
}

#Getting homologous modules
 homo_mods<-function(c1,c2,c3,c4,c5){
  #Combing all the GO annotations into a single df and filter out the terms that are present in at least 4 species.
   #Please not that these are not specific annotations of the modules, but rather generalisations.
   rbind(clown_orth_mods_proc[[c1]][[1]],
         grouper_orth_mods_proc[[c2]][[1]],
         manini_orth_mods_proc[[c3]][[1]],
         goby_orth_mods_proc[[c4]][[1]],
         danio_orth_mods_proc[[c5]][[1]]) %>% 
     summarise(n=n(), .by = "GOBPID") %>% filter(n>=4) %>% 
     pull(GOBPID)->common_GOs

   
   #Use the common GO terms to filter out annotations and filter the FDR sig ones and pick the top 10 sig
   #These top 10 terms should largely be relevant across modules since we picked common terms in at least 4 species
   rbind(clown_orth_mods_proc[[c1]][[1]] %>% mutate(sp="clown"),
         grouper_orth_mods_proc[[c2]][[1]] %>% mutate(sp="grouper"),
         manini_orth_mods_proc[[c3]][[1]] %>% mutate(sp="manini"),
         goby_orth_mods_proc[[c4]][[1]] %>% mutate(sp="goby"),
         danio_orth_mods_proc[[c5]][[1]] %>% mutate(sp="danio")) %>% filter(GOBPID %in% common_GOs) %>% 
     filter(FDR < 0.05) %>% mutate(col.clown=c1) %>% arrange(FDR) %>% head(30) %>% return()
 }
 
# Homologous modules line plots
  homo_mods_lines<-function(mod.clown,mod.manini,mod.goby,mod.grouper,mod.danio,GO_dat){
    len=10
#mod.clown=mods_combs$col.clown[1]
#mod.manini=mods_combs$col.manini[1]
#mod.goby=mods_combs$col.goby[1]
#mod.grouper=mods_combs$col.grouper[1]
#mod.danio=mods_combs$col.danio[1]

    d.clown<-clown_orth_net_in$dge_keep$counts[(clown_orth_net$moduleColors)$moduleColors==mod.clown,] %>% edgeR::cpm(log = T,prior.count = 1e-5,normalized.lib.sizes = T) %>% 
      reshape2::melt() %>% 
      dplyr::rename(orths = Var1,stages = Var2, logCPM = value) %>% 
      summarise(mean_logCPM = mean(logCPM),.by = c(orths,stages)) %>% group_by(orths) %>% 
      mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
    d.clown<-d.clown %>% mutate(z.clipped = case_when(
      z.score > quantile(d.clown$z.score,0.95) ~ quantile(d.clown$z.score,0.95),
      z.score < -quantile(d.clown$z.score,0.95) ~ -quantile(d.clown$z.score,0.95),
      T ~ z.score
    )) %>% mutate(sp="clown")
    
    d.manini<-manini_orth_net_in$dge_keep$counts[(manini_orth_net$moduleColors)$moduleColors==mod.manini,] %>% edgeR::cpm(log = T,prior.count = 1e-5,normalized.lib.sizes = T) %>% 
      reshape2::melt() %>% 
      dplyr::rename(orths = Var1,stages = Var2, logCPM = value) %>% 
      summarise(mean_logCPM = mean(logCPM),.by = c(orths,stages)) %>% group_by(orths) %>% 
      mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
    d.manini<-d.manini %>% mutate(z.clipped = case_when(
      z.score > quantile(d.manini$z.score,0.95) ~ quantile(d.manini$z.score,0.95),
      z.score < -quantile(d.manini$z.score,0.95) ~ -quantile(d.manini$z.score,0.95),
      T ~ z.score
    )) %>% mutate(sp="manini")
    
    d.goby<-goby_orth_net_in$dge_keep$counts[(goby_orth_net$moduleColors)$moduleColors==mod.goby,] %>% edgeR::cpm(log = T,prior.count = 1e-5,normalized.lib.sizes = T) %>%
      reshape2::melt() %>% 
      dplyr::rename(orths = Var1,stages = Var2, logCPM = value) %>% 
      summarise(mean_logCPM = mean(logCPM),.by = c(orths,stages)) %>% group_by(orths) %>% 
      mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
    d.goby<-d.goby %>% mutate(z.clipped = case_when(
      z.score > quantile(d.goby$z.score,0.95) ~ quantile(d.goby$z.score,0.95),
      z.score < -quantile(d.goby$z.score,0.95) ~ -quantile(d.goby$z.score,0.95),
      T ~ z.score
    )) %>% mutate(sp="goby")
    
    d.grouper<-grouper_orth_net_in$dge_keep$counts[(grouper_orth_net$moduleColors)$moduleColors==mod.grouper,] %>% edgeR::cpm(log = T,prior.count = 1e-5,normalized.lib.sizes = T) %>% 
      reshape2::melt() %>% 
      dplyr::rename(orths = Var1,stages = Var2, logCPM = value) %>% 
      summarise(mean_logCPM = mean(logCPM),.by = c(orths,stages)) %>% group_by(orths) %>% 
      mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
    d.grouper<-d.grouper %>% mutate(z.clipped = case_when(
      z.score > quantile(d.grouper$z.score,0.95) ~ quantile(d.grouper$z.score,0.95),
      z.score < -quantile(d.grouper$z.score,0.95) ~ -quantile(d.grouper$z.score,0.95),
      T ~ z.score
    )) %>% mutate(sp="grouper")
    
    d.danio<-danio_orth_net_in$dge_keep$counts[(danio_orth_net$moduleColors)$moduleColors==mod.danio,] %>% edgeR::cpm(log = T,prior.count = 1e-5,normalized.lib.sizes = T) %>% 
      reshape2::melt() %>% 
      dplyr::rename(orths = Var1,stages = Var2, logCPM = value) %>% 
      summarise(mean_logCPM = mean(logCPM),.by = c(orths,stages)) %>% group_by(orths) %>% 
      mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
    d.danio<-d.danio %>% mutate(z.clipped = case_when(
      z.score > quantile(d.danio$z.score,0.95) ~ quantile(d.danio$z.score,0.95),
      z.score < -quantile(d.danio$z.score,0.95) ~ -quantile(d.danio$z.score,0.95),
      T ~ z.score
    )) %>% mutate(sp="danio")
    
    
    d.all<-rbind(d.clown, d.grouper, d.goby, d.manini, d.danio)
    
    
    d.all %>% filter(sp=="clown") %>%  
      ggplot(aes(x=stages, y=z.clipped))+
      ggbeeswarm::geom_quasirandom(shape = 21, color = "grey90", alpha = 0.8, size = 1.5,
                                   fill = "#e77515")+
      geom_boxplot(outliers = FALSE, fill = "grey90",alpha=0.5)+
      geom_line(data = (d.all %>% filter(sp=="clown")
                        %>% summarise(mean_exp = mean(z.score),.by = "stages")),
                aes(x = stages, y = mean_exp, group = 1),linewidth = 1, color = "black")+
      annotate(geom = "rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
               fill = "#594BE4", alpha = 0.2) +  
      ylim(-2.2,2.2)+
      ylab("Z-score")+
      theme_linedraw()+
    theme(
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank())->p1
    
    d.all %>% filter(sp=="manini") %>%  
      ggplot(aes(x=stages, y=z.clipped))+
      ggbeeswarm::geom_quasirandom(shape = 21, color = "grey90", alpha = 0.8, size = 1.5,
                                   fill = "#9acd32")+
            geom_boxplot(outliers = FALSE, fill = "grey90",alpha=0.5)+

      geom_line(data = (d.all %>% filter(sp=="manini")
                        %>% summarise(mean_exp = mean(z.score),.by = "stages")),
                aes(x = stages, y = mean_exp, group = 1),linewidth = 1, color = "black") + 
      ylim(-2.2,2.2)+
      ylab("Z-score")+
      theme_linedraw()+
      theme(
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank())->p2
    
    d.all %>% filter(sp=="goby") %>%  
      ggplot(aes(x=stages, y=z.clipped))+
      ggbeeswarm::geom_quasirandom(shape = 21, color = "grey90", alpha = 0.8, size = 1.5,
                                   fill = "#BA4566")+
            geom_boxplot(outliers = FALSE, fill = "grey90",alpha=0.5)+
      geom_line(data = (d.all %>% filter(sp=="goby")
                        %>% summarise(mean_exp = mean(z.score),.by = "stages")),
                aes(x = stages, y = mean_exp, group = 1),linewidth = 1, color = "black") + 
      ylim(-2.2,2.2)+
      ylab("Z-score")+
      theme_linedraw()+
     theme(
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank())->p3
    
    d.all %>% filter(sp=="grouper") %>%  
      ggplot(aes(x=stages, y=z.clipped))+
      ggbeeswarm::geom_quasirandom(shape = 21, color = "grey90", alpha = 0.8, size = 1.5,
                                   fill = "#D9C343")+
            geom_boxplot(outliers = FALSE, fill = "grey90",alpha=0.5)+
      geom_line(data = (d.all %>% filter(sp=="grouper")
                        %>% summarise(mean_exp = mean(z.score),.by = "stages")),
                aes(x = stages, y = mean_exp, group = 1),linewidth = 1, color = "black")+
      annotate(geom = "rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
               fill = "#594BE4", alpha = 0.15) + 
      annotate(geom = "rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
               fill = "#594BE4", alpha = 0.15) +
      ylim(-2.2,2.2)+
      ylab("Z-score")+
      theme_linedraw()+
      theme(
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank())->p4
    
    d.all %>% filter(sp=="danio") %>%  
      ggplot(aes(x=stages, y=z.clipped))+
      ggbeeswarm::geom_quasirandom(shape = 21, color = "grey90", alpha = 0.8, size = 1.5,
                                   fill = "#594BE4")+
            geom_boxplot(outliers = FALSE, fill = "grey90",alpha=0.5)+
      geom_line(data = (d.all %>% filter(sp=="danio")
                        %>% summarise(mean_exp = mean(z.score),.by = "stages")),
                aes(x = stages, y = mean_exp, group = 1),linewidth = 1, color = "black")+
      ylim(-2.2,2.2)+
      ylab("Z-score")+
      theme_linedraw()+
      theme(
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank())->p5
    
    #GO terms
    #Use different go results when needed
    GO_dat %>% filter(mod==mod.clown) %>% slice_head(n=10) %>% 
      ggplot(aes(label=Term, y=c(1:10)))+
      geom_text(x=1,hjust = "left",size = 3)+
      xlim(1,len)+
      ylab(NULL)+
      labs(subtitle = "Enriched GO terms" )+
      theme_minimal()+
     theme(
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank())->p6
    
    wrap_plots(p1,p2,p3,p4,p5,p6,
               nrow = 1, 
               ncol = 6, 
               widths = c(0.5,0.5,0.5,0.5,0.5,1),
               heights = c(1,1,1,1,1,0.05)) +
      theme(plot.background = element_rect(fill='transparent'),
            legend.background = element_rect(fill = 'transparent'))+
  plot_annotation(title = paste0(mod.clown))->wp
    wp %>% return()
    #change the destination per usage
    #ggsave(plot = wp, filename = paste0("./Figures/FigS5/",mod.clown,".png"), device = "png", width = 12, height = 2, units = "in",bg = "transparent")
  }


correlation_test <- function(x, y) {
  test <- cor.test(x, y, method = "spearman")
  tibble(
    corr = test$estimate,
    p_value_adj = p.adjust(test$p.value,method = "BH",)
  ) 
}

#Used in the validation section
three_system_col_plots<-function(orth_net, col){
  
  orth_net$connectivity %>% dplyr::select(gene,module) %>% 
    filter(gene %in% full_vision_orths$Orthogroup) %>% 
    mutate(sys="vision system") %>% dplyr::count(sys, module) %>% 
    mutate(per_of_gene_set=(n/sum(n))*100) %>% 
    rbind(
      orth_net$connectivity %>% dplyr::select(gene,module) %>% 
        filter(gene %in% full_thyroid_orths$Orthogroup) %>% 
        mutate(sys="thyroid hormone regulation") %>% dplyr::count(sys, module) %>% 
        mutate(per_of_gene_set=(n/sum(n))*100)
    ) %>% 
    rbind(
      orth_net$connectivity %>% dplyr::select(gene,module) %>% 
        filter(gene %in% full_lipid_orths$Orthogroup) %>% 
        mutate(sys="lipid metabolism") %>% dplyr::count(sys, module) %>% 
        mutate(per_of_gene_set=(n/sum(n))*100)
    ) %>% 
    mutate(coexp=case_when(n >= 2 ~"yes", .default = "no")) %>% 
    mutate(module=reorder(module, -per_of_gene_set)) %>% 
    ggplot(aes(y=per_of_gene_set, x=module, fill=coexp))+
    facet_wrap(~sys, scales="free", nrow=3)+
    geom_col(color="black", linewidth=0.2, alpha=0.9)+
    geom_text(aes(y=per_of_gene_set-1, label = paste0(round(per_of_gene_set,1),"%")), 
              size=2, color="grey95")+
    scale_fill_manual(values = c("#25151b",col), name="contains multiple\n system genes")+
    labs(y="Percentage of gene set")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90,size=8,hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size=8),
          strip.background = element_rect(fill = "grey95"))
}

#Makes plots for module eigengene correlations, line plots of correlated modules, and GO plots
ME_GO_plots<-function(orth_net_MEs,
                      orth_net_in,
                      orth_net,
                      orth_mods_proc,
                      c1,c2,
                      col,
                      xmin,
                      xmax){

names(orth_net_MEs)<-str_remove(names(orth_net_MEs), "ME")
  
  orth_net_MEs %>% dplyr::select(any_of(c(c1, c2))) %>% 
  ggpubr::ggscatter( x = c1, y = c2,
          color = "grey80", shape = 21, size = 3,
          fill=col,# Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 0, label.sep = "\n"))+
    xlab(paste0(c1," module eigengene"))+
    ylab(paste0(c2," module eigengene"))+
  theme_linedraw()->p1

d.goby<-orth_net_in$dge_keep$counts[(orth_net$moduleColors)$moduleColors==c1,] %>% reshape2::melt() %>% 
  dplyr::rename("orths" = "Var1","stages" = "Var2", "logCPM" = "value") %>% 
  summarise(mean_logCPM = mean(logCPM),.by = c(orths,stages)) %>% group_by(orths) %>% 
  mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
d.goby<-d.goby %>% mutate(z.clipped = case_when(
  z.score > quantile(d.goby$z.score,0.95) ~ quantile(d.goby$z.score,0.95),
  z.score < -quantile(d.goby$z.score,0.95) ~ -quantile(d.goby$z.score,0.95),
  T ~ z.score
)) %>% mutate(sp="goby")

d.goby %>% filter(sp=="goby") %>%  
  ggplot(aes(x=stages, y=z.clipped, group=orths))+
  ggbeeswarm::geom_quasirandom(shape = 21, color = "grey90", alpha = 0.8, size = 1.5,
                               fill = col)+
  geom_line(data = (d.goby %>% filter(sp=="goby")
                    %>% summarise(mean_exp = mean(z.score),.by = "stages")),
            aes(x = stages, y = mean_exp, group = 1),linewidth = 1, color = "black")+
  annotate(geom = "rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
           fill = "#594BE4", alpha = 0.15) + 
  ylim(-2.2,2.2)+
  ylab("Z-score")+
  theme_linedraw()+
  theme(axis.text.x = element_blank())->p2


d.goby<-orth_net_in$dge_keep$counts[(orth_net$moduleColors)$moduleColors==c2,] %>% reshape2::melt() %>% 
  dplyr::rename("orths" = "Var1","stages" = "Var2", "logCPM" = "value") %>% 
  summarise(mean_logCPM = mean(logCPM),.by = c(orths,stages)) %>% group_by(orths) %>% 
  mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>% ungroup()
d.goby<-d.goby %>% mutate(z.clipped = case_when(
  z.score > quantile(d.goby$z.score,0.95) ~ quantile(d.goby$z.score,0.95),
  z.score < -quantile(d.goby$z.score,0.95) ~ -quantile(d.goby$z.score,0.95),
  T ~ z.score
)) %>% mutate(sp="goby")

d.goby %>% filter(sp=="goby") %>%  
  ggplot(aes(x=stages, y=z.clipped, group=orths))+
  ggbeeswarm::geom_quasirandom(shape = 21, color = "grey90", alpha = 0.8, size = 1.5,
                               fill = col)+
  geom_line(data = (d.goby %>% filter(sp=="goby")
                    %>% summarise(mean_exp = mean(z.score),.by = "stages")),
            aes(x = stages, y = mean_exp, group = 1),linewidth = 1, color = "black")+
  annotate(geom = "rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
           fill = "#594BE4", alpha = 0.15) + 
  ylim(-2.2,2.2)+
  ylab("Z-score")+
  theme_linedraw()+
  theme(axis.text.x = element_blank())->p3

rbind(orth_mods_proc[[c1]]$result %>% head(10) %>%  mutate(mod=paste0(c1)),
orth_mods_proc[[c2]]$result %>% head(10) %>% mutate(mod=paste0(c2))) %>% mutate(Term=reorder(Term,Count)) %>% 
  ggplot(aes(x=Term, y=-log(FDR), size = log10(Count)))+
  facet_wrap(~mod, nrow = 2,scales = "free")+
  geom_segment( aes(xend=Term, yend=0),size=0.4)+
  geom_point(fill=col,color="grey90",shape=21) +
  scale_x_discrete(position = "top")+
  ggtitle("Top 10 enriched GO terms")+
  coord_flip()+
  theme_light()+
  theme(plot.title = element_text(size = 10))->p4


wrap_plots(wrap_plots(p3,p1,p2, nrow = 3,ncol = 1, heights = c(0.5,1,0.5)), p4,guides = "collect") %>% 
  return()
}

#Make the expression matrix data frame
no_of_tissue_samples_in_each_fish<-function(x){
  bgee_dat %>% dplyr::filter(Anatomical.entity.name == x) %>% dplyr::group_by(Species) %>% dplyr::summarise(Libraries = length(unique(Library.ID)))
}
make_exp_data_frames<-function(.data, ...){
  .data  %>% 
    dplyr::select(Gene.ID, TPM, Library.ID, Detection.flag) %>% 
    filter(Detection.flag == "present") %>% dplyr::select(-c(Detection.flag)) %>% 
    pivot_wider(names_from = Library.ID,values_from = TPM) %>% 
    replace(is.na(.), 0)
}



validation_lines<-function(dat1,dat2,fill1,fill2,types){

  #dat1 = clown_orth_net_in
  #dat2 = full_thyroid_orths
  
  d.clown<-dat1$dge_keep$counts %>% 
    edgeR::cpm(log = T,prior.count = 1e-5,normalized.lib.sizes = T)  %>% 
    reshape2::melt() %>% 
    dplyr::rename(orths = Var1,stages = Var2, logCPM = value) %>% 
    filter(orths %in% dat2$Orthogroup) %>% group_by(orths) %>% 
    mutate(z.score = (logCPM - mean(logCPM))/sd(logCPM)) %>% 
    ungroup()
  
  d.clown<-d.clown %>% mutate(z.clipped = case_when(
    z.score > quantile(d.clown$z.score,0.95) ~ quantile(d.clown$z.score,0.95),
    z.score < -quantile(d.clown$z.score,0.95) ~ -quantile(d.clown$z.score,0.95),
    T ~ z.score)) %>% mutate(sp="clown")
  
  d.clown %>% mutate(orths=as.character(orths)) %>% 
    left_join(dat2, by=c("orths"="Orthogroup"), relationship = "many-to-many") %>% group_by(orths,stages,external_gene_name) %>% slice_max(logCPM)->d.clown

d.clown %>% 
  ggplot(aes(x=stages, y=z.score, group=external_gene_name))+
  facet_wrap(vars(external_gene_name), nrow =5)+
  geom_point(shape=21, size=2, color="grey10", fill=fill1)+
  geom_line(aes(x = stages, y = z.score, group = 1),linewidth = 1, color = fill2)+
  ylim(-2.2,2.2)+
  ylab("Z-score")+
  labs(subtitle =paste0(types))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1),
        strip.background = element_rect(fill = "grey95")) %>% return()
}

#Get the gene copy data for shapes species
gene_copy_comp_trait_comp<-function(col){
  #col="black"
orth_gene_copy %>% filter(Orthogroup %in% (cafe_genefamily_outs %>% filter(mod==col) %>% 
                                             pull(genes))) %>% 
  filter(sp %in% no_WGD_tree$tip.label) %>% 
  dplyr::select(-c(Family, Order, SpecCode)) %>% 
  pivot_wider(names_from = Orthogroup, values_from = no_of_gene_copies) %>% 
  replace(is.na(.), 0) %>% filter(sp %in% shapes_dat$sp)->sp_genes_wide
colnames(sp_genes_wide)<-c("sp", paste0("gene",seq(1:(ncol(sp_genes_wide)-1))))

expr_scaled<-sp_genes_wide %>% column_to_rownames("sp")->pca_dat_gene
scale(pca_dat_gene, center = FALSE,scale = TRUE)->pca_dat_gene
pca_gene <- phyl.pca(no_WGD_tree, pca_dat_gene, method = "lambda", mode = "corr")
pca_gene <- as.data.frame(pca_gene$S) %>% dplyr::rename_with(~paste0(.,".cpy")) %>% 
  rownames_to_column("sp")

shapes_dat %>% dplyr::select(sp, length, body_depth, fish_width) %>% column_to_rownames("sp")->pca_dat
scale(pca_dat,center = FALSE,scale = TRUE)->pca_dat
pca_res <- phyl.pca(no_WGD_tree, pca_dat, method = "lambda", mode = "corr")
pca_shapes <- as.data.frame(pca_res$S) %>% dplyr::rename_with(~paste0(.,".trait")) %>% 
  rownames_to_column("sp")

left_join(pca_gene,pca_shapes)->genes_and_morph

shapes_dat %>% dplyr::select(sp, lower_jaw_length, mounth_width) %>% column_to_rownames("sp")->pca_dat
scale(pca_dat,center = FALSE,scale = TRUE)->pca_dat
pca_res <- phyl.pca(no_WGD_tree, pca_dat, method = "lambda", mode = "corr")
pca_jaw <- as.data.frame(pca_res$S) %>% dplyr::rename_with(~paste0(.,".trait")) %>% 
  rownames_to_column("sp")

left_join(pca_gene,pca_jaw)->genes_and_jaw

shapes_dat %>% dplyr::select(sp, min_caudalpeduncle_depth, min_caudalpeduncle_width) %>% column_to_rownames("sp")->pca_dat
scale(pca_dat,center = FALSE,scale = TRUE)->pca_dat
pca_res <- phyl.pca(no_WGD_tree, pca_dat, method = "lambda", mode = "corr")
pca_pedancle <- as.data.frame(pca_res$S) %>% dplyr::rename_with(~paste0(.,".trait")) %>% 
  rownames_to_column("sp")

left_join(pca_gene,pca_pedancle)->genes_and_pedancle

library(patchwork)
library(ggpubr)
genes_and_morph %>% 
  ggplot(aes(x=PC1.trait, y=PC1.cpy))+
  geom_point(fill="#8BAE06",size=2, shape=21, alpha=0.7)+
  labs(
    x="trait phyPC1",
    y="copy number phyPC1")+
  geom_smooth(method = "lm", se = TRUE,colour = "#F0710F")+
  stat_cor(method = "spearman",cor.coef.name = "rho")+
  theme_linedraw()->point
genes_and_morph %>% 
  ggplot(aes(x=PC1.trait))+
  geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
  theme_void()->trait_den
genes_and_morph %>% 
  ggplot(aes(x=PC1.cpy))+
  geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
  coord_flip()+
  theme_void()->gene_den

p_blank <- ggplot() + theme_void()

(trait_den+p_blank+point+gene_den) +
  plot_layout(
    ncol   = 2,
    nrow   = 2,
    widths = c(4, 1),
    heights= c(1, 4))+
  plot_annotation(subtitle = "Morphology")->p1.1

genes_and_morph %>% 
  ggplot(aes(x=PC1.trait, y=PC2.cpy))+
  geom_point(fill="#8BAE06",size=2, shape=21, alpha=0.7)+
  labs(
    x="trait phyPC1",
    y="copy number phyPC2")+
  geom_smooth(method = "lm", se = TRUE,colour = "#F0710F")+
  stat_cor(method = "spearman",cor.coef.name = "rho")+
  theme_linedraw()->point2
genes_and_morph %>% 
  ggplot(aes(x=PC1.trait))+
  geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
  theme_void()->trait_den
genes_and_morph %>% 
  ggplot(aes(x=PC2.cpy))+
  geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
  coord_flip()+
  theme_void()->gene_den2

p_blank <- ggplot() + theme_void()

(trait_den+p_blank+point2+gene_den2) +
  plot_layout(
    ncol   = 2,
    nrow   = 2,
    widths = c(4, 1),
    heights= c(1, 4))+
  plot_annotation(subtitle = "Morphology")->p1.2
  
genes_and_jaw %>% 
    ggplot(aes(x=PC1.trait, y=PC1.cpy))+
  geom_point(fill="#8BAE06",size=2, shape=21, alpha=0.7)+
  labs(
    x="trait phyPC1",
    y="copy number phyPC1")+
  geom_smooth(method = "lm", se = TRUE,colour = "#F0710F")+
  stat_cor(method = "spearman",cor.coef.name = "rho")+
  theme_linedraw()->point
genes_and_jaw %>% 
    ggplot(aes(x=PC1.trait))+
    geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
    theme_void()->trait_den
genes_and_jaw %>% 
    ggplot(aes(x=PC1.cpy))+
    geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
    coord_flip()+
    theme_void()->gene_den
  
  p_blank <- ggplot() + theme_void()
  
  (trait_den+p_blank+point+gene_den) +
    plot_layout(
      ncol   = 2,
      nrow   = 2,
      widths = c(4, 1),
      heights= c(1, 4))+
    plot_annotation(subtitle = "Head&Jaw")->p2.1
  
  genes_and_jaw %>% 
    ggplot(aes(x=PC1.trait, y=PC2.cpy))+
    geom_point(fill="#8BAE06",size=2, shape=21, alpha=0.7)+
    labs(
      x="trait phyPC1",
      y="copy number phyPC2")+
    geom_smooth(method = "lm", se = TRUE,colour = "#F0710F")+
    stat_cor(method = "spearman",cor.coef.name = "rho")+
    theme_linedraw()->point2
  genes_and_jaw %>% 
    ggplot(aes(x=PC1.trait))+
    geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
    theme_void()->trait_den
  genes_and_jaw %>% 
    ggplot(aes(x=PC2.cpy))+
    geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
    coord_flip()+
    theme_void()->gene_den2
  
  (trait_den+p_blank+point2+gene_den2) +
    plot_layout(
      ncol   = 2,
      nrow   = 2,
      widths = c(4, 1),
      heights= c(1, 4))+
    plot_annotation(subtitle = "Head&Jaw")->p2.2

genes_and_pedancle %>% 
    ggplot(aes(x=PC1.trait, y=PC1.cpy))+
  geom_point(fill="#8BAE06",size=2, shape=21, alpha=0.7)+
  labs(
    x="trait phyPC1",
    y="copy number phyPC1")+
  geom_smooth(method = "lm", se = TRUE,colour = "#F0710F")+
  stat_cor(method = "spearman",cor.coef.name = "rho")+
  theme_linedraw()->point
genes_and_pedancle %>% 
    ggplot(aes(x=PC1.trait))+
    geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
    theme_void()->trait_den
genes_and_pedancle %>% 
    ggplot(aes(x=PC1.cpy))+
    geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
    coord_flip()+
    theme_void()->gene_den
  
  p_blank <- ggplot() + theme_void()
  
  (trait_den+p_blank+point+gene_den) +
    plot_layout(
      ncol   = 2,
      nrow   = 2,
      widths = c(4, 1),
      heights= c(1, 4))+
    plot_annotation(subtitle = "Caudal")->p3.1
  
  genes_and_pedancle %>% 
    ggplot(aes(x=PC1.trait, y=PC2.cpy))+
    geom_point(fill="#8BAE06",size=2, shape=21, alpha=0.7)+
    labs(
      x="trait phyPC1",
      y="copy number phyPC2")+
    geom_smooth(method = "lm", se = TRUE,colour = "#F0710F")+
    stat_cor(method = "spearman",cor.coef.name = "rho")+
    theme_linedraw()->point2
  genes_and_pedancle %>% 
    ggplot(aes(x=PC1.trait))+
    geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
    theme_void()->trait_den
  genes_and_pedancle %>% 
    ggplot(aes(x=PC2.cpy))+
    geom_density(fill="#8BAE06", alpha=0.4, linewidth=0.3)+
    coord_flip()+
    theme_void()->gene_den2
  
  p_blank <- ggplot() + theme_void()
  
  (trait_den+p_blank+point2+gene_den2) +
    plot_layout(
      ncol   = 2,
      nrow   = 2,
      widths = c(4, 1),
      heights= c(1, 4))+
    plot_annotation(subtitle = "Caudal")->p3.2
  

ggpubr::ggarrange(p1.1,p1.2,
                  p2.1,p2.2,
                  p3.1,p3.2,
                  nrow = 3,ncol = 2) %>% 
  annotate_figure(bottom = text_grob(paste0(col), face = "bold", size = 14)) %>% 
  return()
}


#CAFE processing functions
cafe_trait_tree<-function(tree, change_dat, model, prune=FALSE, with_tips=FALSE){
  library(ggtree)
  library(ggnewscale)
  
    concat_trees<-function(trees_in,name,tips,prune){
      tree<-trees_in
     
      
      #This is used to get the gene copy number at the root
      (str_split(tree$node.label[1],"_",simplify = TRUE))[,2] %>% 
        as.numeric()->root_node_copy
      
      #Prune the tree to only include species with trait data when option set to TRUE
      pruned_tree.model <- if (isTRUE(prune)) {
        drop.tip(
          tree,
          tree$tip.label[-match(
            tips,
            str_split(tree$tip.label, "</*", simplify = TRUE)[,1]
          )]
        )
      } else {
        tree
      }
      
      
      #intialise the tree
      ggtree::ggtree(pruned_tree.model,layout="rectangular")->g
      gene_copy_levels <- c("0","1","2","3","4","5",">5")#Copy number variable
      
      #Change the CAFE naming format
      g$data %>%
        tidyr::extract(
          label, into = c("species", "angle_tag", "suffix"),
          regex = "^([^<]*)(<\\d+>)(?:\\*)?(_\\d+)$", remove = FALSE) %>%
        mutate(species=na_if(species, ""),
               gene_copy=dplyr::case_when(
                 stringr::str_detect(label, "<135>") ~ root_node_copy,
                 TRUE ~ readr::parse_number(suffix)),
               #collapse to seven display levels
               gene_copy_disp = case_when(
                 is.na(gene_copy) ~ NA_character_,
                 gene_copy > 5    ~ ">5",
                 TRUE ~ as.character(gene_copy)),
               gene_copy_disp = factor(gene_copy_disp, 
                                       levels = gene_copy_levels),
               isSig = ifelse(stringr::str_detect(label, "\\*"), "yes", "no")) %>% 
        mutate(join_col=paste0(species,angle_tag) %>% str_remove_all("NA")) %>% 
        dplyr::select(-c(suffix, angle_tag)) %>% 
        mutate(genes=paste0(name))->g$data
      g$data %>% return()

} 

    model$data$sp->include.tips.model
  
purrr::imap(.x= tree,
            ~concat_trees(trees_in = .x, name = .y, 
                               tips = include.tips.model,
                          prune = prune)) %>% 
   list_rbind()-> tree_data_list
 

change_dat %>% pivot_longer(-FamilyID,names_to = "node", values_to = "change") %>% 
  mutate(genes = str_remove_all(FamilyID, "OG")) %>% 
  filter(genes %in% cafe_genefamily_outs$genes) %>% mutate(
    sign = case_when(
      change >  0 ~ "expansion",
      change < 0 ~ "contraction",
      TRUE          ~ "remains"
    )) %>%
  dplyr::count(node, sign, genes) %>%
  pivot_wider(
    names_from  = sign,
    values_from = n,
    values_fill = 0) %>% 
  mutate(type_of_change=case_when(
    expansion > 0 ~ "expansion", 
    contraction > 0 ~ "contraction", 
    .default = "remains"))->exp_cont_dat

ggtree::ggtree(tree_data_list,layout="rectangular")->g

  left_join(tree_data_list,exp_cont_dat %>%
            dplyr::select(node, genes, type_of_change) %>% 
              mutate(genes=paste0("OG",genes)) %>% 
              dplyr::rename("join_col"="node"),
          by=c("genes","join_col"), relationship = "many-to-many")->g$data


#define  palette 
fill_values_cp <- c(
  "0"  = "grey",
  "1"  = "#A4BFA8",
  "2"  = "#2A8C5E",
  "3"  = "#F2CB05",
  "4"  = "#735702",
  "5"  = "#F25252",   
  ">5" = "#F25252"
  )

  g <- if (isTRUE(with_tips)) {
    ggtree::ggtree(g$data, layout = "rectangular") + geom_tree() +
    geom_tiplab(aes(label=str_replace_all(species,pattern = "_",replacement = " ")),as_ylab = FALSE,offset = 5, size=3,geom = "text", fontface="italic")+ hexpand(.7)
    #g +geom_tiplab(size=3)
    } else {
    g
  }
  
  g+
    facet_wrap(~genes)+
    geom_tree(color = "grey80", show.legend = FALSE) +
    geom_tree(data = function(x) dplyr::filter(x, !isTip),
              aes(color = type_of_change),
              show.legend = TRUE)+
    scale_color_manual(values = c("#5430D9","#D90467","grey"))+
    ggnewscale::new_scale_color()+
    geom_point(data = function(x) dplyr::filter(x, !isTip),
               aes(x=x, y=y, 
                 fill=gene_copy_disp,
                 color=isSig,
                 stroke=ifelse(isSig == "yes", 0.5, 0),
                 size=ifelse(isSig == "yes", 2, 1.5)),
             shape=21)+
  scale_size_identity()+
  scale_fill_manual(values=fill_values_cp)+
  scale_color_manual(values = c("white","black")) %>% return()

}

cafe_get_age_of_change<-function(tree, change_dat, cafe_genefamily_res){
  library(ggtree)
  library(ggnewscale)
  
  trees <- if (inherits(tree, "phylo")) list(tree) else tree
  if (is.null(names(trees))) names(trees) <- paste0("OG", seq_along(trees))
  
  concat_trees<-function(trees_in,name){
    tr<-trees_in
    
    
    #This is used to get the gene copy number at the root
    (str_split(tr$node.label[1],"_",simplify = TRUE))[,2] %>% 
      as.numeric()->root_node_copy
    
    
    #intialise the tree
    ggtree::ggtree(tr,layout="rectangular")->g
    gene_copy_levels <- c("0","1","2","3","4","5",">5")#Copy number variable
    
    #Change the CAFE naming format
    g$data %>%
      tidyr::extract(
        label, into = c("species", "angle_tag", "suffix"),
        regex = "^([^<]*)(<\\d+>)(?:\\*)?(_\\d+)$", remove = FALSE) %>%
      mutate(species=na_if(species, ""),
             gene_copy=dplyr::case_when(
               stringr::str_detect(label, "<135>") ~ root_node_copy,
               TRUE ~ readr::parse_number(suffix)),
             #collapse to seven display levels
             gene_copy_disp = case_when(
               is.na(gene_copy) ~ NA_character_,
               gene_copy > 5    ~ ">5",
               TRUE ~ as.character(gene_copy)),
             gene_copy_disp = factor(gene_copy_disp, 
                                     levels = gene_copy_levels),
             isSig = ifelse(stringr::str_detect(label, "\\*"), "yes", "no")) %>% 
      mutate(join_col=paste0(species,angle_tag) %>% str_remove_all("NA")) %>% 
      dplyr::select(-c(suffix, angle_tag)) %>% 
      mutate(genes=paste0(name))->g$data
    g$data %>% return()
    
  } 
  
  purrr::imap(.x= trees,
              ~concat_trees(trees_in = .x, name = .y)) %>% 
    list_rbind()-> tree_data_list
  
  
  change_dat %>% pivot_longer(-FamilyID,names_to = "node", values_to = "change") %>% 
    mutate(genes = str_remove_all(FamilyID, "OG")) %>% 
    filter(genes %in% cafe_genefamily_res$genes) %>% mutate(
      sign = case_when(
        change >  0 ~ "expansion",
        change < 0 ~ "contraction",
        TRUE          ~ "remains"
      )) %>%
    dplyr::count(node, sign, genes) %>%
    pivot_wider(
      names_from  = sign,
      values_from = n,
      values_fill = 0) %>% 
    mutate(type_of_change=case_when(
      expansion > 0 ~ "expansion", 
      contraction > 0 ~ "contraction", 
      .default = "remains"))->exp_cont_dat
  
  ggtree::ggtree(tree_data_list,layout="rectangular")->g
  
  left_join(tree_data_list,exp_cont_dat %>%
              dplyr::select(node, genes, type_of_change) %>% 
              mutate(genes=paste0("OG",exp_cont_dat$genes)) %>% 
              dplyr::rename("join_col"="node"),
            by=c("genes","join_col"), relationship = "many-to-many")->g$data
  
  g$data %>% return()
  
  
}


phylo_pca_plot<-function(pp, traitlabel){
  
  #pp<-model_brain$pink
  eig  <- pp$pca_gene_pca$Eval
  load <- pp$pca_gene_pca$Evec
  scores <- pp$pca_gene_pca$S
  trait_dat <- pp$model_brms_brain$data
  
  # ensure names for PCs
  pc_names <- paste0("PC", seq_along(eig))
  colnames(load)   <- colnames(load)   %||% pc_names
  colnames(scores) <- colnames(scores) %||% pc_names
  
  var_exp<-eig/sum(eig) * 100  # % variance explained
  
  scores_df <- as.data.frame(scores) %>% 
    tibble::rownames_to_column("species") %>% 
    dplyr::select(species, PC1, PC2)
  
  load_df <- as.data.frame(load) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::select(gene, PC1, PC2)
  
  # Scale eigenvector arrows to fit on the same axes as scores (simple, effective)
  s_range <- max(diff(range(scores_df$PC1)), diff(range(scores_df$PC2)))
  l_range <- max(diff(range(load_df$PC1)),   diff(range(load_df$PC2)))
  arrow_scale <- 0.7 * s_range / l_range
  arrow_df <- load_df %>% 
    mutate(PC1 = PC1 * arrow_scale,
           PC2 = PC2 * arrow_scale)
  
 scores_df<-left_join(trait_dat %>% dplyr::select(2,3), scores_df,by = c("sp"="species"))
  ## 1) Species scores (PC1 vs PC2) -------------------------------------------
 scores_df %>% mutate(sp=str_replace_all(sp,pattern = "_",replacement = " ")) %>% 
   ggplot(aes(PC1, PC2, fill=BrainWeight)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(size = 4, alpha = 0.8, shape=21,
               position = position_jitter(width = 0.8, height = 0.8, seed = 1)) +
    geom_text_repel(aes(label = sp), size = 3, max.overlaps = Inf,fontface = "italic") +
   scale_fill_viridis_c(option = "G",direction = -1, name = paste0(traitlabel))+
   geom_segment(inherit.aes = FALSE,
     data = arrow_df,
     aes(x = 0, y = 0, xend = PC1, yend = PC2),
     arrow = arrow(length = unit(0.02, "npc")),
     linewidth = 0.4,alpha=0.5)+
   geom_text_repel(inherit.aes = FALSE,data = arrow_df,
                   aes(PC1, PC2, label = gene),
                   size = 3, max.overlaps = Inf, 
                   color="red",alpha=0.8)+
   labs(x=paste0("phyPC1 (",var_exp[1,1] %>% round(1),"%)"),
        y=paste0("phyPC2 (",var_exp[2,2] %>% round(1),"%)"))+
    theme_linedraw()+theme(legend.position = "top") %>% return()
}
#Functions to make plots of brms model fit
brms_significance_plot<-function(model_list, submodel, trait, var_type) {
    
    
    ##Use to test:
    #model_list<- model_aspect_ratio
    #submodel<-"model_brms_aspect_ratio"
    #trait<-"aspect_ratio"
    #var_type<-"continuous"
    
   model_plot_dat <- if (var_type=="categorical") {
     map_dfr(
       model_list,
       ~{fit <- purrr::pluck(.x, submodel, "fit", .default = NULL)
       summary(fit, probs = c(0.025, 0.975))$summary %>%
         as.data.frame() %>%
         rownames_to_column("model_vars")},
       .id = "model_id") %>%
       filter(str_detect(model_vars, "\\.cpy")) %>% 
       mutate(feature=str_split(model_vars, pattern = "_",simplify = T)[,2],
              eigengene=str_split(model_vars, pattern = "_",simplify = T)[,3]) %>%
       mutate(feature=reorder(feature,`97.5%`)) 
     } else if (var_type=="continuous") {
       map_dfr(
         model_list,
         ~{fit <- purrr::pluck(.x, submodel, "fit", .default = NULL)
         summary(fit, probs = c(0.025, 0.975))$summary %>%as.data.frame() %>%
           rownames_to_column("model_vars")},
         .id = "model_id") %>%
        filter(str_detect(model_vars, "\\.cpy")) %>% 
        mutate(feature=trait,
               eigengene=str_split(model_vars, pattern = "_",simplify = T)[,2]) %>%
        mutate(feature=reorder(feature,`97.5%`))
     } else {
       stop("var_type must be 'categorical' or 'continuous'")
       }
    
   model_plot_dat %>%
      mutate(eigengene  = factor(eigengene),
             model_id   = factor(model_id),
             sig95      = (`2.5%` >= 0 & `97.5%` >= 0 | `2.5%` <= 0 & `97.5%` <= 0)# CI doe not overlap 0
             ) %>% mutate(feature=reorder(feature, sig95)) %>% 
     filter(feature != "Intercept")->plot_dat
    
    left_join(plot_dat,module_annot,by = c("model_id"="mod"))->plot_dat 
    plot_dat %>% dplyr::select(annot,sig95,eigengene,`2.5%`,`97.5%`,n_eff,Rhat) %>% arrange(desc(sig95)) %>% return()
    
    #Can use the code below for a plot
    #ggplot(plot_dat, aes(y = mean, x = feature, group = eigengene, fill = sig95)) +
      #facet_grid(rows = vars(annot), cols = vars(eigengene), scales = "free", space = "free_x") +
      #geom_hline(yintercept = 0, linetype = 2, color="red") +
      #geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
       #             width = 0.15, 
        #            alpha = 0.9) +
      #geom_point(aes(shape = sig95, color=sig95), size = 2) +
      #scale_color_manual(values = c("#6B9650","#E674EA"))+
      #labs(x = "Model estimate",
       #    y = NULL) +
      #theme_classic() +
      #theme(panel.spacing = unit(5, "pt"),
       #     strip.background = element_rect(fill = "grey95"),
        #    strip.text.x.top =  element_text(angle = 0, size = 8),
         #   strip.text.y.right =  element_text(angle = 360, size = 6,hjust = -0.01),
          #  legend.position = "bottom",
           # axis.text.x = element_text(angle = 90,
            #                           hjust = 1, 
             #                          vjust = 0.3,
              #                         size = 8))->p
      
  }


make_pca_trait_model_plot<-function(model,pca_dat,trait,arrow_scale,fill_label){
    
    #Use for testing
    #model = model_aspect_ratio$salmon$model_brms_aspect_ratio
    #pca_dat = model_aspect_ratio$salmon$pca_gene_pca
    #trait = "aspect_ratio"
    #arrow_scale = 30
    #fill_label = "Scaled aspect ratio"
    
    as_tibble(diag(pca_dat$Eval)/sum(pca_dat$Eval)*100) %>% 
      mutate(PCs=rownames(pca_dat$Eval)) %>% 
      pivot_wider(names_from = PCs, values_from = value)->var_exp
    
    #pca_dat$Evec %>% as_tibble() %>% 
    # mutate(genes=rownames(pca_dat$Evec))->vec_dat
    
    pca_dat$Evec %>% as.tibble() %>% 
      mutate(genes=rownames(pca_dat$Evec)) %>% 
      left_join(danio_orth_mods,by=c("genes"="gene")) %>% 
      distinct(genes, .keep_all = TRUE) %>% 
      dplyr::rename("gene.name"=`Gene name`)->vec_dat
    
    #this is the data for the geom segment
      arrow_tip_dat <- vec_dat %>%
        dplyr::mutate(
          xend = PC1 * arrow_scale,
          yend = PC2 * arrow_scale,
        #offset for label
        r = sqrt(xend^2 + yend^2),
        x_lab = xend + (xend / r) * 0.6,   # 0.6 controls how far from tip
        y_lab = yend + (yend / r) * 0.6)
    
    #remove '_'
    trait_label<-str_replace_all(trait,pattern = "_", replacement = " ")
    
    if (is.numeric(model$data[[trait]])) {
      
      model$data %>% 
        mutate(z=(model$data[[trait]] - mean(model$data[[trait]]))/sd(model$data[[trait]])) %>%
        mutate(z.clipped=case_when(
          z > quantile(z, 0.95) ~ quantile(z, 0.95), 
          z < -quantile(z, 0.95) ~ -quantile(z, 0.95),
          .default = z)) %>% 
        mutate(sp_shrt=str_match(sp, "^([A-Za-z]{2})[A-Za-z]*_([A-Za-z]{3})[A-Za-z]*$") %>%
                 { paste0(.[,2], .[,3]) })->scaled_dat
      
      qnt<-quantile(scaled_dat$z,0.95) #used for plotting
      
      #Plot properties for labels
      pal_fun <- scales::col_numeric(palette = c("#F5D33D", "#e4e4e4", "#346168"),
                                     domain  = c(-qnt, qnt),na.color = "grey80")
      
      #choose text color based on perceived luminance
      pick_text_contrast <- function(hex, threshold = 0.6) {
        rgb <- grDevices::col2rgb(hex) / 255
        rgb_lin <- ifelse(
          rgb <= 0.03928,
          rgb / 12.92,
          ((rgb + 0.055) / 1.055)^2.4)
        
        #relative luminance
        lum <- 0.2126 * rgb_lin[1, ] + 0.7152 * rgb_lin[2, ] + 0.0722 * rgb_lin[3, ]
        ifelse(lum < threshold, "white", "black")
      }
      
      scaled_dat %>%
        mutate(
          z_squished   = pmax(pmin(z.clipped, qnt), -qnt),   # match oob = squish
          label_fill   = pal_fun(z_squished),
          label_text   = pick_text_contrast(label_fill)) %>% 
        
        ggplot(aes(x=PC1.cpy, y= PC2.cpy))+
        geom_segment(aes(0,0, xend=PC1*arrow_scale, yend=PC2*arrow_scale), 
                     data = vec_dat,
                     arrow = arrow(type = "closed", length = unit(0.2, "cm")), 
                     color="#D3004F", 
                     size=0.5, 
                     show.legend = FALSE, 
                     inherit.aes = FALSE)+
        geom_label(aes(label = sp_shrt, fill  = z.clipped, color = label_text),
                   position = position_jitter(width = 3, height = 3, seed = 42),
                   alpha=0.9,
                   size = 2.5,                      # small text
                   linewidth = 0.12,               # thin border
                   label.padding = unit(0.06, "lines")) +
        scale_fill_gradientn(colours = c("#F5D33D","#e4e4e4", "#346168"),
                             limits = c(-qnt,qnt),
                             breaks = c(-qnt, 0, qnt), 
                             labels  = c(paste0("<", round(-qnt)), "0", paste0(">", round(qnt))),
                             oob = scales::squish)+
        scale_color_identity()+   #precomputed color values
        geom_text(aes(x = x_lab, y = y_lab, label = gene.name),
                  data = arrow_tip_dat,
                  color="#D3004F",
                  size = 3,
                  inherit.aes = FALSE)+
        labs(x=paste0("PC1 gene copy","(",round(var_exp$PC1),"%)"), 
             y=paste0("PC2 gene copy","(",round(var_exp$PC2),"%)"),
             fill=paste0(fill_label))+
        theme_bw()->pc_plot
      
      
      
      scaled_dat  %>% 
        ggplot(aes(y=PC1.cpy,x=scaled_dat[[trait]]))+
        geom_point(aes(fill=z.clipped),shape=21,size=3,stroke=0.3)+
        geom_smooth(method = "glm", color="#D3004F",alpha = 0.3,fill = "grey70")+
        scale_fill_gradientn(colours = c("#F5D33D","#e4e4e4", "#346168"),
                             limits = c(-qnt,qnt),
                             breaks = c(-qnt, 0, qnt), 
                             labels  = c(paste0("<", round(-qnt)), "0", paste0(">", round(qnt))),
                             oob = scales::squish)+
        labs(x=paste0(trait_label))+
        theme_bw()+
        theme(legend.position = "none")->pc1_genes_plot
      
      
      scaled_dat %>% 
        ggplot(aes(y=PC2.cpy,x=scaled_dat[[trait]]))+
        geom_point(aes(fill=z.clipped),shape=21,size=3,stroke=0.3)+
        geom_smooth(method = "glm", color="#D3004F",alpha = 0.3,fill = "grey70")+
        scale_fill_gradientn(colours = c("#F5D33D","#e4e4e4", "#346168"),
                             limits = c(-qnt,qnt),
                             breaks = c(-qnt, 0, qnt), 
                             labels  = c(paste0("<", round(-qnt)), "0", paste0(">", round(qnt))),
                             oob = scales::squish)+
        labs(x=paste0(trait_label))+
        theme_bw()+
        theme(legend.position = "none")->pc2_genes_plot
      #get the legend
      leg_grob <- leg_grob <- get_legend(
        pc_plot +
          theme(
            legend.title = element_text(size=9.5),
            legend.position   = "bottom",
            legend.background = element_blank(),
            legend.key.width  = unit(0.8, "cm"),
            legend.key.height = unit(0.25, "cm"),
            legend.text       = element_text(size = 9)
          ) +
          guides(fill = guide_colorbar(
            barwidth  = unit(4, "cm"),
            barheight = unit(0.35, "cm"),
            title.position = "right",
            label.position = "bottom",
            label.theme = element_text(size = 8, margin = margin(t = -8)) # move labels into bar
          )))
      
      #remove legend from PCA plot
      pc_plot_noleg  <- pc_plot  + theme(legend.position = "none")
      
      #get top row
      top_row <- ggarrange(pc1_genes_plot, pc2_genes_plot, ncol = 2)
      
      #combine
      ggarrange(
        top_row,
        leg_grob,
        pc_plot_noleg,
        ncol = 1,
        heights = c(1, 0.1, 2)   # adjust middle number to change legend row height
      ) 
    } else {
      
      model$data %>%  
        mutate(sp_shrt=str_match(sp, "^([A-Za-z]{2})[A-Za-z]*_([A-Za-z]{3})[A-Za-z]*$") %>%
                 { paste0(.[,2], .[,3]) }) %>% 
        ggplot(aes(x=PC1.cpy, y= PC2.cpy))+
        geom_segment(aes(0,0, xend=PC1*arrow_scale, yend=PC2*arrow_scale), data = vec_dat,
                     arrow = arrow(type = "closed", 
                                   length = unit(0.2, "cm")), 
                     color="#D3004F", size=0.5, show.legend = FALSE, inherit.aes = FALSE)+
        geom_label(aes(label = sp_shrt, fill = model$data[[trait]]),
                   position = position_jitter(width = 3, height = 3, seed = 42),
                   color = "white",
                   alpha=0.9,
                   size = 3,                     
                   linewidth = 0.10,               
                   label.padding = unit(0.06, "lines")) +
        ggsci::scale_fill_jama()+
        geom_text(aes(x = x_lab, y = y_lab, label = gene.name),
                  data = arrow_tip_dat,
                  color="#D3004F",
                  size = 3,
                  inherit.aes = FALSE)+
        labs(x=paste0("PC1 gene copy","(",round(var_exp$PC1),"%)"), 
             y=paste0("PC2 gene copy","(",round(var_exp$PC2),"%)"), 
             fill=paste0(fill_label))+
        theme_bw()
      
    }
  }


#Processing Bgee data
Bgee_dge<-function(.data){
  library(edgeR)
#.data=bgee_dat$Danio_rerio
  
design_matrix<-.data  %>% 
                 dplyr::select(Experiment.ID,Library.ID,Anatomical.entity.name) %>% distinct()
  
  dge<-edgeR::DGEList(.data %>% make_exp_data_frames() %>% column_to_rownames("Gene.ID"), group = factor(design_matrix$Anatomical.entity.name))
    return(dge$counts %>% as_tibble() %>% mutate(Gene.ID=rownames(dge$counts)))
    
}
#Building the tau data frame 
build_exp_data_frame<-function(.data, Species) {
  left_join(.data %>% pivot_longer(!Gene.ID,names_to = "Library.ID",values_to = "cpm"),
            bgee_dat$Danio_rerio %>% 
              dplyr::select(Library.ID,Anatomical.entity.name),by = "Library.ID") %>% 
    group_by(Gene.ID,Anatomical.entity.name) %>% 
    summarise(cpm = mean(cpm)) %>% 
    spread(Anatomical.entity.name,cpm) %>% ungroup() %>% return()
}

bgee_orth_json<-function(orth){
  library(jsonlite)
  f<-fromJSON(paste0("./11.Bgee_comparative_transcriptomics/",orth,".json "),flatten = TRUE)
  f$data$orthologsByTaxon$genes %>% list_rbind() %>% 
    mutate(danio_gene=orth, sp=paste0(species.genus,"_",species.speciesName)) %>% return()
}

build_data_frame_from_bgee<-function(dat, species){
  library(edgeR)
  library(DESeq2)
  make_exp_data_frames_Read.count<-function(.data, ...){
    .data  %>% 
      dplyr::select(Gene.ID, Read.count, Library.ID) %>% 
      pivot_wider(names_from = Library.ID,values_from = Read.count) %>%
      mutate(across(everything(), ~replace_na(.x, 0)))
  }
  
  design_matrix<-dat %>% 
    dplyr::select(Experiment.ID,Library.ID,Anatomical.entity.name) %>% distinct()
  
  dge<-edgeR::DGEList(dat %>% make_exp_data_frames_Read.count() %>% 
                        column_to_rownames("Gene.ID"), 
                      group = factor(design_matrix$Anatomical.entity.name))
  keep<-filterByExpr(cpm(dge), min.count = 1, group = factor(design_matrix$Experiment.ID))
  dge<-dge[keep, ,keep.lib.sizes =T]
  dge<-calcNormFactors(dge,method = "TMM")
  dat_genes<-rownames(dge$counts)
  dge<-edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  dge %>% as_tibble() %>% 
    mutate(genes=dat_genes) %>% 
    pivot_longer(cols = c(-genes), names_to = "Library.ID", values_to = "cpm") %>%
    left_join(design_matrix, by = "Library.ID") %>% 
    filter(Anatomical.entity.name !="multicellular organism") %>% 
    mutate(sp=species) %>% return()
  }

 phylo_bgee_plot<-function(col){
 bgee_combined_long  %>% 
     distinct(sp) %>% pull(sp)->sp.to.get
 tibble(tips=tree$tip.label) %>% 
   filter(tips %in% sp.to.get) ->include.tips.bgee
 pruned_tree_bgee<-drop.tip(tree,tree$tip.label[-match(include.tips.bgee$tips, tree$tip.label)])
 #plot(pruned_tree_bgee)

 orths_in<-cafe_sig %>% filter(mod==col)
 
 genes_to_use <- bgee_orth_res_to_use %>%
   filter(Orthogroup %in% orths_in$genes) %>% pull(geneId)
 
  bgee_combined_long  %>% distinct(Anatomical.entity.name,sp) %>% 
   dplyr::count(Anatomical.entity.name) %>% 
   filter(n >=5)->tissues.to.keep

 
 # Make plot data
 plot_dat <- bgee_combined_long %>%
   filter(Gene.ID %in% genes_to_use)  %>% 
   group_by(sp, Anatomical.entity.name) %>%
   summarise(m_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
   mutate(log10_m_TPM = log10(m_TPM)) %>%
   group_by(sp) %>%
   ungroup()
 

 # Species order from tree
 #sp_order <- pruned_tree_bgee$tip.label %>%
   #str_replace_all("_", " ")
 
 plot_dat2 <- plot_dat %>%
   filter(Anatomical.entity.name %in% tissues.to.keep$Anatomical.entity.name) %>%
   dplyr::select(sp, Anatomical.entity.name, log10_m_TPM) %>%
   complete(
     sp = sp,
     Anatomical.entity.name = tissues.to.keep$Anatomical.entity.name
   ) %>%
   group_by(Anatomical.entity.name) %>%
   mutate(tissue_mean_expr = mean(log10_m_TPM, na.rm = TRUE)) %>%
   ungroup() %>%
   mutate(
     missing_cell = is.na(log10_m_TPM),
     Anatomical.entity.name =
       fct_reorder(Anatomical.entity.name, tissue_mean_expr, .desc = TRUE))
 
 
 lims <- quantile(plot_dat2$log10_m_TPM, c(0.05, 0.95), na.rm = TRUE)
 
 
 species_order <- c(
   "Nothobranchius_furzeri","Oryzias_latipes","Astatotilapia_calliptera",
   "Neolamprologus_brichardi","Scophthalmus_maximus","Gasterosteus_aculeatus",
   "Gadus_morhua","Salmo_salar","Esox_lucius","Danio_rerio",
   "Astyanax_mexicanus","Anguilla_anguilla","Lepisosteus_oculatus"
 )
 
 plot_dat2 <- plot_dat2 %>% left_join(
   bgee_combined_long %>%
     filter(Gene.ID %in% genes_to_use) %>% 
     group_by(sp) %>% summarise(genes_per_sp=n_distinct(Gene.ID)) %>% 
     ungroup(),
     by = "sp") %>%
   mutate(
     sp = factor(sp, levels = rev(species_order))
   )
 
 ggplot(plot_dat2, aes(x = Anatomical.entity.name, y = sp)) +
   # black circles for missing combinations
   geom_point(
     data = ~ dplyr::filter(.x, missing_cell),
     shape = 21,
     size = 2,
     colour = "grey20",
     stroke = 0.2
   ) +
   # filled circles for observed combinations
   geom_point(
     data = ~ dplyr::filter(.x, !missing_cell),
     aes(fill = log10_m_TPM, size = genes_per_sp),
     shape = 21,
     colour = "black",
     stroke = 0.1
   ) +
   scale_size_binned_area(max_size = 12)+
   scale_fill_gradientn(colours = c("#EEE4E4", "#7A5F80", "#8E3B73", "#D96F5E", "#F4C95D"),
     name = expression(log[10] * " mean TPM"),
     limits = lims,
     oob = scales::squish,
     breaks = pretty(lims, n = 5)
   ) + theme_minimal() +
   theme(
     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
     axis.title = element_blank(),
     panel.grid = element_blank())
 }
 
comp_dev_danio_heatmap<-function(system_proc){
  #read_csv("./08.ZFIN/Danio_biomart.csv")->Danio_biomart #Obtained from bioMart package
  #Add pseudocount and log
  ((danio_orth_exp_mat[,3:22])+1e-5) %>% log10->danio_orth_exp_mat_log
  danio_orth_exp_mat_log$orths=danio_orth_exp_mat$orths.danio
  #Subset to match data from the input for the zebrafih expression matrix
  
  danio_genes_modules_annot %>% 
    left_join(Danio_biomart %>% dplyr::select(`Gene stable ID`, `Gene name`),
              by = c(orths="Gene stable ID")) %>% 
    distinct(orths, .keep_all = TRUE) %>% 
    left_join(danio_orth_exp_mat_log,by="orths") %>% drop_na()->joined_danio
  
  #then only pick those genes in the above filtered data to visualise;
  #So we can compare expression at embryogenesis and post embryonic development
  
  danio_genes_modules_annot<-danio_genes_modules_annot %>% 
    filter(orths %in% joined_danio$orths)
  
  danio_muscle<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Muscle system") %>% pull(orths)))
  danio_matrix<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Extracellular matrix and cartilage development 2") %>% pull(orths)))
  danio_rna1<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="RNA processing 1") %>% pull(orths)))
  danio_rna2<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="RNA processing 2") %>% pull(orths)))
  danio_dna<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="DNA replication") %>% pull(orths)))
  danio_immune<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Immune system") %>% pull(orths)))
  danio_nutrient<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Nutrient absorption") %>% pull(orths)))
  danio_brush<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Brush border epithelium development") %>% pull(orths)))
  danio_synapse1<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Synaptic signalling 1") %>% pull(orths)))
  danio_synapse2<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Synaptic signalling 2") %>% pull(orths)))
  danio_peroxi<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Peroxisome process") %>% pull(orths)))
  danio_cellcylce<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Cell cycle process") %>% pull(orths)))
  danio_resp<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Aerobic respiration") %>% pull(orths)))
  danio_neuron<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Neuron migration and placode development") %>% pull(orths)))
  danio_tranl1<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Translation") %>% pull(orths)))
  danio_tranl2<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Translation 2") %>% pull(orths)))
  danio_ribo<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Ribosome related") %>% pull(orths)))
  danio_protcat<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Protein catabolism") %>% pull(orths)))
  danio_hormone<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Hormone response") %>% pull(orths)))
  danio_blood<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Blood related") %>% pull(orths)))
  danio_protfold<-Danio_biomart %>% filter(`Gene stable ID` %in% (danio_genes_modules_annot %>% filter(annot=="Protein folding") %>% pull(orths)))
  
  rowData(zhub_dat)$marker_class <- ifelse(rownames(zhub_dat) %in% danio_muscle$`Gene name`, "Muscle system",
                                         ifelse(rownames(zhub_dat) %in% danio_matrix$`Gene name`, "Extracellular matrix and cartilage development 2",
                                                ifelse(rownames(zhub_dat) %in% danio_rna1$`Gene name`, "RNA processing 1",
                                                       ifelse(rownames(zhub_dat) %in% danio_rna2$`Gene name`, "RNA processing 2",
                                                              ifelse(rownames(zhub_dat) %in% danio_dna$`Gene name`, "DNA replication",
                                                                     ifelse(rownames(zhub_dat) %in% danio_immune$`Gene name`, "Immune system",
                                                                            ifelse(rownames(zhub_dat) %in% danio_nutrient$`Gene name`, "Nutrient absorption",
                                                                                   ifelse(rownames(zhub_dat) %in% danio_brush$`Gene name`, "Brush border epithelium development",
                                                                                          ifelse(rownames(zhub_dat) %in% danio_synapse1$`Gene name`, "Synaptic signalling 1",
                                                                                                 ifelse(rownames(zhub_dat) %in% danio_synapse2$`Gene name`, "Synaptic signalling 2",
                                                                                                        ifelse(rownames(zhub_dat) %in% danio_peroxi$`Gene name`, "Peroxisome process",
                                                                                                               ifelse(rownames(zhub_dat) %in% danio_cellcylce$`Gene name`, "Cell cycle process",
                                                                                                                      ifelse(rownames(zhub_dat) %in% danio_resp$`Gene name`, "Aerobic respiration",
                                                                                                                             ifelse(rownames(zhub_dat) %in% danio_neuron$`Gene name`, "Neuron migration and placode development",
                                                                                                                                    ifelse(rownames(zhub_dat) %in% danio_tranl1$`Gene name`, "Translation",
                                                                                                                                           ifelse(rownames(zhub_dat) %in% danio_tranl2$`Gene name`, "Translation 2",
                                                                                                                                                  ifelse(rownames(zhub_dat) %in% danio_ribo$`Gene name`, "Ribosome related",
                                                                                                                                                         ifelse(rownames(zhub_dat) %in% danio_protcat$`Gene name`, "Protein catabolism",
                                                                                                                                                                ifelse(rownames(zhub_dat) %in% danio_hormone$`Gene name`, "Hormone response",
                                                                                                                                                                       ifelse(rownames(zhub_dat) %in% danio_blood$`Gene name`, "Blood related",
                                                                                                                                                                              ifelse(rownames(zhub_dat) %in% danio_protfold$`Gene name`, "Protein folding",
                                                                                                                                                                                     "other")))))))))))))))))))))
#Get pseudobulk
timepoint_mean <- try(aggregateAcrossCells(as(zhub_dat, "SingleCellExperiment"),  
                                           ids = zhub_dat$timepoint, 
                                           statistics = "sum",
                                           use.assay.type = "counts", 
                                           subset.row = rownames(zhub_dat)[rowData(zhub_dat)$marker_class == system_proc]))
#Set size factors
sizeFactors(timepoint_mean)<-colSums(assay(timepoint_mean, "counts"))/mean(colSums(assay(timepoint_mean, "counts")))
timepoint_mean<-logNormCounts(timepoint_mean)
rbind(
  joined_danio %>% filter(annot==paste0(system_proc)) %>% dplyr::select(`Gene name`, contains("danio")) %>% 
    dplyr::rename("genes"="Gene name") %>% mutate(dat_type="bulk") %>% 
    pivot_longer(cols = contains("danio"),names_to = "stages", values_to = "logcounts") %>% 
    mutate(stages=str_remove(stages,"\\.[0-9]")) %>% mutate(stages=str_remove(stages,".danio")) %>% 
    mutate(dat_type="post-embryonic") %>% group_by(genes) %>% 
    mutate(z.score=scale(logcounts)) %>% ungroup(),
  
  assay(timepoint_mean) %>%
    as.matrix() %>%as.data.frame() %>%
    rownames_to_column("genes") %>% dplyr::select(-`10dpf`) %>% 
    pivot_longer(cols = contains("pf"),
                 names_to = "stages", values_to = "logcounts") %>% 
    mutate(dat_type="embryonic") %>% group_by(genes) %>%
    mutate(z.score=scale(logcounts)) %>% ungroup()
  )->d

ggarrange(
  d %>% filter(dat_type=="embryonic") %>% group_by(genes) %>% 
    mutate(z.score=scale(logcounts)) %>% ungroup() %>% 
    ggplot(aes(x=stages, y=genes, fill=z.score))+
    facet_wrap(~dat_type, scales="free")+
    geom_tile()+
    scale_fill_gradient2(
      low = "#0072B2",
      mid = "#F7F7F7",
      high = "#D55E00",
      midpoint = 0,
      limits = c(-2, 2),
      oob = scales::squish,
      name = "Z-score")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)),
  
  d %>% filter(dat_type=="post-embryonic") %>% group_by(genes) %>%
    mutate(z.score=scale(logcounts)) %>% ungroup() %>%
    ggplot(aes(x=stages, y=genes, fill=z.score))+
    facet_wrap(~dat_type, scales="free")+
    geom_tile()+
    #scale_fill_continuous(palette = viridis::viridis(100,direction = -1, option = "B"))+
    scale_fill_gradient2(
      low = "#0072B2",
      mid = "#F7F7F7",
      high = "#D55E00",
      midpoint = 0,
      limits = c(-2, 2),
      oob = scales::squish,
      name = "Z-score")+
    labs(y=NULL)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)),
  ncol = 2,
  common.legend = TRUE, legend = "right")->plots
  
annotate_figure(plots, top = text_grob(paste0(system_proc),face = "bold", size = 12)) %>% return()
}

meta_exp_cafe_genes <- function(species, orthogroup){
  #dge_keep = grouper$dge_keep
  
  make_plot <- function(dge_keep, orthogroup, species, rects = list()) {
    danio_genes_modules_annot %>% left_join(Danio_biomart %>% dplyr::select(`Gene stable ID`, `Gene name`), by=c("orths"="Gene stable ID")) %>% distinct(orths, .keep_all = TRUE)->gene_annot
    
    samp <- dge_keep$samples %>%
      tibble::rownames_to_column("samples")
    
    df <- dge_keep$counts %>%
      edgeR::cpm(log = TRUE, prior.count = 1e-5, normalized.lib.sizes = TRUE) %>%
      reshape2::melt() %>%
      dplyr::rename(orths = Var1, samples = Var2, logCPM = value) %>%
      left_join(samp, by = "samples") %>%
      dplyr::rename(stages = group) %>%
      summarise(mean_logCPM = mean(logCPM, na.rm = TRUE),
                       .by = c(orths, stages)) %>% group_by(orths) %>%
      mutate(z.score = (mean_logCPM - mean(mean_logCPM))/sd(mean_logCPM)) %>%
      ungroup() %>% left_join(orth %>% filter(sp == species) %>%
          dplyr::select(Orthogroup, orths) %>% mutate(orths = stringr::str_remove(orths, "\\.[0-9]*")),
        by = "orths")
  
    p <- df %>% left_join(gene_annot, by = c("Orthogroup")) %>%
      dplyr::filter(Orthogroup == orthogroup) %>%
      ggplot(ggplot2::aes(x = stages, y = z.score)) +
      facet_wrap(~`Gene name`) +
      geom_line(aes(group = orths.x, color = orths.x), linewidth = 0.6) +
      ggsci::scale_color_nejm() +
      theme_bw() +
      theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            strip.text.x.top = element_text(size=7), legend.position = "none")
    
    # add the rectangles for the metamorphosis transitions
    for (r in rects) {
      p <- p + ggplot2::annotate(
        geom = "rect",
        xmin = r$xmin, xmax = r$xmax, ymin = -Inf, ymax = Inf,
        fill = r$fill %||% "#594BE4",
        alpha = r$alpha %||% 0.2
      )
    }
    
    p
  }
  

  if (species == "Amphiprion_ocellaris") {
    make_plot(
      dge_keep = clown$dge_keep,
      orth = orthogroup,
      species = species,
      rects = list(list(xmin = 4.5, xmax = 5.5, alpha = 0.2, fill = "#594BE4"))
    )
    
  } else if (species == "Epinephelus_malabaricus") {
    make_plot(
      dge_keep = grouper$dge_keep,
      orth = orthogroup,
      species = species,
      rects = list(
        list(xmin = 1.5, xmax = 2.5, alpha = 0.15, fill = "#594BE4"),
        list(xmin = 5.5, xmax = 6.5, alpha = 0.15, fill = "#594BE4")
      )
    )
    
  } 
}

 pheno_terms_plot<-function(col){
   
#col="salmon"
   
   pheno_orth_mods<-pheno_orth_mods %>% 
     distinct(`Gene Symbol`,`Phenotype Keyword ID`, .keep_all = TRUE)
   
   gs_dat <- pheno_orth_mods %>%
     filter(mod == col) %>%
     dplyr::count(`Gene Symbol`, 
                  `Affected Structure or Process 1 superterm Name`, 
                  name = "n")
   
   top_genes <- gs_dat %>%
     dplyr::count(`Gene Symbol`, wt = n, name = "total") %>%
     arrange(desc(total)) %>%
     slice_head(n = 100) %>%
     pull(`Gene Symbol`)
   
   top_structures <- gs_dat %>%
     dplyr::count(`Affected Structure or Process 1 superterm Name`, wt = n, name = "total") %>%
     arrange(desc(total)) %>%
     slice_head(n = 20) %>%
     pull(`Affected Structure or Process 1 superterm Name`)
   
   gs_plot <- gs_dat %>%
     filter(
       `Gene Symbol` %in% top_genes,
       `Affected Structure or Process 1 superterm Name` %in% top_structures
     ) %>%
     complete(
       `Gene Symbol` = top_genes,
       `Affected Structure or Process 1 superterm Name` = top_structures,
       fill = list(n = 0)
     ) %>%
     mutate(
       structure = factor(
         str_wrap(`Affected Structure or Process 1 superterm Name`, width = 20),
         levels = rev(str_wrap(top_structures, width = 20))
       ),
       gene = factor(`Gene Symbol`, levels = top_genes),
       n_plot = ifelse(n == 0, NA, n)
     )
   
   p1 <- ggplot(gs_plot, aes(x = gene, y = structure, fill = n_plot)) +
     geom_tile(color = "white", linewidth = 0.25) +
     scale_fill_viridis_c(na.value = "#eeeeee") +
     labs(x = "Gene", y = "Affected structure", fill = "Count") +
     theme_minimal() +
     theme(
       panel.grid = element_blank(),
       axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
   p1 
   
 }

 
 pglmm_copy_bar_plot<-function(col){
   orth %>% filter(sp %in% no_WGD_tree$tip.label) %>% 
     filter(Orthogroup %in% 
              (cafe_sig %>% filter(mod== col) %>% pull(genes))) %>% 
     dplyr::count(Orthogroup, sp)->pglmm_orths_count
   
   danio_orth_genes_modules %>% 
     left_join(Danio_biomart %>% dplyr::select(`Gene stable ID version`,`Gene name`), 
               by=c("orths"="Gene stable ID version"))->orths_gene_name
   
   left_join(pglmm_orths_count,orths_gene_name, by=c("Orthogroup"="gene")) %>% 
     complete(Orthogroup=Orthogroup, sp=sp) %>% 
     replace_na(list(n=0))%>% 
     distinct(Orthogroup,sp,.keep_all = TRUE) %>% 
     dplyr::select(Orthogroup, sp, n) %>% 
     left_join(orths_gene_name, by=c("Orthogroup"="gene")) %>% 
     distinct(Orthogroup,sp,.keep_all = TRUE) %>% 
     mutate(sp=factor(sp, levels = rev(no_WGD_tree$tip.label)),
            `Gene name`=reorder(`Gene name`, -n)) %>% 
     mutate(cn_group = case_when(
       n == 0~"0",
       n == 1~"1",
              n == 2~"2",
              n == 3~"3",
              n == 4~"4",
              n == 5~"5",
              n  > 5~">5"),
            cn_group = factor(cn_group, levels = c(as.character(0:5), ">5")))->heat_dat_pglmm_cp
   
  green_palette <- colorRampPalette(brewer.pal(9, "YlGnBu")[3:8])(5)
  highlight   <- "#D55E00"  #orange for >5
  fill_values <- c("white",green_palette, highlight)# this is for 0, 1-5, and >5 gene copies
  names(fill_values) <- levels(heat_dat_pglmm_cp$cn_group)
 
  heat_dat_pglmm_cp %>% 
   ggplot(aes(y=sp,x=`Gene name`, fill=factor(cn_group)))+
   geom_tile(alpha=0.8, color="black", linewidth = 0.1)+
    scale_fill_manual(values = fill_values,drop=FALSE,
                      name   = "Copy Number")+
   theme_bw()+
   theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1)) %>% return()
 }
 
 bgee_cafe_genes_plots <- function(orths_in, system, slice_size,ncol,nrow) {
   library(stringr)
   library(tidytext)
   library(ggh4x)
   library(randomcoloR)
   
   # Predefined anatomical-entity colors
   anatomical_cols <- c(
     "ovary" = "#355C7D",
     "muscle tissue" = "#E07A5F",
     "mesonephros" = "#D9A441",
     "heart" = "#C06C84",
     "intestine" = "#6C5B7B",
     "blastula" = "#7FB069",
     "camera-type eye" = "#4E9F8A",
     "tail" = "#F2CC8F",
     "testis" = "#8D99AE",
     "head kidney" = "#A44A3F",
     "nose" = "#5B8E7D",
     "bone element" = "#E6B8A2",
     "caudal fin" = "#8C6F9C",
     "zone of skin" = "#7A9E9F",
     "liver" = "#B86F52",
     "brain" = "#A7C957",
     "animal zygote" = "#B8A1D9",
     "adult organism" = "#C9ADA7",
     "pharyngeal gill" = "#6D597A",
     "spleen" = "#90BE6D",
     "actinopterygian pyloric caecum" = "#4D908E"
   )
   
   genes_to_use <- bgee_orth_res_to_use %>%
     filter(Orthogroup %in% orths_in) %>%
     pull(geneId)
   
   # Make plot data
   plot_dat <- bgee_combined_long %>%
     filter(Gene.ID %in% genes_to_use) %>%
     filter(!Anatomical.entity.name %in% tissues_to_remove)  %>%
     group_by(sp, Anatomical.entity.name) %>%
     summarise(m_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
     mutate(m_log10_TPM = log10(m_TPM)) %>%
     group_by(sp) %>%
     slice_max(order_by = m_log10_TPM, n = slice_size, with_ties = FALSE) %>%
     ungroup() %>%
     mutate(sp = str_replace(sp, "_", " "))
   
   # All tissues present in this plot
   tissue_levels <- plot_dat %>%
     distinct(Anatomical.entity.name) %>%
     arrange(Anatomical.entity.name) %>%
     pull(Anatomical.entity.name) %>%
     as.character()
   
   # Keep predefined colors for matching tissues
   matched_cols <- anatomical_cols[names(anatomical_cols) %in% tissue_levels]
   
   # Assign random colors to tissues not in the predefined palette
   missing_tissues <- setdiff(tissue_levels, names(anatomical_cols))
   
   if (length(missing_tissues) > 0) {
     set.seed(123)
     random_cols <- randomcoloR::distinctColorPalette(length(missing_tissues))
     names(random_cols) <- missing_tissues
   } else {
     random_cols <- character(0)
   }
   
   # Final named color vector
   my_colors <- c(matched_cols, random_cols)
   
   # Reorder tissues within each species facet, keeping original names for fill
   plot_dat <- plot_dat %>%
     mutate(
       tissue_ord = tidytext::reorder_within(
         Anatomical.entity.name,
         m_log10_TPM,
         sp
       ),
       Anatomical.entity.name = factor(Anatomical.entity.name, levels = tissue_levels)
     )
   
   plot_dat %>%
     filter(sp != "Poecilia reticulata") %>%
     ggplot(aes(y = tissue_ord, x = m_log10_TPM+1, fill = Anatomical.entity.name)) +
     ggh4x::facet_wrap2(
       ~ sp,
       scales = "free",
       strip.position = "right",
       ncol = ncol,
       nrow = nrow
     ) +
     geom_col(color = "black", linewidth = 0.2, alpha = 0.9) +
     tidytext::scale_y_reordered() +
     scale_fill_manual(values = my_colors, drop = FALSE) +
     labs(
       subtitle = system,
       y = NULL,
       x = "Mean expression (log10 TPM)"
     ) +
     theme_bw() +
     theme(
       legend.position = "none",
       strip.background = element_rect(fill = "grey96"),
       strip.text = element_text(size = 7)
     )
 }
 