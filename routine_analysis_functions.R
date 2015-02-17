
#############
# FUNCTIONS #
#############

load_taxon_table = function(tab_fp, map_fp, filter_cat, filter_vals, keep_vals){
  require(tools)
  require(biom)
  # load data
  if(file_ext(tab_fp) == 'biom'){
    data = read_biom(tab_fp)
    data = as.data.frame(as.matrix(biom_data(data)))
  }
  else if(file_ext(tab_fp) == 'txt'){
    data = read.table(tab_fp,sep='\t',skip=1,comment.char='',header=T,check.names=F,row.names=1)
    if(names(data)[ncol(data)] == 'taxonomy'){
      data$taxonomy = NULL
    }
  }
  else stop('Input file must be either biom (.biom) or tab-delimited (.txt) format.')
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  # optionally, subset data
    # cant subset if trying to filter out certain values and keep certain values
    # use one or the other
  if(!missing(filter_vals) & !missing(keep_vals)){
  }
  # filter out values within mapping category
  else if(!missing(filter_cat) & !missing(filter_vals)){
    map.f = map[!map[,filter_cat] %in% filter_vals, ]
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map.f = map[map[,filter_cat] %in% keep_vals, ]
  }
  else map.f = map
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(names(data),row.names(map.f))
  data.use = data[,match(samplesToUse,names(data))]
  data.use = data.use[rowSums(data.use)!=0,]
  map.use = map.f[match(samplesToUse,row.names(map.f)),]
  # output
  list(data_loaded = data.use, map_loaded = map.use)
}


load_ts_table = function(tab_fp, map_fp, filter_cat, filter_vals, keep_vals){
  require(tools)
  require(biom)
  # load data
  if(file_ext(tab_fp) == 'biom'){
    data = read_biom(tab_fp)
    data = as.data.frame(as.matrix(biom_data(data)))
  }
  else if(file_ext(tab_fp) == 'txt'){
    data = read.table(tab_fp, header=TRUE, sep="\t", row.names=1, comment.char="", check.names=FALSE)
  }
  else stop('Input file must be either biom (.biom) or tab-delimited (.txt) format.')
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  # optionally, subset data
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_vals) & !missing(keep_vals)){
  }
  # filter out values within mapping category
  else if(!missing(filter_cat) & !missing(filter_vals)){
    map.f = map[!map[,filter_cat] %in% filter_vals, ]
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map.f = map[map[,filter_cat] %in% keep_vals, ]
  }
  else map.f = map
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(names(data),row.names(map.f))
  data.use = data[,match(samplesToUse,names(data))]
  data.use = data.use[rowSums(data.use)!=0,]
  map.use = map.f[match(samplesToUse,row.names(map.f)),]
  # output
  list(data_loaded = data.use, map_loaded = map.use)
}


load_dm = function(dm_fp, map_fp, filter_cat, filter_vals, keep_vals){
  dm = read.table(dm_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  # optionally, subset data
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_vals) & !missing(keep_vals)){
  }
  # filter out values within mapping category
  else if(!missing(filter_cat) & !missing(filter_vals)){
    map.f = map[!map[,filter_cat] %in% filter_vals, ]
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map.f = map[map[,filter_cat] %in% keep_vals, ]
  }
  else map.f = map
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(names(dm), row.names(map.f))
  dm.use = as.dist(dm[match(samplesToUse,names(dm)), match(samplesToUse,names(dm))])
  map.use = map.f[match(samplesToUse,row.names(map.f)), ]
  # output
  list(dm_loaded = dm.use, map_loaded = map.use)
}

load_2_dms = function(dm1_fp, dm2_fp, map_fp, filter_cat, filter_vals, keep_vals){
  dm1 = read.table(dm1_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  dm2 = read.table(dm2_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  # optionally, subset data
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_vals) & !missing(keep_vals)){
  }
  # filter out values within mapping category
  else if(!missing(filter_cat) & !missing(filter_vals)){
    map.f = map[!map[,filter_cat] %in% filter_vals, ]
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map.f = map[map[,filter_cat] %in% keep_vals, ]
  }
  else map.f = map
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(intersect(names(dm1), row.names(map.f)), names(dm2))
  dm1.use = as.dist(dm1[match(samplesToUse,names(dm1)), match(samplesToUse,names(dm1))])
  dm2.use = as.dist(dm2[match(samplesToUse,names(dm2)), match(samplesToUse,names(dm2))])
  map.use = map.f[match(samplesToUse,row.names(map.f)), ]
  # output
  list(dm1_loaded = dm1.use, dm2_loaded = dm2.use, map_loaded = map.use)
}

filter_data = function(data, filter_cat, filter_vals, keep_vals){
  # input is list from 'load_data' function
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_vals) & !missing(keep_vals)){
  }
  # filter out values within mapping category
  else if(!missing(filter_cat) & !missing(filter_vals)){
    map.f = data$map_loaded[!data$map_loaded[,filter_cat] %in% filter_vals, ]
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map.f = data$map_loaded[data$map_loaded[,filter_cat] %in% keep_vals, ]
  }
  else map.f = data$map_loaded
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(names(data$data_loaded),row.names(map.f))
  data.use = data$data_loaded[,match(samplesToUse,names(data$data_loaded))]
  data.use = data.use[rowSums(data.use)!=0,]
  map.use = map.f[match(samplesToUse,row.names(map.f)),]
  # output
  list(data_loaded = data.use, map_loaded = map.use)
}

filter_taxa = function(table, filter_thresh, taxa_to_keep, taxa_to_remove){
  # filter taxa from otu table or taxa summary table based on mean abundance
  # optionally, specify additional taxa to keep
  means = apply(table[, 1:ncol(table)], 1, function(x){mean(x,na.rm=TRUE)})
  number_retained = sum((means >= filter_thresh) *1)
  taxa_keep = names(means[means >= filter_thresh])
  if(!missing(taxa_to_keep)) {taxa_keep = taxa_keep[taxa_keep %in% taxa_to_keep]}
  if(!missing(taxa_to_remove)) {taxa_keep = taxa_keep[!taxa_keep %in% taxa_to_remove]}
  table[row.names(table) %in% taxa_keep, ]
}
# table = arch_data$data_loaded
  
  
export_otu_table = function(tab, tax_fp, seq_fp, outfp){
  tax = read.table(tax_fp,sep='\t',comment.char='',header=F,check.names=F,row.names=1)
  seqs = read.table(seq_fp,sep='\t',comment.char='',header=F,check.names=F,row.names=1)
  otus = row.names(tab)
  tab.out = data.frame(tab, taxonomy = tax[match(otus, row.names(tax)),1],
                       sequence = seqs[match(otus, row.names(seqs)),1])
  write.table(tab.out, outfp, sep='\t', col.names=NA)
}

calc_dm = function(tab){
  require(vegan)
  # transform otu table (square root transformation)
  otuTable.xform = t(sqrt(tab))
  # create dissimilarity matrix from otu table
  otuTable.dist = vegdist(otuTable.xform, method='bray')
  otuTable.dist
}

plot_nmds = function(data, color_cat, shape_cat){
  require(ggplot2)
  if(missing(color_cat)) stop('Must include a mapping category to color by.')
  # format data and do NMDS
  dm = as.dist(data$dm_loaded)
  dm.mds = metaMDS(dm, k=2)
  # plot w shape
  if(!missing(shape_cat)){
    points = data.frame(dm.mds$points, cat = data$map_loaded[,color_cat], 
                        cat2 = data$map_loaded[,shape_cat])
    ggplot(points,aes(MDS1,MDS2,color=cat,shape=cat2)) +
      geom_point(size=3,alpha=.8) + theme_bw() +
      scale_color_discrete('') + scale_shape_discrete('') 
  }
  # plot without shape
  else{
    points = data.frame(dm.mds$points, cat = data$map_loaded[,color_cat])
    ggplot(points,aes(MDS1,MDS2,color=cat)) +
      geom_point(size=3,alpha=.8) + theme_bw() +
      scale_color_discrete('') + scale_shape_discrete('') 
  }
}

# functions to convert and manipulate dissimilarity matrices in 3 column format

convert_dm_to_3_column = function(dm){
  if(class(dm) == 'data.frame'){
    dmat = as.dist(dm)
  }
  else{
    dmat = dm
  }
  dmat.clmns = data.frame(t(combn(unlist(labels(dmat)),2)),as.numeric(dmat))
  names(dmat.clmns) = c('x1','x2','dist')
  dmat.clmns
}

add_metadata_to_df = function(dmat_clmns, map, cat){
  cat1 = map[match(dmat_clmns$x1,row.names(map)),cat]
  cat2 = map[match(dmat_clmns$x2,row.names(map)),cat]
  dmat_clmns_wCat = cbind(dmat_clmns, cat1, cat2)
  names(dmat_clmns_wCat) = c(names(dmat_clmns), paste(cat, "_1", sep=''), paste(cat, "_2", sep=''))
  dmat_clmns_wCat
}

cats_equal = function(x, col1, col2){
  if(x[col1] == x[col2]){"same"} else{"different"}
}

get_combination_category = function(x, accepted_categories){
  if(paste(x, collapse='__') %in% accepted_categories) {return(paste(x, collapse='__'))}
  else if(paste(rev(x), collapse='__') %in% accepted_categories) {return(paste(rev(x), collapse='__'))}
  else {return("Not accepted category")}
}

id_treatment_combination = function(col1, col2){
  # get the list of unique treatments
  unique_levels = unique(c(as.character(col1), as.character(col2)))
  # get all possible combinations
  combinations = rbind(t(combn(unique_levels, 2)), t(as.data.frame(lapply(unique_levels, FUN=rep, times=2))))
  combinations = paste(combinations[,1], combinations[,2], sep='__')
  # identify the combination for each pair of categories testing each order
  comparison_types = apply(data.frame(col1, col2), 1, get_combination_category, accepted_categories=combinations)
  comparison_types
}

get_mean_dissimilarities = function(data, comparison_cat, within_cat, comparison_cat_val){
  require(dplyr)
  dm_clmns = convert_dm_to_3_column(data$dm_loaded)
  # list x1 and x2 categories in new clmns
  dm_clmns_wCat = add_metadata_to_df(dm_clmns, data$map_loaded, within_cat)
  dm_clmns_wCat = add_metadata_to_df(dm_clmns_wCat, data$map_loaded, comparison_cat)
  # add column to indicate filter cat values (eg, only within sites)
  filter_cat_vals = apply(dm_clmns_wCat, 1, cats_equal, col1=4, col2=5)
  dm_clmns_wCat = cbind(dm_clmns_wCat, filter_cat_vals)
  # add column to say between what treatments
  comparison_types = id_treatment_combination(dm_clmns_wCat[ ,6], dm_clmns_wCat[ ,7])
  dm_clmns_wCat = cbind(dm_clmns_wCat, comparison_types)
  # calc mean dissimilarity of a treatment to control within site
  dissim_filt = filter(dm_clmns_wCat, filter_cat_vals=='same')
  means = summarize(group_by(dissim_filt, comparison_types, site_code_1), mean_dist = mean(dist))
  # filter to retain only categories of interest
  of_interest = comparison_cat_val
  means_filt = filter(means, comparison_types==of_interest)
  means_filt
}


