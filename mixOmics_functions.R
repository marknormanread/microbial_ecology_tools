library(matrixStats)


#### Filter low value ASVs out of the dataset. 
# Adapted from http://mixomics.org/mixmc/pre-processing/
# Function to perform pre-filtering. 
# To be retained, an ASV must comprise at least the given percentage of all reads found across all samples (pooled). 
low_count_removal = function(data,  # ASV count data frame of size n (sample) x p (OTU)
                             percent=0.01)  # Cutoff chosen
{
  keep_asv = which(colSums(data) * 100 / (sum(colSums(data))) > percent)
  data_filter = data[, keep_asv]
  return(list(data_filter = data_filter, keep_asv = keep_asv))
}


#### Total sum scaling (TSS) normalisation
# Each variable read count is divided by the total number of read counts. 
# Use on single columns at a time, NOT the whole table. 
# E.g. table_relab = t(apply(table, 1, TSS_normalisation)) where tables place features as columns. 
TSS_normalisation = function(x) { x/sum(x) }


#### Useful to run before choosing filter cutoff. 
# What is the distribution of percentages of total reads captured by the OTUs?
plot_feature_abundances = function(feature_table,  # Features as columns, samples as rows. Counts (not rel. abundance).
                                   xlab = 'Features, sorted by % total counts',
                                   ylab = 'Percent of total counts (%)',
                                   plot_log = TRUE,
                                   plot_ggplot = FALSE,
                                   filter_relative_abundance = -Inf,  # Exclude features with percentate less than this. 
                                   filter_n_abundant = Inf  # Include the top n.
                                   )
{
  feature_percent_all_counts = (100 * colSums(feature_table)) / sum(colSums(feature_table))
  feature_percent_all_counts = sort(feature_percent_all_counts, decreasing=TRUE)
  # Filter out ASV capturing less than this given relative abundance
  feature_percent_all_counts = feature_percent_all_counts[feature_percent_all_counts >= filter_relative_abundance]
  # Show only top n ASV by relative abundance. 
  feature_percent_all_counts = feature_percent_all_counts[1: min(length(feature_percent_all_counts), filter_n_abundant)]
  if (plot_ggplot)
  {
    df = data.frame(feature_percent_all_counts, rank=1:length(feature_percent_all_counts))
    p = ggplot(df, aes(x=rank, y=feature_percent_all_counts)) +
        geom_point() +
        ylab(ylab) +
        xlab(xlab)
    return(p)
  }
  else {
    plot(feature_percent_all_counts, xlab=xlab, ylab=ylab)
    if (plot_log)
      plot(feature_percent_all_counts, xlab=xlab, ylab=ylab, log='y')
  }
}

#### Transforms a raw feature-sample (as rows-cols) table into transposed relative abundance (total sum scaling).
# mixOmics adheres to traditional format of features as columns, however many microbial ecology bioinformatics 
# pipelines place them as rows. 
# Additionally, this prepares a pseudocount, needed for centred/isometric log ratio transforms, and filters. 
# Consider if filtering is important before performing it. 
# In targetted applications (looking at specific bugs or pathways) filtering may not be needed. 
# For exploratory analyses with many features it can be beneficial in disuading overfitting. 
prepare_mixomics = function(table,  # Supply as samples across columns and taxa/genes/predictors/features as rows. 
                            pseudocount = 1,  # Select this to be small w.r.t. the magnitude of data in `table`. 
                            filter_percent = 0.01  # Remove features representing < this value of all counts.
                            )
{
  table_mixomics = t(table)  # Transpose such that samples = rows (mixomics format, also more conventional)
  table_pseudocount = table_mixomics + pseudocount  # Can't log(0), which is used in CLR transform: add pseudocount.
  
  filter_result = low_count_removal(table_pseudocount, percent=filter_percent)
  table_filtered = filter_result$data_filter
  
  table_relab = t(apply(table_filtered, 1, TSS_normalisation))
  return(table_relab)
}


#### Genereric function for plotting a PCA/sPLS-DA ordination. 
ordinate = function(projection,  # A PCA or sPLS-DA
                    group,  # Factor vector descirbing the group of each sample included in `projection`.
                    plot_title,
                    filename = NULL,  # Location to save plot. Remember to append '.pdf'. 
                    name_samples = FALSE,  # If TRUE, write name of each sample on the plot.
                    sample_plot_characters = 16,  # Single char, else vector of n=len(samples). Used if name_samples=TRUE
                    ellipse = FALSE,
                    background = NULL,
                    distance = "mahalanobis.dist", 
                    components = c(1, 2),
                    col_per_group_map = NULL,  # Named vector. `group` names used as look up for colors. 
                    color_scheme = NULL,  # e.g., use 'rainbow'
                    show_legend = TRUE,
                    dimensions_mm = c(NA, NA),
                    fontsize = NULL,
                    title_fontsize = 18,  # Overriden by fontsize, if it is specified. 
                    spartan_plot = FALSE,  # Disable axis tick marks, ticks lines,  
                    point_size = 3)
{
  if (is.null(col_per_group_map))
  {
    if (length(levels(group)) > 10 | identical(color_scheme, 'rainbow'))
      col_per_group = rainbow(n=length(levels(group)))
    else
      col_per_group = color.mixo(1:length(levels(group)))
  } else {
    col_per_group = col_per_group_map[levels(group)]  # Extract the colours pertinent to the groups. 
  }

  if (! is.null(fontsize))
    title_fontsize = fontsize
  
  if (identical(background, TRUE) && identical(components, c(1, 2))) {
    # Not clear how to generate a background for anything but the first two components. 
    background = background.predict(projection, comp.predicted=2, dist=distance)
  }
  
  if (name_samples)
  {
    res = plotIndiv(projection, 
                    comp = components, # Select which components to plot on each axis. 
                    # This needs to be absent from the call to plot sample names. 
                    # pch = 16,  # Plot character to use
                    cex = point_size,
                    ind.names = name_samples,  # Do not name each item on the plot
                    group = group, 
                    col.per.group = col_per_group,
                    legend = show_legend,
                    ellipse = ellipse, 
                    background = background,
                    title = plot_title,
                    size.title = title_fontsize,
                    style = 'ggplot2'  # Set to return ggplot object as res$graph. Supposed to be default but doens't work. 
                    )
  } else {
    res = plotIndiv(projection,
                    comp = components, # Select which components to plot on each axis.
                    pch = sample_plot_characters,  # Plot character to use
                    cex = point_size,
                    ind.names = name_samples,  # Do not name each item on the plot
                    group = group,
                    col.per.group = col_per_group,
                    legend = show_legend,
                    ellipse = ellipse,
                    background = background,
                    title = plot_title,
                    size.title = title_fontsize,
                    style = 'ggplot2'  # Set to return ggplot object as res$graph. Supposed to be default but doens't work.
                    )
  }
  # MixOmics labels axes "X-variate". Much more meaningful to say "Component". 
  res$graph$labels$x = gsub('X-variate', 'Component', res$graph$labels$x)
  res$graph$labels$y = gsub('X-variate', 'Component', res$graph$labels$y)
  
  if (spartan_plot)
  {
    res$graph = res$graph + theme(legend.position="none",
                  text = element_text(size=fontsize),
                  axis.text = element_blank(),  # Turn off tick label titles.
                  axis.ticks = element_blank()  # Turn off black tick marks.
                  # plot.title = element_text(size=fontsize)
                  # Haven't been able to turn off the plot title completely. Stuck with the grey bar at top.
                  )
    
  }

  if (! is.null(filename)) {
    ggsave(filename, width=dimensions_mm[1], height=dimensions_mm[2], units='mm')
  }
  return(res$graph)
}


ordinate_1d = function(splsda_model, 
                       shape_vector = NULL,
                       col_per_group = NULL,
                       y_label = 'Groups\nseparated',
                       graph_path = NULL,  # Location on filesystem. Remember to provide file extension. 
                       width_mm = 60, height_mm = 30
                       )
{
  # if (is.null(shape_vector)) {
  #   shape_vector = rep(16, length(splsda_model$variates$X))
  # }
  
  data_df = data.frame(comp = splsda_model$variates$X, component=rep('1', length(splsda_model$variates$X)),
                       group_labels = splsda_model$Y, 
                       shape = shape_vector)
  names(data_df)[which(names(data_df) == 'comp.1')] = 'value'  # Rename column. 
  if (! is.null(shape_vector)) {
    # There is shape information; separate gtoups by this too. 
    p = ggplot(data_df, aes(x=value, y=component, color=group_labels, shape=shape)) +
      # cex adjusts point spacing within group dodge.width adjusts spacing between groups. 
      geom_beeswarm(cex = 5.0, groupOnX = FALSE, size = 1, dodge.width=1.3, aes(color=group_labels, shape=shape))
  } else {
    
    p = ggplot(data_df, aes(x=value, y=component, group=group_labels)) +
      geom_beeswarm(cex = 1.2, groupOnX = FALSE, size = 1, dodge.width=1.5, aes(color=group_labels, shape=group_labels))
  }
  p = p +
    theme(text = element_text(size=8)) +
    theme(legend.position="none") +
    ylab(y_label) +
    theme(axis.title.y = element_text(size=7)) +
    theme(axis.text.y = element_blank()) +  # Turn off tick label titles.
    xlab(paste('PC1:', round(100 * splsda_model$explained_variance$X), 'expl. var')) +
    theme(axis.title.x = element_text(size=7)) +
    theme(axis.text.x = element_blank()) +  # Turn off tick label titles.
    theme(axis.ticks = element_blank())  # Turn off black tick marks.
  if (! is.null(col_per_group)) {
    p = p + scale_color_manual(values=col_per_group)
  }

  if (! is.null(graph_path)) {
    ggsave(graph_path, width=width_mm, height=height_mm, units='mm')
  }
  return(p)
}


#### Gauge the statistical significance of a given real sPLS-DA model's performance by permutation testing. 
# The real sample-class assignments are permuted and the model is trained again. 
# Use the supplied arguments in exactly the same manner as the real model was built on (otherwise it is not a fair comparison)
statisitcal_significance_permutation_test = function(seed=123, repetitions=50,
                                                     feature_table,  
                                                     select_max_ncomp, 
                                                     select_distance,  # e.g. 'mahalanobis.dist'
                                                     select_test_keepX,  # c(<integer values>)
                                                     select_error_mode,  # 'BER' or 'overall'
                                                     select_validation,  # e.g. 'loo'
                                                     select_logratio = 'CLR',
                                                     select_cpus = 8, 
                                                     graph_path = NULL,
                                                     graph_title = '',
                                                     real_model_performance,
                                                     real_model_classes  # vector of classes used for the real model.
                                                     )
{
  set.seed(seed)  # Ensure reproducibility
  repetitions = repetitions
  lowest_errors = rep(NA, repetitions)
  for (repetition in 1:repetitions)
  {
    randomised_categories = factor(sample(real_model_classes))  # Shuffle the class labels assigned to our samples.
    tune_splsda = tune.splsda(X=feature_table, 
                              Y=randomised_categories,
                              ncomp=select_max_ncomp, 
                              logratio=select_logratio, 
                              test.keepX=select_test_keepX,  
                              validation=select_validation,
                              dist=c(select_distance),  
                              measure=select_error_mode,  # Our classes are unbalanced, this is a better measure.  
                              scale=TRUE,
                              near.zero.var=TRUE,
                              progressBar=FALSE,
                              cpus=select_cpus
    )
    lowest_errors[repetition] = min(tune_splsda$error.rate)
  }

  #### Make graph for publication
  # Paper is written in terms of accuracy, not error. 
  real_model_best_accuracy = 100 * (1 - tail(real_model_performance$error.rate[[select_error_mode]][, paste(select_distance, '.dist', sep='')], n=1))
  highest_accuracies = (1 - lowest_errors) * 100
  
  df = data.frame(accuracies = highest_accuracies)
  p = ggplot(df, aes(accuracies)) + stat_ecdf(geom = "step") + 
    # Represent where the real model's accuacy lay. 
    geom_vline(xintercept=real_model_best_accuracy, color='red') +
    ggtitle(graph_title) +
    xlab('Highest accuracies following tuning (%)') + 
    ylab('Cumulative distribution') +
    xlim(c(0, 100)) + ylim(c(0, 1)) +
    theme(text = element_text(size=8))
  if (! is.null(graph_path)) {
    ggsave(graph_path, width=7, height=5, units='cm')
  }
  
  beaten_by_random = sum(sort(highest_accuracies) >= real_model_best_accuracy)
  p_val_resolution = 1 / repetitions
  p_value_upper = (beaten_by_random + 1) * p_val_resolution
  p_value_lower = beaten_by_random * p_val_resolution
  cat('Highest accuracy from the real model was', real_model_best_accuracy, '%\n')
  cat('Random group re-assignments beat the performance of the real data (accuracy)', beaten_by_random, 'times\n')
  cat('This gives a p-value in range:', 
      p_value_upper, '< p <=', p_value_lower, 
      '(', repetitions, 'repetitions, giving p-resolution of', p_val_resolution, ')')

  return(invisible(list(plot=p, highest_accuracies=sort(highest_accuracies), lowest_errors=sort(lowest_errors))))
}


#### Wrapper that handles a select number of components for a splsda model. 
extract_important_features_all_components = function(
  model,  # an splsda model
  model_perf,  # Same `model``, but run through the `perf` function.
  max_components,  # How many components to explore. 
  csv_file_prefix,  # Specify directory and prefix (to identify model and groups explored) for csv files. 
  plot_loadings = TRUE,
  # 1D matrix, OTU/feature indexes as row names, the single column containing
  # a human-readable string (column named "name").
  feature_map = NULL)
{
  all_component_features = list()  # Populate in following loop.
  # Analyse each component in turn. Plot loadings, and output table that also describes stability in selection under LOO.
  for (component in 1:max_components) {
    contributions = plotLoadings(model, comp=component, 
                                 method='median',  # Advised for microbiome data
                                 contrib='max', 
                                 plot=FALSE  # Allows user to extract contribution matrix. Extraction disabled under TRUE. 
    )
    comp_features = extract_important_features(
      model_perf, contributions,
      component=component,
      filename=paste(csv_file_prefix, '_loadingsDim', component, '_selectedFeatures.csv', sep=''),
      feature_map=feature_map)
    all_component_features[[paste('component', component, sep='')]] = comp_features
    if (plot_loadings)
      pl = plotLoadings(model, comp=component, method='median',  # Advised for microbiome data
                        contrib='max')  
  }
  # How many features used in each component, and for each group. 
  for (component in 1:select_ncomp) {
    used_features = all_component_features[[component]][complete.cases(all_component_features[[component]]), ]
    cat('\nComponent', component, 'employed', dim(used_features)[1], 'features\n')
    groups = unique(used_features$GroupContrib)
    for (group in groups) {
      group_features = used_features[used_features$GroupContrib == group, ]
      cat(' ', dim(group_features)[1], 'taxa were associated with group', group, '\n')
    }
  }
  
  return(invisible(all_component_features))
}

#### Writes a table, captures features used in training a splsda model, but also queries the performance object 
# that evaluates it to capture how frequently each feature was used in training a model under cross validation. 
extract_important_features = function(splsda_perf, 
                                      # As generated by the plotLoadings function (w plot=FALSE)
                                      feature_component_contributions,  
                                      # 1D matrix, OTU/feature indexes as row names, the single column containing
                                      # a human-readable string (column named "name"). 
                                      feature_map = NULL, 
                                      component, filename)
{
  # How often (in cross validation) was each feature selected in model building?
  # Behaves differently if there is 1 vs >1 items in `splsda_perf$features$stable[[component]]`.
  # Try to handle both cases with the if statement conditional on length. 
  stability = as.data.frame(splsda_perf$features$stable[[component]])  
  if (length(splsda_perf$features$stable[[component]]) == 1)
    colnames(stability) = c('Freq')  # This isn't named correctly for length = 1
  else
    stability$Var1 = NULL  # Delete this column if length > 1
  rownames(stability) = names(splsda_perf$features$stable[[component]])
  
  # Adds new rows (containing NA values) for features used in model building, but not in the final model.
  feature_component_contributions = merge(feature_component_contributions, stability, 
                                          by.x='row.names', by.y='row.names', all=TRUE)
  rownames(feature_component_contributions) = feature_component_contributions$Row.names
  feature_component_contributions$Row.names = NULL
  
  # Rename columns 'importance' to 'loading'
  colnames(feature_component_contributions)[colnames(feature_component_contributions) == 'importance'] = 'loading'
  # feature_component_contributions$freq = freq
  if (! is.null(feature_map)) {
    # Put phylogenetic names in table. 
    feature_component_contributions$human_readable = feature_map[rownames(feature_component_contributions), 'name']
  }
  
  # Sort by frequency in cross validation loops, and then by absolute magnitude of loading. 
  feature_component_contributions = feature_component_contributions[
    order(feature_component_contributions$Freq, # Re-order rows
          abs(feature_component_contributions$loading),
          decreasing=TRUE)
    , ]  # All columns
  
  write.csv(feature_component_contributions, file=filename)
  return(feature_component_contributions)
}

#### Find the taxa (in terms of taxonomic names) that possess a given enzyme. Operates on PICRUSt2 data.
extract_taxa_possessing_enzyme = function(
  enzyme,  # E.g. "EC:3.2.1.26"
  enzyme_function='',  # E.g. "beta-fructofuranosidase"
  # Dataframe. PICRUSt2-generated file (e.g. EC_predicted.tsv). 
  # ASV sequences as rows, features (e.g. EC:3.2.1.26) as cols.
  ec_asv_tab, 
  # Dataframe. DADA2-generated file. First column is ASV sequences. Subsequent columns are:
  # "Kingdom,Phylum,Class,Order,Family,Genus,Species". Maps sequences to human readable taxonomy. 
  asv_taxonomy,
  filename=NULL,  # Where to save CSV
  return_table=FALSE
  )
{
  # Create a mask, extract rownames. Returns the ID, which is rRNA sequence. 
  asv_ids_possessing_ec = rownames(ec_asv_tab[ec_asv_tab[, enzyme] != 0, ])  
  asv_possessing_ec = asv_taxonomy[asv_ids_possessing_ec, ]  # Retrieve corresponding taxonomy. 
  asv_possessing_ec$Sequence_16S = rownames(asv_possessing_ec)  # Sequence as last column. 
  if (! is.null(filename))
    write.csv(asv_possessing_ec, file=filename, quote=FALSE, row.names=FALSE)
  cat('Found', length(asv_ids_possessing_ec), 'capable of producing', enzyme, '-', enzyme_function, 
      'out of', dim(ec_asv_tab)[1],'\n')
  if (return_table)
    return(asv_possessing_ec)
}

