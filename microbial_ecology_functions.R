library(stats)
library(mixOmics)

#### How well separated are the samples from two different experimental groups?
# Two groups, A and B. 
# Calculates all pairwise distances between samples in A with themselves. Same with samples from B. 
# These compbine to give the 'within group' distribution. 
# Then calculate the pairwise distances between samples from A with those from B: 'between group' distribution. 
# These are plotted as probability density functions. 
# Kolmogorov smirnov test conducted on the distributions to gauge statistical significance of any separations. 
#
# Mark N. Read, 2019. 
within_between_group_distances = function(
  # Euclidean distance between CLR-transformed data recommended, as it copes well with compositional data. 
  distance_matrix,  # n x n matrix of distances between items. 
  group_vector,  # Eg column from metadata file. Order of items must match those in distance matrix.
  plotsave_prefix='',  # To save to a different directory, provide path through this prefix.
  # If TRUE, plot separate distributions for within groups A and B distances. Otherwise, combine them.
  separate_within_distributions=FALSE,
  titles=TRUE,
  xlimits=NULL  # or c(lower_limit, upper_limit)
)
{
  group_set = as.vector(unique(group_vector))  # Contrast only the groups in the supplied vector
  # Store pairwise group comparisons in terms of kolmogorov smirnov D and associated  p values in these nxn dataframes.
  group_differences_ks = data.frame(matrix(nrow=length(group_set), ncol=length(group_set)))
  group_differences_p = data.frame(matrix(nrow=length(group_set), ncol=length(group_set)))
  
  rownames(group_differences_ks) = group_set
  colnames(group_differences_ks) = group_set
  rownames(group_differences_p) = group_set
  colnames(group_differences_p) = group_set
  
  pairwise_combinations = combn(unique(group_set), 2)
  num_comparisons = ncol(pairwise_combinations)
  
  for (comparison in 1:num_comparisons)
  {
    group_a = pairwise_combinations[1, comparison]
    group_b = pairwise_combinations[2, comparison]
    
    a_vs_a = distance_matrix[group_vector == group_a, # Each of these produces a boolean vector.
                           group_vector == group_a]
    mask = lower.tri(a_vs_a)  # Retruns a boolean mask of dim(a_vs_a). Removes the diagonal too.
    a_vs_a_vect = a_vs_a[mask]  # Transforms to a vector.
    
    b_vs_b = distance_matrix[group_vector == group_b, group_vector == group_b]
    b_vs_b_vect = b_vs_b[lower.tri(b_vs_b)]  
    
    a_vs_b = distance_matrix[group_vector == group_a, group_vector == group_b]
    # DO NOT take the lower triangle only, all distances here need to be considered - no duplicate values. 
    between_group_distances = a_vs_b[TRUE]  # Transforms to a vector. 
    within_group_distances = c(a_vs_a_vect, b_vs_b_vect)  # Concatenate vectors
    
    ks = ks.test(x=within_group_distances, y=between_group_distances, 
                 alternative=c('greater'))  # Hypothesise that within (x) not greater than between (y)
    # Adjust for multiple comparisons. 
    p_value = p.adjust(ks$p.value, method='bonferroni', n=num_comparisons)
    
    if (separate_within_distributions == FALSE) {
      # ggplot needs data in dataframe format
      df = data.frame(Distance=c(within_group_distances, between_group_distances), 
                      Condition=c(rep('Within', length(within_group_distances)),
                                  rep('Between', length(between_group_distances))))
      p = ggplot(df, aes(x=Distance, fill=Condition)) + 
          # There are other options; this shows subtble differences in distribution.  
          # http://www.cookbook-r.com/Graphs/Plotting_distributions_(ggplot2)/
          geom_density(alpha=.3) +
          theme(plot.background = element_rect(fill = "transparent", color = NA)) + # Transparent plot background. 
          coord_cartesian(xlim = xlimits)
    } else {
      # ggplot needs data in dataframe format
      df = data.frame(Distance=c(within_group_distances, between_group_distances), 
                      Condition=c(rep(paste('Within ', group_a, sep=''), length(a_vs_a_vect)),
                                  rep(paste('Within ', group_b, sep=''), length(b_vs_b_vect)),
                                  rep('Between', length(between_group_distances))))
      p = ggplot(df, aes(x=Distance, color=Condition)) + geom_density(alpha=.3) +
          theme(plot.background = element_rect(fill = "transparent", color = NA))  # Transparent plot background. 
    }
    
    if (! is.null(xlimits)) {
      p + xlim(xlimits[1], xlimits[2])
    }
    if (titles) {
      p + ylab('Probability density') +
            xlab('Distance') +
            ggtitle(paste('Groups: ', group_a, ' and ', group_b,
                          # Include stats on title.
                          '; D=', round(ks$statistic, digits=2), ', p=', round(p_value, digits=2), sep=''))
      ggsave(paste(plotsave_prefix, 'group_distances_', group_a, '_vs_', group_b, '.pdf', sep=''))
    } else {
      # This completely removes ALL space and items around the polot itself. 
      p + theme(legend.position="none",
                axis.title = element_blank(),  # Turn off axis titles.
                axis.text = element_blank(),  # Turn off tick label titles.
                axis.ticks = element_blank(),  # Turn off black tick marks.
                text = element_blank(),
                panel.border = element_blank(),
                plot.margin=grid::unit(c(0,0,0,0),"cm"),
                axis.ticks.length = unit(0, "pt")
                )
      ggsave(paste(plotsave_prefix, 'group_distances_', group_a, '_vs_', group_b, '.pdf', sep=''),
             width=34, height=22, units='mm')
    }

    
    # Write to dataframe of results
    group_differences_ks[group_b, group_a] = ks$statistic
    group_differences_p[group_b, group_a] = p_value
  }
  return(list(group_differences_ks=group_differences_ks, 
              group_differences_p=group_differences_p))
}


#### Compositional data-robust distance calculations. 
# Conventional distance metrics (Bray Curtis, UniFrac, etc) have been criticised for introducing artefacts under 
# compisitional data, as microbiome sequencing data is. 
# Gloor et al advocate the 'Aitchison' distance metric instead.
# (https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full) 
# This is Euclidean distance in centred logratio transformed data.
aitchison_distance = function(taxa_matrix,  # Taxa/Features as columns, samples as rows
                              lr='CLR'  # Can also be 'ILR' for Isometric Logratio Transform. 
                              ) 
{
  taxa_clr = mixOmics::logratio.transfo(taxa_matrix,  # Samples as rows, Taxa/Features as columns. 
                                        logratio=lr, offset=0)
  
  clr_dist_raw = dist(taxa_clr, method='euclidean', diag=FALSE)  # Returns some weird 'dist' class ...
  # ... turn it into something better behaved, with functioning colnames and rownames. 
  clr_dist = as.matrix(clr_dist_raw) 
  return (clr_dist)
}
