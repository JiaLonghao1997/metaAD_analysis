library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)

options(stringsAsFactors = F)
inf <- '/share/home1/jialh/brain/metaAD2023/preprocessing/01_processing/readcounts.tsv'
outf <- '/share/home1/jialh/brain/metaAD2023/preprocessing/01_processing/readcounts.pdf'
# inf <- snakemake@input[[1]]
# outf <- snakemake@output[[1]]

readcounts <- read.table(inf, sep='\t', quote='', header=T)
readcounts$raw_count_frac <- 1.0
readcounts$raw_len_frac <- 1.0

readcounts <- readcounts[1:min(50, nrow(readcounts)), ]

print(colnames(readcounts))
readcounts.melt.count <- melt(readcounts[,c('Sample', 'raw_count',
                                            'dedup_count', 'trimmed_count',
                                            'host_removed_count')], id.vars = 'Sample')
readcounts.melt.count_frac <- melt(readcounts[,c('Sample', 'raw_count_frac',
                                            'dedup_count_frac', 'trimmed_count_frac',
                                            'host_removed_count_frac')], id.vars = 'Sample')


readcounts.melt.len <- melt(readcounts[,c('Sample', 'raw_len',
                                            'dedup_len', 'trimmed_len',
                                            'host_removed_len')], id.vars = 'Sample')
readcounts.melt.len_frac <- melt(readcounts[,c('Sample', 'raw_len_frac',
                                            'dedup_len_frac', 'trimmed_len_frac',
                                            'host_removed_len_frac')], id.vars = 'Sample')


g.count <- ggplot(readcounts.melt.count, aes(x=Sample, y=value, fill=variable)) + 
    geom_bar(stat='identity' , position='dodge') + 
    theme_bw() + 
    scale_fill_brewer(palette = 'Set2') +
    labs(title='Readcounts at each processing step', 
         y = 'Read count') + 
    guides(fill=guide_legend(title="Processing level")) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust=0.5))
# g.count

g.count_frac <- ggplot(readcounts.melt.count_frac, aes(x=Sample, y=value, fill=variable)) +
    geom_bar(stat='identity' , position='dodge') +
    geom_hline(aes(yintercept=0.5), color="red",size=2) +
    geom_hline(aes(yintercept=0.8), color="orange",size=2) +
    theme_bw() + 
    scale_fill_brewer(palette = 'Set2') +
    labs(title='Readcounts at each processing step', 
         y = 'Fraction of Raw ') + 
    guides(fill=guide_legend(title="Processing level")) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust=0.5))
# g.frac

####<--------------------------------------------------->####
g.len <- ggplot(readcounts.melt.len, aes(x=Sample, y=value, fill=variable)) +
    geom_bar(stat='identity' , position='dodge') +
    theme_bw() +
    scale_fill_brewer(palette = 'Set2') +
    labs(title='Read lens at each processing step',
         y = 'Read len') +
    guides(fill=guide_legend(title="Processing level")) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust=0.5))
# g.len

g.len_frac <- ggplot(readcounts.melt.len_frac, aes(x=Sample, y=value, fill=variable)) +
    geom_bar(stat='identity' , position='dodge') +
    geom_hline(aes(yintercept=0.5), color="red",size=2) +
    geom_hline(aes(yintercept=0.8), color="orange",size=2) +
    theme_bw() +
    scale_fill_brewer(palette = 'Set2') +
    labs(title='Read lens at each processing step',
         y = 'Fraction of Raw ') +
    guides(fill=guide_legend(title="Processing level")) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust=0.5))

# set plot width to be a function of the number of samples
nsamp <- nrow(readcounts)
plot.width <- 4 + (nsamp/4)
plot.height <- 6

pdf(outf, height = plot.height, width=plot.width)
print(g.count)
print(g.count_frac)
print(g.len)
print(g.len_frac)
dev.off()