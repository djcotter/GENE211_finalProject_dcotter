suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(stringr))

# parse arguments
option_list = list(
  make_option(c('-i', '--input'), type='character', default=NULL,
              help="path to input file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--pop'), type='character', default=NULL,
              help="Population to include in figure"),
  make_option(c('--width'), type='double', default=4.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=4.0,
              help='height for figure'),
  make_option(c('--units'), type='character', default='in',
              help='units for the figure', metavar="['in', 'cm', 'mm']"),
  make_option(c('--maxHeight'), type='double', default=1,
              help='max height for the plot (ylim).')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$units = tolower(opt$units)

if(is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Input and Output files must be specified.", call.=FALSE)
}

if(is.null(opt$pop)) {
  print_help(opt_parser)
  stop("A population must be specified.", call.=FALSE)
}

if( !(opt$units == 'in' || opt$units == 'mm' || opt$units == 'cm') ) {
  print_help(opt_parser)
  stop("Units must be provided as 'in', 'mm', or 'cm'.", call.=FALSE)
}

if( !(opt$width > 0 && opt$height > 0) ){
  print_help(opt_parser)
  stop("Width and Height must be a number greater than 0.")
}

# define function ------------------------------------------------------------------
pop_names <- function(x) {
  if (x %in% c("LWK", "GWD", "ACB", "ASW", "MSL", "YRI", "ESN")) {
    result <- 'AFR'
  } else if (x %in% c("CEU", "IBS", "TSI", "FIN", "GBR")) {
    result <- 'EUR'
  } else if (x %in% c("GIH", "BEB", "ITU", "STU", "PJL")) {
    result <- 'SAS'
  } else if (x %in% c("JPT", "CDX", "CHB", "CHS", "KHV")) {
    result <- 'EAS'
  } else if (x %in% c("PUR", "CLM", "MXL", "PEL")) {
    result <- 'AMR'
  } else {
    result <- NA
  }
  return(result)
}

# begin script --------------------------------------------------------------------------------------------------------

data <- read.delim(opt$input, header=F)
data <- dplyr::rename(data, subpop=V1, chr=V2, value=V3, val_low=V4, val_high=V5)
data$pop <- sapply(data$subpop, pop_names)
data$val_low <- sapply(data$val_low, function(x){if(x<=0){return(0)}else{return(x)}})

plot_data <- data %>% dplyr::filter(subpop == opt$pop)

p1 <- ggplot(plot_data, aes(x=reorder(chr,-value), y=value)) + 
  geom_errorbar(aes(ymin=val_low, ymax=val_high), color='black') + 
  geom_col(fill='black') + theme_pubr() + xlab('Genomic Region') +
  ylab('Average Heterozygosity') + 
  theme(legend.position = 'right')

ggsave(plot = p1, file=opt$output, height=opt$height, width=opt$width, units=opt$units)
