library(gTrack)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(Homo.sapiens)
library(zoo)
library(ggplot2)
library(ArchR)
source("data/custom_colors.R")

#--------------------------------------------
# Functions
#--------------------------------------------

hic_to_gTrack_matrix <- function(hic_file, intervals, res, chr_order, max_value, juicer_path, java_path) {
	options(scipen = 99)
	intervals_gr <- GRanges(seqnames = intervals$V1
		, ranges = IRanges(round_any(intervals$V2 ,res)
			, round_any(intervals$V3, res, f = ceiling)))
	output_files <- list()
	for (i in 1:length(intervals_gr)) {
		chr1 <- as.character(seqnames(intervals_gr))[i]
		start1 <- start(intervals_gr)[i]
		end1 <- end(intervals_gr)[i]
		for (j in 1:length(intervals_gr)) {
			chr2 <- as.character(seqnames(intervals_gr))[j]
			start2 <- start(intervals_gr)[j]
			end2 <- end(intervals_gr)[j]
			output <- paste0(gsub('.allValidPairs.hic', '', hic_file), '.',chr1,':',start1, ":",end1, "."
				, chr2,':',start2, ":",end2, ".", res,'.txt')
			if (file.exists(paste0(gsub('.allValidPairs.hic', '', hic_file), '.',chr2,':',start2, ":",end2, ".",
				chr1,':',start1, ":",end1, ".", res,'.txt'))) {
				next
			}
			system(paste0(java_path, ' -jar ', juicer_path,' dump observed ','KR ',hic_file, ' '
				, chr1, ':',start1, ":",end1,' ',chr2,':',start2, ":",end2,' BP ', res,' ',output))
			output_files[[paste0(i,"_",j)]] <- output
		}
	}
	
	matrix_all <- data.frame(matrix(nrow = 0, ncol = 0))
	for (i in 1:length(output_files)) {
		int1 <- as.numeric(gsub("_.*", "", names(output_files[i])))
		int2 <- as.numeric(gsub(".*_", "", names(output_files[i])))
		range1 <- seq(start(intervals_gr)[int1], end(intervals_gr)[int1], by = res)
		range2 <- seq(start(intervals_gr)[int2], end(intervals_gr)[int2], by = res)
		matrix <- read.delim(output_files[[i]], header = FALSE, stringsAsFactors = FALSE)
		matrix$V1 <- factor(matrix$V1, levels = range1)
		matrix$V2 <- factor(matrix$V2, levels = range2)
		matrix <- matrix %>% complete(V1, V2, fill = list(V3 = 0))
		matrix <- matrix[order(as.numeric(matrix$V1), as.numeric(matrix$V2)), 1:3]
		matrix$V1 <- paste(as.character(seqnames(intervals_gr))[int1], matrix$V1, sep ="_")
		matrix$V2 <- paste(as.character(seqnames(intervals_gr))[int2], matrix$V2, sep ="_")
		rownames(matrix) <- NULL
		matrix_all <- rbind(matrix_all, matrix)
		system(paste0("rm ", output_files[i]))
	}
	matrix_all$chr1 <- factor(gsub("_.*", "", matrix_all$V1), levels = chr_order)
	matrix_all$chr2 <- factor(gsub("_.*", "", matrix_all$V2), levels = chr_order)
	matrix_all <- matrix_all[order(matrix_all$chr1, matrix_all$chr2),1:3]
	matrix_all_spread <- data.frame(spread(matrix_all, key = V2, value = V3, fill = 0))
	rownames(matrix_all_spread) <- matrix_all_spread$V1
	matrix_all_spread <- data.frame(matrix_all_spread[,c(2:length(colnames(matrix_all_spread)))])
	matrix_all_spread <- matrix_all_spread[unique(matrix_all$V1),unique(matrix_all$V2)]
	matrix_all_spread[matrix_all_spread > max_value] <- max_value
	matrix_all_mat <- as.matrix(matrix_all_spread)
	m <- dim(matrix_all_mat)[1]
	for (i in 1:m) {
		for (j in 1:m) {
			n <- max(matrix_all_mat[i,j], matrix_all_mat[j,i], na.rm = TRUE)
			matrix_all_mat[i,j] <- n
			matrix_all_mat[j,i] <- n
		}
	}
	rownames(matrix_all_mat) <- NULL
	colnames(matrix_all_mat) <- NULL
	matrix_all_mat
}

# from ArchR package (Granja, Corces et al. Nature Genetics 2021)
getArchDF <- function(lp, r = 100){
	angles <- seq(pi, 2*pi,length.out=100)
	rx <- (end(lp)-start(lp))/2
	rscale <- r * (rx/max(rx))
	cx <- start(lp) + rx
	if(is.null(mcols(lp)$value)){
		mcols(lp)$value <- 1
	}
	df <- lapply(seq_along(cx), function(z){
		xz <- rx[z]*cos(angles)+cx[z]
		dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), value = mcols(lp)$value[z])
	}) %>% Reduce("rbind",.)
	return(df)
}

getPromotersOfGene <- function(gene, upstream = 100, downstream = 100, genome = 'hg38'){
	require(org.Hs.eg.db)
	if (genome == 'hg38'){
		require(TxDb.Hsapiens.UCSC.hg38.knownGene)
		txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	} else if (genome == 'hg19'){
		require(TxDb.Hsapiens.UCSC.hg19.knownGene)
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
	} else { stop('genome does not exist') }
	
	promoters <- suppressWarnings( promoters(txdb, upstream = upstream, downstream = downstream) )
	promoters <- GenomeInfoDb::keepStandardChromosomes(promoters, pruning.mode = "coarse")
	eid <- suppressMessages( AnnotationDbi::select(org.Hs.eg.db, gene, 'ENTREZID', 'SYMBOL') )
	transcript_id <- suppressMessages( AnnotationDbi::select(txdb, eid$ENTREZID, "TXNAME", "GENEID") )
	gene_key <- merge(eid, transcript_id, by.x = 'ENTREZID', by.y = 'GENEID')
	promoters_of_gene <- promoters %>% subset(tx_name %in% gene_key$TXNAME)
	promoters_of_gene$symbol <- gene_key[match(promoters_of_gene$tx_name, gene_key$TXNAME),'SYMBOL']
	return(promoters_of_gene)   
}

plot_v4c <- function(hic_files, names, norm, bedpe, window, anchor_gene, color
		, gene_list, ymax,res = 10000, juicer_path, java_path) {
	options(scipen = 999)
	loops <- read.delim(bedpe)
	chr <- gsub(":.*", "", window)
	chrstart_plot <- as.numeric(gsub("-.*", "", gsub(".*:", "", window)))
	chrstart <- as.numeric(gsub("-.*", "", gsub(".*:", "", window))) - res*10
	chrend_plot <- as.numeric(gsub(".*-", "", window))
	chrend <- as.numeric(gsub(".*-", "", window)) + res *10
	prmtrs <- getPromotersOfGene(anchor_gene, genome = "hg19")
	anchor_chr <- unique(as.character(seqnames(prmtrs)))
	anchor <- min(round((end(prmtrs) + start(prmtrs))/2))
	region <- GRanges(seqnames = chr, ranges = IRanges(start = chrstart_plot, end = chrend_plot))
	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
	species <- Homo.sapiens
	genes <- suppressMessages(data.frame(genes(txdb)))
	genenames <- suppressMessages(AnnotationDbi::select(species, keys = unique(genes$gene_id), columns = 'SYMBOL', keytype = 'GENEID'))
	genes <- merge(genes, genenames, by.x = "gene_id", by.y = "GENEID")
	geneinfobed <- genes[,c("seqnames", "start", "end", "SYMBOL", "strand")]
	geneinfobed <- geneinfobed[geneinfobed$SYMBOL %in% gene_list,]
	colnames(geneinfobed) <- c("chrom", "start", "stop", "gene", "strand")
	geneinfobed$score <- "."
	geneinfobed$pos <- seq(1:length(rownames(geneinfobed)))
	geneinfobed$labels <- "Genes"
	vp_plot <- data.frame(matrix(nrow = 0, ncol = 4))

	for (i in 1:length(hic_files)) {
		hic_file <- hic_files[i]
		dump_file <- paste0("~/tmp/", basename(hic_file), "_", chr,"_",res, ".dump")
		system(command = paste0(java_path, " -jar ", juicer_path," dump observed NONE ", hic_file, " "
			, anchor_chr, ":", round_any(anchor, res, f = floor),":"
			,round_any(anchor, res, f = ceiling)," ", chr, " BP ", res, " ", dump_file))
		chr_norm <- suppressMessages(readr::read_delim(dump_file,delim = "\t", col_names = FALSE))
		anchor_bin <- unique(chr_norm$X1[anchor >= chr_norm$X1 & anchor < chr_norm$X1 + res])
		if( length(anchor_bin) == 0 ) {
			anchor_bin <- unique(chr_norm$X2[anchor >= chr_norm$X2 & anchor < chr_norm$X2 + res])
		}
		loops_plot <- loops[loops$chr1 == anchor_chr & (abs(rowMeans(loops[,c("x1", 'x2')]) - anchor) < 25000 |
			abs(rowMeans(loops[,c("y1", 'y2')]) - anchor) < 25000),]
		loops_plot$Q.Value_Bias[loops_plot$Q.Value_Bias == 0] <- min(loops_plot$Q.Value_Bias[
			loops_plot$Q.Value_Bias != 0], na.rm = TRUE)
		loops_plot <- loops_plot[order(loops_plot$Q.Value_Bias, decreasing = TRUE),]
		loops_gr <- GRanges(seqnames = loops_plot$chr1
			, ranges = IRanges(start = rowMeans(loops_plot[,c("x1", 'x2')])
			, end = rowMeans(loops_plot[,c("y1", 'y2')])))
		loops_gr$value <- -log10(loops_plot$Q.Value_Bias)
		loops_list <- SimpleList(Loops = loops_gr)
		
		loopO <- lapply(seq_along(loops_list), function(x){
			subLoops <- subsetByOverlaps(loops_list[[x]]
				, region, ignore.strand = TRUE, type = "within") 
			if(length(subLoops)>0){
				dfx <- getArchDF(subLoops)
				dfx$name <- Rle(paste0(names(loops_list)[x]))
				dfx
			}else{
				NULL
			}})
		chr_norm_vp <- chr_norm[(chr_norm$X1 == anchor_bin & chr_norm$X2 > chrstart & chr_norm$X2 < chrend) ,]
		tmp <- chr_norm[(chr_norm$X2 == anchor_bin & chr_norm$X1 > chrstart & chr_norm$X1 < chrend) ,]
		colnames(tmp) <- c("X2", "X1", "X3")
		chr_norm_vp <- rbind(chr_norm_vp, tmp)
		chr_norm_vp <- chr_norm_vp[,c("X1", "X2","X3")]
		if (!(chrstart %in% chr_norm_vp$X2)) {
			chr_norm_vp[nrow(chr_norm_vp) + 1,] = list(anchor_bin,chrstart, 0)
		}
		
		if (!(chrend %in% chr_norm_vp$X2)) {
			chr_norm_vp[nrow(chr_norm_vp) + 1,] = list(anchor_bin,chrend, 0)
		}
		chr_norm_vp <- chr_norm_vp[order(chr_norm_vp$X2),]
		chr_norm_vp_plot <- chr_norm_vp
		chr_norm_vp_plot$sample <- names[i]
		chr_norm_vp_plot$X3 <- chr_norm_vp_plot$X3/ (norm[i]) *100000
		vp_plot <- rbind(vp_plot, chr_norm_vp_plot)
	}
	vp_plot$sample <- factor(vp_plot$sample, levels = names)
	vp_plot$X3 <- rollmean(vp_plot$X3, k = 3, fill = 0)
	if (is.null(ymax)) {
		ymax <- max(vp_plot$X3[abs(vp_plot$X2 - vp_plot$X1) > res*5])
	}
	p1 <- ggplot() + 
		geom_vline(xintercept = anchor, color = "grey", linetype = "longdash", size = 1) + 
		geom_line(data = vp_plot, aes(x = X2, y = X3, color = sample)) + theme_classic() + 
		coord_cartesian(ylim=c(0, ymax), xlim = c(chrstart_plot, chrend_plot)) +
		scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
		ggtitle(label = paste0(chr, ":", chrstart_plot, "-", chrend_plot, "; ", res /1000, " kb resolution")) + 
		theme(legend.position = "top", axis.title.x = element_blank(), axis.line = element_blank(),
			panel.border = element_rect(colour = "black", fill=NA, size = 1)) + ylab("Norm EIS") + 
			scale_y_continuous(breaks = scales::pretty_breaks(n = 2))
	if (is.null(color)) {
		p1 <- p1 + scale_color_manual(values = paletteDiscrete(values = levels(vp_plot$sample)))
	} else {
		p1 <- p1 + scale_color_manual(values = color)
	}
	p1 <- p1 + facet_grid(sample ~ .) + theme(legend.position = "none",  strip.background = element_blank())
	if (is.null(color)) {
		pal <- colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7","black"))(100)
	} else {
		pal <- colorRampPalette(c("white",color))(100)
	}
	if (!is.null(loopO[[1]])) {
		valueMin <- min(loopO[[1]]$value)
		valueMax <- max(loopO[[1]]$value)
		p2 <- ggplot(data = data.frame(loopO), aes(x = x, y = y, group = id, color = value)) + 
			geom_line() +
			facet_grid(name ~ .) + theme_classic() +
			ylab("") + xlab("") + 
			coord_cartesian(ylim = c(-100,0)) +
			scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
			scale_color_gradientn(colors = pal, limits = c(valueMin, valueMax), name = "-log10(Adjusted P-value)") +
			theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.line = element_blank()
				, strip.background = element_blank(), axis.text.x = element_blank(), 
				, axis.ticks.x = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size = 1)
				, legend.box.background = element_rect(color = NA)) +
			guides(color= guide_colorbar(barwidth = 0.75, barheight = 3)) 
		} else {
			df <- data.frame(facet = "Loops", start = 0, end = 0, strand = "*", symbol = "none")
			p2 <- ggplot(data = df, aes(start, end)) + 
				geom_point() +
				facet_grid(facet~.) + theme_classic() +
				scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
				theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
				theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()
					, axis.line = element_blank()
					, strip.background = element_blank(), axis.text.x = element_blank()
					, axis.ticks.x = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size = 1))
				}
	
	p3 <- ggplot(data = geneinfobed) + 
		geom_segment(data = geneinfobed, aes(x = start, xend = stop, y = pos, yend = pos, color = strand),size = 2) + 
		geom_label(data = geneinfobed, aes(x = start, y = pos, label = gene, color = strand),nudge_y = -0.5) +
		facet_grid(labels ~ .) + theme_classic() +
		scale_color_manual(values = paletteDiscrete(values = c("+","-")), guide = "none") + 
		ylab("") + xlab("") + 
		ylim(min(geneinfobed$pos) -1, max(geneinfobed$pos) +1) + 
		scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) + 
		theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank()
			, axis.ticks.x = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size = 1)
			, strip.background = element_blank(), axis.line = element_blank())
			
	height <- 1
	if (length(levels(vp_plot$sample)) > 3) {
		height <- min(1*length(levels(vp_plot$sample)), 6) 
	}
	plot_return <- patchwork::wrap_plots(p1,p2,p3, ncol =1, heights = c(height,0.5, 0.25))
	plot_return
}

#--------------------------------------------
# Plotting
#--------------------------------------------

# Get genes for gTrack plotting
gt.ge = track.gencode(genes =c("MYC", "PVT1", "CASC8", "FAM84B", "DUSP22", "PLUT", "PCAT1")
	, bg.col = "black", cds.col = "black", utr.col = "black"
	, st.col = "black", en.col = "black", height = 50)

# Read in genomic intervals and specify plotting windows for gTrack
intervals <- read.delim("data/COLO320DM_intervals.txt", header = FALSE)
intervals_gr_vis <- GRanges(seqnames = intervals$V1
	, ranges = IRanges(intervals$V2, intervals$V3+res-1))
intervals_gr_vis_windows <- unlist(slidingWindows(intervals_gr_vis, width = res, step = res))
seqlevels(intervals_gr_vis) <- chr_order
seqlevels(intervals_gr_vis_windows) <- chr_order
intervals_gr_vis <- sort(intervals_gr_vis)
intervals_gr_vis_windows <- sort(intervals_gr_vis_windows)

# Dump matrix from hic file and reformat for gTrack
hic_file <- 'data/COLO320DM_K27ac_HiChIP.allValidPairs.hic' # available on GEO at GSE159985
res <- 10000 
chr_order <- c("chr6", "chr16", "chr13", "chr8")
max_value <- 4
juicer_path <- '~/tools/juicer/scripts/common/juicer_tools.jar'
java_path <- '~/anaconda3/bin/java'

matrix_all_mat <- hic_to_gTrack_matrix(hic_file, intervals, res, chr_order, max_value, juicer_path, java_path)

# Plot HiChIP matrix
plot.gTrack(c(gt.ge, gTrack(intervals_gr_vis_windows, mdata = matrix_all_mat, stack.gap = 0.5
	, colormaps = colorRampPalette(c("white", brewer.pal(n = 9, name = "Reds")))(100), height = 500))
	, windows = intervals_gr_vis, legend.params = list(plot = FALSE))

# Read in ecDNA reconstruction and make graph for gTrack plotting
segments <- read.delim("data/COLO320DM_OM_segments_ordered.txt")
gr <- GRanges(seqnames = segments$chrom, ranges = IRanges(start = segments$start,end = segments$end), strand = segments$strand)
gr$color = as.character(seqnames(gr))
gr$y = c(1:length(gr))
colors <- c("#54AC68", "#483D8B", "#054907", "#8756E4")
names(colors) <- c("chr6", "chr8", "chr13", "chr16")
graph = data.frame(from = c(seq(1, length(gr)-1), 22), to = c(seq(2, length(gr)), 1))

# Plot ecDNA reconstruction and ChIP coverage tracks
# bigwig files available on GEO at GSE159972
plot.gTrack(c(gt.ge, gTrack('data/COLO320DM_WGS.bw', col = my_custom_palettes$splash_of_salmon[4]
			, yaxis.pretty = 1, ylab = "WGS",cex.label = 0.5,height = 100, bar = TRUE, y1 = 100)
		, gTrack('data/COLO320DM_Brd4.bw', col = my_custom_palettes$splash_of_salmon[4], yaxis.pretty = 1
			, ylab = "BRD4 ChIP-seq",cex.label = 0.5,height = 100, bar = TRUE, y1 = 500)
		, gTrack('data/COLO320DM_DMSO_H3K27ac.bw', col = my_custom_palettes$splash_of_salmon[4]
			, yaxis.pretty = 1, ylab = "H3K27ac ChIP-seq",height = 100, bar = TRUE, y1 = 500)
		, gTrack(gr, edges = graph , stack.gap = 5, y.field = "y",height = 500, yaxis = FALSE
			, gr.colorfield = "color", colormaps = colors, y1 = length(gr) + 1))
	, windows = intervals_gr_vis, legend.params = list(plot = FALSE))

# Plot virtual 4C signal at PVT1 promoter
names(hic_file) <- gsub("\\..*", "",gsub(".*data/", "", hic_file))
validPairs <- 19162671
names(validPairs) <- names(hic_file)
loops_file <- "data/COLO320DM_K27ac_HiChIP_FitHiChIP.interactions_Q0.1_MergeNearContacts.bedpe"
names(loops_file) <- names(hic_file)
mycolors <- c(my_custom_palettes$splash_of_salmon, "#E6C1BC")
gene_list <- c("MYC", "PVT1", "CASC8", "PCAT1")

plot_v4c(hic_file, color = mycolors[1], names = names(hic_file)[1], res = res
	, norm = validPairs, bedpe = loops_file, window = "chr8:127420000-129050000"
	, anchor_gene = "PVT1", gene_list = gene_list, ymax = 2, java_path = java_path
	, juicer_path = juicer_path)