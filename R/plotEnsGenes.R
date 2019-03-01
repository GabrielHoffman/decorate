
#' Get coordinates of exons
#' 
#' Get coordinates of exons from ENSEMBL database
#' 
#' @param ensdb ENSEMBL database object like EnsDb.Hsapiens.v86 
#' @param query GRranges ofject of one interval.  "chr20" should be coded as "20"
#' @param biotypes gene biotypes to return
#'
#' @return GRanges object of exon locations 
#' 
#' @examples
#' library(EnsDb.Hsapiens.v86)
#' library(GenomicRanges) 
#' 
#' # gene database
#' ensdb = EnsDb.Hsapiens.v86
#' 
#' # interval
#' query = GRanges("20", IRanges(62045027,62164563))
#' 
#' # get GRanges object of exon locations
#' get_exon_coords( ensdb, query)
#' 
#' @export
#' @importFrom GenomicFeatures exonsByOverlaps
#' @import GenomicRanges
get_exon_coords = function( ensdb, query, biotypes = c("protein_coding") ){

    if( ! is(ensdb, 'EnsDb') ){
        stop("ensdb must be an ENSEMBL databsed of type EnsDb")
    }

	# get gene info
	gr_exon = exonsByOverlaps( ensdb, query, columns=c( "exon_seq_start", "exon_seq_end", "symbol", 'gene_biotype', 'tx_seq_start', 'tx_seq_end', 'tx_cds_seq_start', 'tx_cds_seq_end'))# "gene_id", "exon_id",

    if( !is.na(biotypes) ){
    	# filter by biotype
    	gr_exon = gr_exon[gr_exon$gene_biotype %in% biotypes]
    }

	gr_exon
}


#' Plot ENSEMBL genes
#' 
#' Plot ENSEMBL genes in region
#' 
#' @param ensdb ENSEMBL database object like EnsDb.Hsapiens.v86 
#' @param minRange start genome coordinate
#' @param maxRange end genome coordinate
#' @param chromosome chrom
#' @param plot_lines_distance veritcal distance between genes
#' @param vp viewport
#' @param splice_variants if TRUE, show multiple transcripts from the same gene
#' @param non_coding if TRUE, also show non-coding genes
#'
#' @return GRanges object of exon locations 
#' 
#' @examples
#' library(EnsDb.Hsapiens.v86)
#' library(GenomicRanges) 
#' library(grid) 
#' 
#' # gene database
#' ensdb = EnsDb.Hsapiens.v86
#' 
#' # interval
#' query = GRanges("20", IRanges(62045027,62164563))
#'
#' # plot genes
#' fig = plotEnsGenes( ensdb, start(query), end(query), seqnames(query))
#' grid.draw( fig )
#' 
#' @export
#' @import GenomicRanges
#' @import grid
#' @importFrom data.table data.table
plotEnsGenes = function(ensdb, minRange, maxRange, chromosome, plot_lines_distance = 0.03,
    vp = viewport(x = 0, y = 0.95, just = c("left", "top")),
    splice_variants = TRUE, non_coding = TRUE){

    vp$xscale <- c(minRange, maxRange)
    vp$name <- "transcriptsVP"
    Range = maxRange - minRange
    map_len <- convertX(vp$width, "npc", valueOnly = TRUE)

	gr = GRanges(gsub( "^chr", "", chromosome), IRanges(minRange, maxRange))

    # get gene coordinates
    if( non_coding ){
        biotype = NA
    }else{
        biotype = c("protein_coding")
    }
	gr_exons = get_exon_coords( ensdb, gr, biotype )

    if( length(gr_exons) > 0){
   
    	df = data.table(data.frame(gr_exons))

        gene_biotype = 0
        symbol = 0
        tx_cds_seq_end = 0
        tx_cds_seq_start = 0
        tx_seq_end = 0
        tx_seq_start = 0

    	if( !splice_variants ){
    		# single body per gene
    		suppressWarnings(df_wide <- df[,data.frame(
    			gene_name = unique(symbol),
    			chrom = unique(seqnames),
    			strand = unique(strand), 
    			exonStarts = paste(start, collapse=','), 
    			exonEnds = paste(end, collapse=','),
    			exonCount = length(start),
    			biotype = unique(gene_biotype),
    			txStart = min(tx_seq_start, na.rm=TRUE),
    			txEnd = max(tx_seq_end, na.rm=TRUE),
    			cdsStart = as.integer(min(tx_cds_seq_start, na.rm=TRUE)),
    			cdsEnd = as.integer(max(tx_cds_seq_end, na.rm=TRUE)),
    			stringsAsFactors=FALSE),by=c("symbol")])
    	}else{

    		# multiple transcripts per gene		
            suppressWarnings(df_wide <- df[,data.frame(
    			gene_name = unique(symbol),
    			chrom = unique(seqnames),
    			strand = unique(strand), 
    			exonStarts = paste(start, collapse=','), 
    			exonEnds = paste(end, collapse=','),
    			exonCount = length(start),
    			biotype = unique(gene_biotype),
    			txStart = tx_seq_start,
    			txEnd = tx_seq_end,
                cdsStart = as.integer(min(tx_cds_seq_start, na.rm=TRUE)),
                cdsEnd = as.integer(max(tx_cds_seq_end, na.rm=TRUE)),
    			stringsAsFactors=FALSE),by=c("symbol", 'tx_seq_start', 'tx_seq_end')])
    	}
        t = data.frame(df_wide, stringsAsFactors=FALSE)
        t$plot_line <- 0
        t$plot_line[1] <- 1
    }else{
        t = data.frame()
    }
	
	plot_lines_no <- 1
	# browser()
    if (dim(t)[1] > 1) {
        for (i in 2:dim(t)[1]) {
            gene_name_width <- Range/map_len * convertWidth(grobWidth(textGrob(paste(t[1,"gene_name"], "   "), gp = gpar(fontsize = 7))),
                "npc", valueOnly = TRUE)
            for (j in 1:plot_lines_no) {
                if (max(t[1:(i - 1), ][t[, "plot_line"] == j, "txEnd"], na.rm = TRUE) < t[i, "txStart"] - gene_name_width*2) {
                  t[i, "plot_line"] <- j
                  break
                }
            }
            if (!t[i, "plot_line"])
                t[i, "plot_line"] <- plot_lines_no <- plot_lines_no + 1
        }
    }
    height = plot_lines_no * plot_lines_distance
    vp$height <- unit(height, "npc")
    # txt = "UCSC Genes Based on RefSeq, UniProt, GenBank, CCDS and Comparative Genomics"
    txt=''
    gene_plot_title <- textGrob(txt,
        gp = gpar(fontsize = 7, fontfamily = "mono"),
        just = c("centre", "center"), name = "gene_plot_title",
        default.units = "native")
    Transcripts <- gTree(children = gList(gene_plot_title), name = "transcripts",
        vp = vp)
    for(i in seq_len(nrow(t)) ){
        tx_upArrows <- tx_downArrows <- tx_exons <- tx_cds <- tx_leftArrow <- tx_rightArrow <- NULL
        y = 1 - t[i, "plot_line"]/plot_lines_no
        txStart <- max(minRange, t[i, "txStart"], na.rm=TRUE)
        txEnd <- min(maxRange, t[i, "txEnd"], na.rm=TRUE)
        cdsStart <- max(minRange, t[i, "cdsStart"], na.rm=TRUE)
        cdsEnd <- min(maxRange, t[i, "cdsEnd"], na.rm=TRUE)
        exonStarts <- as.numeric(strsplit(as.character(t[i, "exonStarts"]),
            ",")[[1]])
        exonStarts[exonStarts < minRange] <- minRange
        exonStarts[exonStarts > maxRange] <- maxRange
        exonEnds <- as.numeric(strsplit(as.character(t[i, "exonEnds"]),
            ",")[[1]])
        exonEnds[exonEnds < minRange] <- minRange
        exonEnds[exonEnds > maxRange] <- maxRange

        tx_name <- textGrob(x = txStart, y = y, paste(t[i, "gene_name"],
            " ", sep = ""), just = c("right", "center"), gp = gpar(fontsize = 7,
            fontfamily = "mono"), default.units = "native", name = "gene_name")
        tx_region <- linesGrob(x = c(txStart, txEnd), y = c(y,
            y), gp = gpar(lwd = 1, lineend = "butt"), default.units = "native",
            name = "tx_region")
        distance_between_arrows <- Range * (10/6)/convertX(unit(map_len,
            "npc"), "millimeters", valueOnly = TRUE)
        no_arrows <- floor(1/distance_between_arrows * (txEnd -
            txStart)) - 1
        if (no_arrows > 0) {
            arrow_heads_x0 <- unit(txStart + 1:no_arrows * distance_between_arrows,
                "native")
            arrow_heads_x1 <- unit(txStart + 1:no_arrows * distance_between_arrows,
                "native") + unit(0.5, "millimeters")
            if (t[i, "strand"] == "+") {
                change_y0 = unit(0.5, "millimeters")
                change_y1 <- unit(0, "native")
            }
            else {
                change_y0 <- unit(0, "native")
                change_y1 <- unit(0.5, "millimeters")
            }
            tx_upArrows <- segmentsGrob(x0 = arrow_heads_x0,
                x1 = arrow_heads_x1, y0 = unit(y, "native") +
                  change_y0, y1 = unit(y, "native") + change_y1,
                name = "tx_upArrows")
            tx_downArrows <- segmentsGrob(x0 = arrow_heads_x0,
                x1 = arrow_heads_x1, y0 = unit(y, "native") -
                  change_y0, y1 = unit(y, "native") - change_y1,
                name = "tx_downArrows")
        }
        tx_exons <- segmentsGrob(x0 = exonStarts, x1 = exonEnds,
            y0 = y, y1 = y, gp = gpar(lwd = 5, lineend = "butt"),
            default.units = "native", name = "tx_exons")
        cds0 <- cds1 <- NULL
        for (j in 1:length(exonStarts)) {
            if (cdsStart < exonEnds[j] && cdsEnd > exonStarts[j]) {
                cds0 <- c(cds0, max(exonStarts[j], cdsStart))
                cds1 <- c(cds1, min(exonEnds[j], cdsEnd))
            }
        }
        if (!is.null(cds0))
            tx_cds <- segmentsGrob(x0 = cds0, x1 = cds1, y0 = y,
                y1 = y, gp = gpar(lwd = 10, lineend = "butt"),
                default.units = "native", name = "tx_cds")
        if (t[i, "txStart"] < minRange)
            tx_leftArrow <- polygonGrob(x = unit.c(unit(minRange,
                "native") + unit(1, "millimeters"), unit(minRange,
                "native") + unit(2, "millimeters"), unit(minRange,
                "native") + unit(2, "millimeters"), unit(minRange,
                "native") + unit(2, "millimeters"), unit(minRange,
                "native") + unit(3, "millimeters"), unit(minRange,
                "native") + unit(3, "millimeters")), y = rep(unit.c(unit(y,
                "native"), unit(y, "native") + unit(1, "millimeters"),
                unit(y, "native") - unit(1, "millimeters")),
                2), rep(1:2, each = 3), gp = gpar(fill = "white"),
                name = "tx_leftArrow")
        if (t[i, "txEnd"] > maxRange)
            tx_rightArrow <- polygonGrob(x = unit.c(unit(maxRange,
                "native") - unit(1, "millimeters"), unit(maxRange,
                "native") - unit(2, "millimeters"), unit(maxRange,
                "native") - unit(2, "millimeters"), unit(maxRange,
                "native") - unit(2, "millimeters"), unit(maxRange,
                "native") - unit(3, "millimeters"), unit(maxRange,
                "native") - unit(3, "millimeters")), y = rep(unit.c(unit(y,
                "native"), unit(y, "native") + unit(1, "millimeters"),
                unit(y, "native") - unit(1, "millimeters")),
                2), rep(1:2, each = 3), gp = gpar(fill = "white"),
                name = "tx_rightArrow")
       
        gene_color <- "grey50"
        
        if( t[i,'biotype'] == "protein_coding"){
        	gene_color = rgb(12, 12, 120, maxColorValue = 255)
        }else if( t[i,'biotype'] == "miRNA"){
            gene_color = rgb(80, 80, 160, maxColorValue = 255)
        }else if( t[i,'biotype'] == "lincRNA"){
            gene_color = rgb(130, 130, 120, maxColorValue = 255)
        }else if( t[i,'biotype'] == "processed_transcript"){
            gene_color = rgb(255, 51, 255, maxColorValue = 255)
        }

        transcript <- gTree(children = gList(tx_name, tx_region,
            tx_upArrows, tx_downArrows, tx_exons, tx_cds, tx_leftArrow,
            tx_rightArrow), gp = gpar(col = gene_color))#, vp=vp)#, name = t[i,
            # "name"])
        # browser() 
        Transcripts <- addGrob(Transcripts, transcript)
    }
    attr(Transcripts, 'height') = height + (1-as.numeric(vp$y))
    return(Transcripts)
}