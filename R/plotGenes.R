# Adapted from LDheatmap
#' Plot genes from a specified region of the human genome.
#'
#'  Retrieves genes from the UCSC Genome Browser and generate the genes plot.
#'
#' @usage plotGenes(minRange, maxRange, chromosome, genome = "hg19", plot_lines_distance = 0.03,
#'vp = viewport(x = 0, y = 0.99, just = c("left", "top")), splice_variants = TRUE,
#'non_coding = TRUE)
#' @param minRange The sequence minimum range in base pairs.
#' @param maxRange The sequence maximum range in base pairs.
#' @param chromosome A character string identifying the chromosome.
#' @param genome The genome assembly to use. The default is hg19, the most recent human genome assembly on
#'the UCSC genome browser.
#' @param plot_lines_distance The distance between the lines of genes plotted.
#' @param vp A \code{viewport}.
#' @param splice_variants If \code{FALSE}, exclude gene splice variants.
#' @param non_coding If \code{FALSE}, exclude non-coding genes.
#' @details The genes are color coded as follows:
#'Black -- feature has a corresponding entry in the Protein Data Bank (PDB)
#'Dark blue -- transcript has been reviewed or validated by either the RefSeq, SwissProt or CCDS staff
#'Medium blue -- other RefSeq transcripts
#'Light blue -- non-RefSeq transcripts
#'
#'For assemblies older than hg18, all genes are plotted in grey.
#' @return A \code{grob} of gene plots.
#' @references \url{http://genome.ucsc.edu/cgi-bin/hgTrackUi?g=knownGene}
#' @author Sigal Blay <sblay@sfu.ca> and more
#' @examples \dontrun{
#'grid.newpage()
#'plotGenes(149500000, 150000000, "chr7")
#'}
#' @import rtracklayer
#' @import GenomeInfoDb
#' @importFrom grDevices rgb
#' @importFrom IRanges IRanges
#
# ldheatmap - Plots measures of pairwise linkage disequilibria for SNPs
# Copyright (C) 2004  J.Shin, S. Blay, N. Lewin-Koh, J.Graham, B.McNeney
#
# To cite LDheatmap in publications use:
# Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap: An R
# Function for Graphical Display of Pairwise Linkage Disequilibria
# Between Single Nucleotide Polymorphisms. J Stat Soft, 16 Code Snippet 3
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
###########################################################################
plotGenes = function (minRange, maxRange, chromosome, genome = "hg19", plot_lines_distance = 0.03,
    vp = viewport(x = 0, y = 0.99, just = c("left", "top")),
    splice_variants = TRUE, non_coding = TRUE){
    requireNamespace("grid")
    map_len <- convertX(vp$width, "npc", valueOnly = TRUE)
    Range <- maxRange - minRange
    vp$xscale <- c(minRange, maxRange)
    vp$name <- "transcriptsVP"
    cat("Connecting to UCSC...\n")
    session <- browserSession()
    genome(session) <- genome
    cat("Connection extablished.  Processing query...\n")
    query1 <- ucscTableQuery(session, "knownGene",
        GRangesForUCSCGenome(genome, chromosome,
            IRanges(minRange, maxRange)))
    t <- getTable(query1)
    if (!dim(t)[1]) {
        print("The genetic region of the data does not correspond to any genes in the UCSC genome browser")
        return()
    }
    t[, "name"] <- as.character(t[, "name"])
    if (!non_coding) {
        ind <- NULL
        for (i in 1:dim(t)[1]) if (t[i, "cdsStart"] == t[i, "cdsEnd"])
            ind <- c(ind, i)
        if (!is.null(ind))
            t <- t[-ind, ]
    }
    if (!splice_variants) {
        query2 <- ucscTableQuery(session, "knownGene",
            GRangesForUCSCGenome(genome, chromosome,
                IRanges(minRange, maxRange)), table = "knownCanonical")
        tcanon <- getTable(query2)
        ind <- (t$name %in% as.character(tcanon$transcript))
        t <- subset(t, ind)
    }
    t[, "gene_name"] <- ""
    tbl <- "kgXref"
    query2 <- ucscTableQuery(session, "knownGene",
        GRangesForUCSCGenome(genome, chromosome,
            IRanges(minRange, maxRange)), table = tbl,
        names = t[, "name"])
    t1 <- getTable(query2)
    t1[, "kgID"] <- as.character(t1[, "kgID"])
    t1[, "geneSymbol"] <- as.character(t1[, "geneSymbol"])
    for (i in 1:dim(t)[1]) {
        gene_name <- t1[t1[, "kgID"] == t[i, "name"], "geneSymbol"]
        if (length(gene_name) != 0)
            t[i, "gene_name"] <- gene_name
    }
    if ("kgColor" %in% tableNames(query1)) {
        query3 <- ucscTableQuery(session, "knownGene",
            GRangesForUCSCGenome(genome, chromosome,
                IRanges(minRange, maxRange)), table = "kgColor",
            names = t[, "name"])
        color_tbl <- getTable(query3)
    }
    t$plot_line <- 0
    t$plot_line[1] <- 1
    plot_lines_no <- 1
    if (dim(t)[1] > 1) {
        for (i in 2:dim(t)[1]) {
            gene_name_width <- Range/map_len * convertWidth(grobWidth(textGrob(paste(t[1,
                "gene_name"], "   "), gp = gpar(fontsize = 7))),
                "npc", valueOnly = TRUE)
            for (j in 1:plot_lines_no) {
                if (max(t[1:(i - 1), ][t[, "plot_line"] == j,
                  "txEnd"], na.rm = TRUE) < t[i, "txStart"] -
                  gene_name_width) {
                  t[i, "plot_line"] <- j
                  break
                }
            }
            if (!t[i, "plot_line"])
                t[i, "plot_line"] <- plot_lines_no <- plot_lines_no +
                  1
        }
    }
    height = plot_lines_no * plot_lines_distance
    vp$height <- unit(height, "npc")
    # txt = "UCSC Genes Based on RefSeq, UniProt, GenBank, CCDS and Comparative Genomics"
    txt=''
    gene_plot_title <- textGrob(txt,
        gp = gpar(fontsize = 7, fontfamily = "mono"), y = 1,
        just = c("centre", "center"), name = "gene_plot_title",
        default.units = "native")
    Transcripts <- gTree(children = gList(gene_plot_title), name = "transcripts",
        vp = vp)
    for (i in 1:dim(t)[1]) {
        tx_upArrows <- tx_downArrows <- tx_exons <- tx_cds <- tx_leftArrow <- tx_rightArrow <- NULL
        y = 1 - t[i, "plot_line"]/plot_lines_no
        txStart <- max(minRange, t[i, "txStart"])
        txEnd <- min(maxRange, t[i, "txEnd"])
        cdsStart <- max(minRange, t[i, "cdsStart"])
        cdsEnd <- min(maxRange, t[i, "cdsEnd"])
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
            if (cdsStart < exonEnds[j] & cdsEnd > exonStarts[j]) {
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
        if (exists("color_tbl")) {
            r <- color_tbl[color_tbl[, "kgID"] == t[i, "name"],
                "r"]
            g <- color_tbl[color_tbl[, "kgID"] == t[i, "name"],
                "g"]
            b <- color_tbl[color_tbl[, "kgID"] == t[i, "name"],
                "b"]
            gene_color <- rgb(r, g, b, maxColorValue = 255)
        }
        transcript <- gTree(children = gList(tx_name, tx_region,
            tx_upArrows, tx_downArrows, tx_exons, tx_cds, tx_leftArrow,
            tx_rightArrow), gp = gpar(col = gene_color), name = t[i,
            "name"])
        Transcripts <- addGrob(Transcripts, transcript)
    }
    # grid.draw(Transcripts)
    return(Transcripts)
}

# adapted from LDheatmap:::LDheatmapMapNew.add
addMap = function (nsnps, add.map, genetic.distances, geneMapLocation = 0.15,
    geneMapLabelX = NULL, geneMapLabelY = NULL, distances = "physical",
    vp = NULL, SNP.name = NULL, ind = 0, flip = FALSE)
{
    snp <- ((1:nsnps - 1) + 0.5)/nsnps
    if (add.map) {
        min.dist <- min(genetic.distances)
        max.dist <- max(genetic.distances)
        total.dist <- max.dist - min.dist
        if (flip)
            geneMapLocation <- (-geneMapLocation)
        seq.x <- c(0.5 * geneMapLocation + 1/(nsnps * 2), 1 +
            0.5 * geneMapLocation - 1/(nsnps * 2))
        seq.y <- c(-0.5 * geneMapLocation + 1/(nsnps * 2), 1 -
            0.5 * geneMapLocation - 1/(nsnps * 2))
        diagonal <- linesGrob(seq.x, seq.y, gp = gpar(lty = 1),
            name = "diagonal", vp = vp)
        regionx <- seq.x[1] + ((genetic.distances - min.dist)/total.dist) *
            (seq.x[2] - seq.x[1])
        regiony <- seq.y[1] + ((genetic.distances - min.dist)/total.dist) *
            (seq.y[2] - seq.y[1])
        segments <- segmentsGrob(snp, snp, regionx, regiony,
            name = "segments", vp = vp)
        # if (distances == "physical")
        #     mapLabel <- paste("Physical Length:", round((total.dist/1000),
        #         1), "kb", sep = "")
        # else mapLabel <- paste("Genetic Map Length:", round(total.dist,
        #     1), "cM", sep = "")
        if (!flip) {
            if (is.null(geneMapLabelY))
                geneMapLabelY <- 0.3
            if (is.null(geneMapLabelX))
                geneMapLabelX <- 0.5
        }
        else {
            if (is.null(geneMapLabelY))
                geneMapLabelY <- 0.8
            if (is.null(geneMapLabelX))
                geneMapLabelX <- 0.4
        }
        # title <- textGrob(mapLabel, geneMapLabelX, geneMapLabelY,
        #     gp = gpar(cex = 0.9), just = "left", name = "title")
        geneMap <- gTree(children = gList( segments), name = "geneMap")
        if (!is.null(SNP.name) && (any(ind != 0))) {
            if (flip) {
                length_SNP_name <- max(nchar(SNP.name))
                long_SNP_name <- paste(rep(8, length_SNP_name),
                  collapse = "")
                name_gap <- convertWidth(grobWidth(textGrob(long_SNP_name)),
                  "npc", valueOnly = TRUE)/sqrt(2)
                diagonal <- linesGrob(seq.x, seq.y, gp = gpar(lty = 1),
                  name = "diagonal", vp = vp)
                segments <- segmentsGrob(snp, snp, regionx, regiony,
                  name = "segments", vp = vp)
                symbols <- pointsGrob(snp[ind], snp[ind], pch = "*",
                  gp = gpar(cex = 1.25, bg = "blue", col = "blue"),
                  name = "symbols", vp = vp)
                SNPnames <- textGrob(SNP.name, just = "left",
                  rot = -45, regionx[ind] - sqrt(2 + 0.5) * name_gap,
                  regiony[ind] + sqrt(2 + 0.5) * name_gap, gp = gpar(cex = 0.6,
                    col = "blue"), name = "SNPnames", vp = vp)
                # title <- editGrob(title, y = unit(geneMapLabelY +
                  # name_gap, "npc"))
            }
            else {
                symbols <- pointsGrob(snp[ind], snp[ind], pch = "*",
                  gp = gpar(cex = 1.25, bg = "blue", col = "blue"),
                  name = "symbols", vp = vp)
                SNPnames <- textGrob(paste(" ", SNP.name), just = "left",
                  rot = -45, regionx[ind], regiony[ind], gp = gpar(cex = 0.6,
                    col = "blue"), name = "SNPnames", vp = vp)
            }
            geneMap <- gTree(children = gList( segments,
                symbols, SNPnames), name = "geneMap")
        }
    }
    else if (!add.map && !is.null(SNP.name) && (any(ind != 0))) {
        geneMap <- textGrob(paste(" ", SNP.name), just = "left",
            rot = -45, snp[ind], snp[ind], gp = gpar(cex = 0.6,
                col = "blue"), name = "SNPnames")
        if (flip)
            geneMap <- editGrob(geneMap, vp = vp)
    }
    else geneMap <- NULL
    geneMap
}