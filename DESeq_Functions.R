## prepare deg table
prepare_deg_table <- function(res,
                              padj_cutoff = 0.05,
                              lfc_cutoff = 0) {
  
  deg_table <- res %>%
    data.frame() %>%
    mutate(geneid = rownames(.)) %>%
    filter(!is.na(padj)) %>%
    mutate(
      diffexpressed = case_when(
        padj > padj_cutoff ~ "NS",
        log2FoldChange > lfc_cutoff & padj < padj_cutoff ~ "UP",
        log2FoldChange < -lfc_cutoff & padj < padj_cutoff ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  
  return(deg_table) 
  
}

## volcano plot
library(ggrepel)
volcano_plot <- function(deg_table,
                         ntop_deg_labels = 20,
                         lfc_cutoff = 0,
                         padj_cutoff = 0.05) {
  
  top20 <- deg_table%>%
    filter(diffexpressed != "NS")%>%
    arrange(padj)%>%
    slice_head(n = ntop_deg_labels)
  
  
  p <- ggplot(deg_table, aes(x=log2FoldChange, y=-log10(padj), color = diffexpressed)) +
    geom_point(size = 1) +
    geom_vline(xintercept = c(lfc_cutoff, -lfc_cutoff), linetype= 'dashed', color = 'gray') +
    geom_hline(yintercept = -log10(padj_cutoff),linetype= 'dashed', color = 'gray') +
    scale_color_manual(values = c(DOWN = '#2E6DB4', NS = "gray40", UP="#D6403A"))+ 
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_cartesian() +
    ggtitle('Dex: trt vs. untrt')+
    geom_text_repel(data = top20, aes(label = geneid), size = 3, max.overlaps = Inf) 
  
  return(p)
}
