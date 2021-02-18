### custom NicheNet functions ###
## load packages
suppressPackageStartupMessages(library(nichenetr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(plotrix))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(viridis))

## function to perform ligand activity analysis ##
ligand_activity_analysis <- function(sender_ct, receiver_ct, geneset_oi, 
                            geneset_title,
                            pct_expr_table, pct_thresh = 0.10,
                            pearson_thresh = 0.08,
                            save_top_ligands = TRUE){
    
    
    ## retrieve genes expressed by receiver
    expr_genes_receiver = rownames(pct_expr_table[pct_expr_table[, receiver_ct] > pct_thresh, ])
    background_expr_genes = expr_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    ## retrieve genes expressed by sender
    list_expr_genes_sender = lapply(sender_ct, function(x){rownames(pct_expr_table[pct_expr_table[, x] > pct_thresh, ])})
    expr_genes_sender = list_expr_genes_sender %>% unlist() %>% unique()
    
    ## status message ##
    print(paste0("Using ", length(geneset_oi), " genes differently regulated genes in ", receiver_ct, " (",
                 geneset_title, ")"))
    
    ### STEP1: Ligand activity analysis ###
    ## Define a set of potential ligands and receptors 
    # retrieve ligands and receptors
    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()
    # ligands expressed by sender celltypes
    expr_ligands = intersect(ligands, expr_genes_sender) 
    # receptor expressed by receiver celltypes
    expr_receptors = intersect(receptors, expr_genes_receiver)
    ### status messages ###
    print(paste0("Expressed Ligands ", length(expr_ligands)))
    print(paste0("Expressed Receptors ", length(expr_receptors)))
    
    ## filter ligands
    # only consider ligands with matching receptors (according to NicheNets databases)
    potential_ligands = lr_network %>% filter(from %in% expr_ligands & to %in% expr_receptors) %>%
                        pull(from) %>% unique()
    ### status message ###
    print(paste0("Potential Ligands ", length(potential_ligands)))
    
    ## predict ligand activities
    ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                                  background_expressed_genes = background_expr_genes,
                                                  ligand_target_matrix = ligand_target_matrix,
                                                  potential_ligands = potential_ligands)

    ## rank ligands by pearson correlation coefficient
    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
    
    # filter consider ligands with pearson's correlation >= pearson_tresh
    ligand_activities = ligand_activities %>% filter(pearson >= pearson_thresh)
    
    ### status message ###
    print(paste0("Top ranked ligands ", length(ligand_activities$test_ligand)))
    
    return(ligand_activities)  
}

## function to perform ligand-target analysis ##
ligand_target_analysis <- function(sender_ct, receiver_ct, geneset_oi, 
                            geneset_title,
                            pct_expr_table, pct_thresh = 0.10,
                            avg_expr_table,
                            pearson_thresh = 0.08, n_top_ligands = 25,
                            save_top_ligands = TRUE,
                            n_targets = 100, target_thresh = 0.33,
                            save_ligand_targets = TRUE,
                            save_fig = TRUE){
    
    
    ## retrieve genes expressed by receiver
    expr_genes_receiver = rownames(pct_expr_table[pct_expr_table[, receiver_ct] > pct_thresh, ])
    background_expr_genes = expr_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    ## retrieve genes expressed by sender
    list_expr_genes_sender = lapply(sender_ct, function(x){rownames(pct_expr_table[pct_expr_table[, x] > pct_thresh, ])})
    expr_genes_sender = list_expr_genes_sender %>% unlist() %>% unique()
    
    ## status message ##
    print(paste0("Using ", length(geneset_oi), " genes differently regulated genes in ", receiver_ct, " (",
                 geneset_title, ")"))
    
    ### STEP 1: Ligand activity analysis ###
    ## Define a set of potential ligands and receptors 
    # retrieve ligands and receptors
    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()
    # ligands expressed by sender celltypes
    expr_ligands = intersect(ligands, expr_genes_sender) 
    # receptor expressed by receiver celltypes
    expr_receptors = intersect(receptors, expr_genes_receiver)
    ## status messages ##
    print(paste0("Expressed Ligands ", length(expr_ligands)))
    print(paste0("Expressed Receptors ", length(expr_receptors)))
    
    ## filter ligands
    # only consider ligands with matching receptors (according to NicheNets databases)
    potential_ligands = lr_network %>% filter(from %in% expr_ligands & to %in% expr_receptors) %>%
                        pull(from) %>% unique()
    ## status message ##
    print(paste0("Potential Ligands ", length(potential_ligands)))
    
    ## predict ligand activities
    ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                                  background_expressed_genes = background_expr_genes,
                                                  ligand_target_matrix = ligand_target_matrix,
                                                  potential_ligands = potential_ligands)

    ## rank ligands by pearson correlation coefficient
    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
    
    # filter consider ligands with pearson's correlation >= pearson_tresh
    ligand_activities = ligand_activities %>% filter(pearson >= pearson_thresh)
    ## status message ##
    print(paste0("Top ranked ligands ", length(ligand_activities$test_ligand)))
    
    ## save top predicted ligands table to working directory ##
    if(save_top_ligands == TRUE){
        write.csv(ligand_activities, file = "top_ranked_ligands.csv")
    }
    
    ## filter top ranked ligands
    best_upstream_ligands = ligand_activities %>% top_n(n_top_ligands, pearson) %>% arrange(-pearson) %>% 
    pull(test_ligand) %>% unique()
    ## status message ##
    if(length(best_upstream_ligands) == 0){
        print("Predicted 0 top ligands. Unable to proceed with Target Gene Prediction")
        return 
    }    
    print(paste0("Continue with top ", length(best_upstream_ligands), " highest ranked ligands (pearson's rho)"))
    
    ### STEP 2: Target Gene Predicton
    
    ## identify ligand targets
    active_ligand_target_links_df = best_upstream_ligands %>% 
                                    lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = n_targets) %>% bind_rows() %>% drop_na()
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = target_thresh)
    ## status message ##
    print(paste0("Predicting target genes"))
    
    ## reformat data
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    order_targets = active_ligand_target_links_df$target %>% unique() %>% 
                    intersect(rownames(active_ligand_target_links)) %>% make.names()
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() 
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
    
    ## final ligand-target heatmap
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
    if(dim(vis_ligand_target)[2] == 0){
        print("Predicted 0 Target Genes. Unable to compute ligand-target heatmap")
        return 
    }
    ## save ligand-target heatmap to working directory ##
    if(save_ligand_targets == TRUE){
        write.csv(vis_ligand_target, file = "ligand_target_heatmap.csv")
    }
    ## status message ##
    print(paste0("Computing final ligand-target heatmap"))
    
    ### STEP 3: Data preparation ###
    
    ## left panel: ligand pearson correlation heatmap
    ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% 
    magrittr::set_rownames(ligand_activities$test_ligand)
    rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
    colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

    vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% 
                            magrittr::set_colnames("Pearson")

    ## status message ##
    print(paste0("Printing left panel: ligand pearson correlation heatmap"))
    LigandPearsonCor_df = vis_ligand_pearson %>% data.frame() %>% rownames_to_column("y") %>%
                            tbl_df() %>% gather(x, "score", -y) %>% 
                            mutate(y = factor(y, levels = rownames(vis_ligand_pearson), ordered = TRUE), 
                                               x = factor(x, levels = colnames(vis_ligand_pearson), ordered = TRUE)) %>% 
                            mutate(score = ifelse(score == 0, NA, score))
    
    ## middle panel: scaled of expression of ligands in SENDER celltypes
    # filter for SENDER celltypes only
    top_ligands_avg_expression = avg_expr_table %>% select(all_of(sender_ct))
    
    # prepare data
    top_ligands_avg_expression <- data.frame(t(top_ligands_avg_expression))
    colnames(top_ligands_avg_expression) <- rownames(vis_ligand_pearson)
    top_ligands_avg_expression = top_ligands_avg_expression %>% select(all_of(rownames(vis_ligand_pearson)))
    top_ligands_avg_expression$celltype <- rownames(top_ligands_avg_expression)
    
    # center and scale expression values per gene
    center_and_scale <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
    top_ligands_avg_expression <- top_ligands_avg_expression %>% mutate_at(rownames(vis_ligand_pearson), center_and_scale)
   
    # reformat data
    rownames(top_ligands_avg_expression) <- top_ligands_avg_expression$celltype
    top_ligands_avg_expression$celltype <- NULL
    top_ligands_avg_expression <- t(top_ligands_avg_expression)
    top_ligands_avg_expression_rows <- rownames(top_ligands_avg_expression)
    top_ligands_avg_expression_cols <- colnames(top_ligands_avg_expression)
    top_ligands_avg_expression = top_ligands_avg_expression %>% data.frame()
    rownames(top_ligands_avg_expression) <- top_ligands_avg_expression_rows
    colnames(top_ligands_avg_expression) <- top_ligands_avg_expression_cols

    ## status message ##
    print(paste0("Printing center panel: scaled ligand expression in Sender celltypes"))
    LigandExpr_df = top_ligands_avg_expression %>% rownames_to_column("y") %>% 
                    tbl_df() %>% gather(x, "score", -y) %>% 
                    mutate(y = factor(y, levels = rownames(top_ligands_avg_expression), ordered = TRUE), 
                           x = factor(x, levels = colnames(top_ligands_avg_expression), ordered = TRUE)) %>% 
                    mutate(score = ifelse(score == 0, NA, score))
        
    ## right panel: ligand target heatmap
    ## status message ##
    print(paste0("Printing right panel: ligand target heatmap"))
    LigandTarget_df = vis_ligand_target %>% data.frame() %>% rownames_to_column("y") %>% tbl_df() %>% 
                        gather(x, "score", -y) %>% 
                        mutate(y = factor(y, levels = rownames(vis_ligand_target), ordered = TRUE), 
                                x = factor(x, levels = colnames(vis_ligand_target), ordered = TRUE)) %>% 
                        mutate(score = ifelse(score == 0, NA, score))
    
    ### STEP 4: Data visualization ###
    
    ## left panel: ligand pearson correlation heatmap
    plot_LigandPearsonCor <- LigandPearsonCor_df %>% 
    ggplot(aes(x, y, fill = score)) + 
    geom_tile(color = "white", size = 0.5) + 
    coord_equal() +
    scale_fill_gradient(low = "gold", high = "firebrick", 
                        na.value = "whitesmoke") + 
    theme_minimal() + 
    theme(panel.grid.minor = element_line(color = "transparent"), 
          panel.grid.major = element_line(color = "transparent"),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom", 
          legend.title = element_text(size = 10, angle = 0),
          legend.text = element_text(size = 8, hjust = 1), 
          axis.ticks = element_line(size = 0), 
          axis.text.x.top = element_text(angle = 90, hjust = 0), 
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1), 
          axis.title.y = element_text(size = 14), 
          axis.text.y = element_text(size = 14)) + 
          scale_x_discrete(position = "top", labels = "Pearson") + 
          xlab(paste0("")) + 
          ylab(paste0("Top ligands")) + 
    labs(title = paste0("Ligand activity")) +
    guides(fill = guide_colourbar("Pearson correlation coefficient\ntarget gene prediction ability", 
                                  title.position = "top", angle = 90))
    #plot_LigandPearsonCor
    
    ## middle panel: scaled of expression of ligands in SENDER celltypes
    plot_LigandExpr <- LigandExpr_df %>% 
    ggplot(aes(x, y, fill = score)) + 
    geom_tile(color = "white", size = 0.5) + 
    coord_equal() +
    scale_fill_gradient2(low = 'deepskyblue3', mid = 'white', high = 'red3') +
    theme_minimal() +
    theme(panel.grid.minor = element_line(color = "transparent"), 
          panel.grid.major = element_line(color = "transparent"),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = 10, angle = 0),
          legend.text = element_text(size = 8, hjust = 1),  
          axis.ticks = element_line(size = 0), 
          axis.text.x.top = element_text(angle = 90, hjust = 0), 
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1), 
          axis.title = element_text(size = 14), 
          axis.text.y = element_text(size = 14)) + 
          scale_x_discrete(position = "top") + 
          xlab(paste0("")) + 
          ylab(paste0("")) + 
    labs(title = paste0("Expression in Sender")) +
    guides(fill = guide_colourbar("Relative expression of ligands\n(centered and scaled)", title.position = "top"))
    #plot_LigandExpr
    
    ## right panel: ligand target heatmap
    plot_LigandTarget <- LigandTarget_df %>% 
    ggplot(aes(x, y, fill = score)) + 
    geom_tile(color = "white", size = 0.5) + 
    coord_equal() +
    scale_fill_gradient(low = "#E8EAF6", high = "magenta", na.value = "whitesmoke")  + 
    #scale_fill_material("green", na.value = "whitesmoke") +
    theme_minimal() + 
    theme(panel.grid.minor = element_line(color = "transparent"), 
          panel.grid.major = element_line(color = "transparent"),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = 10, angle = 0),
          legend.text = element_text(size = 8, hjust = 1),  
          axis.ticks = element_line(size = 0), 
          axis.text.x.top = element_text(angle = 90, hjust = 0), 
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1),
          axis.title = element_text(size = 14), 
          axis.text.y = element_text(size = 14)) + 
          scale_x_discrete(position = "top") + 
          xlab(paste0("")) + 
          ylab(paste0("")) + 
    labs(title = paste0(receiver_ct, " geneset: ", geneset_title)) +
    guides(fill = guide_colourbar("Regulatory potential", title.position = "top"))
    #plot_LigandTarget
    
    ### STEP 5: Build Summary fiugure
    ## status message ##
    print(paste0("Composing summary figure"))
    summary_figure <- plot_LigandPearsonCor + plot_LigandExpr + plot_LigandTarget + plot_layout(guides = 'collect')
    ## save summary figure 
    if(save_fig == TRUE){
        png(paste0("summary_figure_1.png"), width=2000, height=1000, units="px")
        print(summary_figure)
        dev.off()
    }
    ## status message ##
    print(paste0("Save summary figure to ", getwd(), "/summary_figure_1.png"))
    return(summary_figure)  
}

## function to perform ligand-receptor analysis ##
receptor_ligand_analysis <- function(ligand_activities, receiver_ct, background_ct,
                                     n_top_ligands = 25,
                                     pct_expr_table, pct_thresh = 0.1,
                                     avg_expr_table,
                                     network_thresh = 0.4,
                                     save_ligand_receptor_network = TRUE,
                                     save_fig = TRUE){
    
    ### STEP 1: receptor prediction ###
    ## filter top N top ranked ligands
    top_ligands = ligand_activities %>% top_n(n_top_ligands, pearson) %>% arrange(-pearson) %>% 
    pull(test_ligand) %>% unique()
    
    ## get expressed receptors ##
    expr_genes_receiver = rownames(pct_expr_table[pct_expr_table[, receiver_ct] > pct_thresh, ])
    background_expr_genes = expr_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    ## retrieve receptors from NicheNet's receptor ligand network
    receptors = lr_network %>% pull(to) %>% unique()
    ## receptor expressed by receiver celltypes
    expr_receptors = intersect(receptors, expr_genes_receiver)
    ## status message ##
    print(paste0("Expressed Receptors ", length(expr_receptors)))
    
    ## predict receptors ##
    lr_network_top = lr_network %>% filter(from %in% top_ligands & to %in% expr_receptors) %>% distinct(from,to)
    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
    
    lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% top_ligands)
    lr_network_top_df_large = lr_network_top_df_large %>% filter(weight > network_thresh )
    
    lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
    lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
    
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
    
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
    
    order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
    order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
    
    ## status message ##
    print(paste0("Computing ligand-receptor network"))
    vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
    rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
    vis_ligand_receptor_network = vis_ligand_receptor_network %>% t()
    if(save_ligand_receptor_network == TRUE){
        write.csv(vis_ligand_receptor_network, file = "ligand_receptor_network.csv")
    }
    
    ### STEP 2: data preparation ###
    
    ## top panel: ligand-receptor heatmap ##
    ## status message ##
    print(paste0("Printing top panel: ligand receptor heatmap"))
    LigandReceptor_df = vis_ligand_receptor_network %>% data.frame() %>% rownames_to_column("y") %>% tbl_df() %>% 
                        gather(x, "score", -y) %>% 
                        mutate(y = factor(y, levels = rownames(vis_ligand_receptor_network), ordered = TRUE), 
                                x = factor(x, levels = colnames(vis_ligand_receptor_network), ordered = TRUE)) %>% 
                        mutate(score = ifelse(score == 0, NA, score))
    
    ## bottom panel: scaled expression of top receptors ##
    top_receptors <- colnames(vis_ligand_receptor_network)
    
    # filter for SENDER celltypes only
    top_receptors_avg_expression = avg_expr_table %>% select(all_of(background_ct))
    
    # prepare data
    top_receptors_avg_expression <- data.frame(t(top_receptors_avg_expression))
    colnames(top_receptors_avg_expression) <- top_receptors
    top_receptors_avg_expression = top_receptors_avg_expression %>% select(all_of(top_receptors))
    top_receptors_avg_expression$celltype <- rownames(top_receptors_avg_expression)
    
    ## center and scale expression values per gene
    center_and_scale <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
    top_receptors_avg_expression <- top_receptors_avg_expression %>% mutate_at(top_receptors, center_and_scale)
    
    # reformat data
    rownames(top_receptors_avg_expression) <- top_receptors_avg_expression$celltype
    top_receptors_avg_expression$celltype <- NULL
    
    ## reformat data 
    ## status message ##
    print(paste0("Printing bottom panel: scaled receptor expression across ALL celltypes"))
    ReceptorExpr_df = top_receptors_avg_expression %>% rownames_to_column("y") %>% 
                    tbl_df() %>% gather(x, "score", -y) %>% 
                    mutate(y = factor(y, levels = rownames(top_receptors_avg_expression), ordered = TRUE), 
                           x = factor(x, levels = colnames(top_receptors_avg_expression), ordered = TRUE)) %>% 
                    mutate(score = ifelse(score == 0, NA, score))
    
    ### STEP 3: data visualization ###
    
    # top panel: ligand-receptor heatmap ##
    plot_LigandReceptor <- LigandReceptor_df %>% 
    ggplot(aes(x, y, fill = score)) + 
    geom_tile(color = "white", size = 0.5) + 
    coord_equal() +
    #scale_fill_gradient(low = "#E8EAF6", high = "magenta", na.value = "whitesmoke")  + 
    scale_fill_material("pink", na.value = "whitesmoke") +
    theme_minimal() + 
    theme(panel.grid.minor = element_line(color = "transparent"), 
          panel.grid.major = element_line(color = "transparent"),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = 10, angle = 0),
          legend.text = element_text(size = 8, hjust = 1),  
          axis.ticks = element_line(size = 0), 
          axis.text.x.top = element_text(angle = 90, hjust = 0), 
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1),
          axis.title = element_text(size = 14), 
          axis.text.y = element_text(size = 14)) + 
          scale_x_discrete(position = "top") + 
          xlab(paste0("")) + 
          ylab(paste0("")) + 
    labs(title = paste0(receiver_ct, " receptors")) +
    guides(fill = guide_colourbar("Regulatory potential", title.position = "top"))
    
    # bottom panel: scaled expression of top receptor across ALL celltypes
    plot_ReceptorExpr <- ReceptorExpr_df %>% 
    ggplot(aes(x, y, fill = score)) + 
    geom_tile(color = "white", size = 0.5) + 
    coord_equal() +
    #scale_fill_gradient2(low = 'deepskyblue3', mid = 'white', high = 'red3') +
    scale_fill_gsea() +
    theme_minimal() +
    theme(panel.grid.minor = element_line(color = "transparent"), 
          panel.grid.major = element_line(color = "transparent"),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = 10, angle = 0),
          legend.text = element_text(size = 8, hjust = 1),  
          axis.ticks = element_line(size = 0), 
          axis.text.x.top = element_text(angle = 90, hjust = 0), 
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1), 
          axis.title = element_text(size = 14), 
          axis.text.y = element_text(size = 14)) + 
          scale_x_discrete(position = "top") + 
          xlab(paste0("")) + 
          ylab(paste0("")) + 
    labs(title = paste0("Expression in tissue niche")) +
    guides(fill = guide_colourbar("Relative expression of receptors\n(centered and scaled)", title.position = "top"))
    
    ## STEP 4: compose summary figure ###
    summary_figure <- plot_LigandReceptor / plot_ReceptorExpr + plot_layout(guides = 'collect')
    if(save_fig == TRUE){
        png(paste0("summary_figure_2.png"), width=2000, height=2000, units="px")
        print(summary_figure)
        dev.off()
    }
    # status message ##
    print(paste0("Save summary figure to ", getwd(), "/summary_figure_2.png"))
    return(summary_figure)  
    
}

## function to assess cell-cell interactions QUANTITATIVELY
quant_interactions <- function(ligand_activities, sender_ct, receiver_ct, 
                               n_top_ligands = 25, pct_expr_table, pct_thresh = 0.1){
     
    # filter top N top ranked ligands
    ligand_pearson = ligand_activities %>% top_n(n_top_ligands, pearson) %>% arrange(-pearson) %>% select(test_ligand) %>% arrange(test_ligand)
    
    ## for each predicted ligand, check in which celltypes they are expresed above pct.tresh
    # build ligand count vector
    ligand_count <- data.frame(row.names = receiver_ct)
    # subset pct epxressed table
    pct_expr = pct_expr_table %>% filter(rownames(pct_expr_table) %in% ligand_pearson$test_ligand) %>% select(all_of(sender_ct)) 
    
    for (ct in colnames(pct_expr)){
        ligand_count[, ct] = sum(pct_expr[ligand_pearson$test_ligand, ct] > pct_thresh)
    }
    
    # scale function 
    scale <- function(x)((x - min(x))/(max(x)-min(x)))
    # create results data.frame
    q_interactions <- data.frame(cell_type = colnames(ligand_count),
                      receiver = rep(receiver_ct, length(colnames(ligand_count))),
                      nr_intact = as.numeric(ligand_count[1,]))
    # scale adjacency matrix
    q_interactions$nr_intact <- scale(q_interactions$nr_intact)*100+1
    # round results
    q_interactions$nr_intact <- round(q_interactions$nr_intact, 1)
    # order output data.frame
    q_interactions = q_interactions %>% arrange(cell_type)
    
    ## return results
    return(q_interactions)
}

## function to assess cell-cell interaction QUALITATIVELY
qual_interactions <- function(ligand_activities, n_top_ligands = 25,
                              sender_ct, receiver_ct, avg_expr_table){
        
    ## step 1: create data.frame ligand x celltype
    
    # filter top N top ranked ligands
    ligand_pearson = ligand_activities %>% top_n(n_top_ligands, pearson) %>% arrange(-pearson) %>%  select(test_ligand, pearson)
    ligand_pearson = ligand_pearson %>% arrange(test_ligand)
    
    ## step 2: compute contribution of each sender celltype to ligand effect
    # average of all top ranked ligands in all cell types at timepoint after treatment
    sender_avg_expr = avg_expr_table %>% select(all_of(sender_ct)) %>% 
                      filter(rownames(avg_expr_table) %in% ligand_pearson$test_ligand)
    sender_avg_expr = sender_avg_expr %>% arrange(rownames(sender_avg_expr))

    # scale average expression per ligand between 0 and 1
    # scaling function
    scale <- function(x) ((x - min(x)) / (max(x) - min(x)))
    # transpose data (because normalizing data is more elegant to do columnwise)
    rownames(sender_avg_expr) = rownames(sender_avg_expr) %>% make.names() # make.names() for heatmap 
    sender_avg_expr = sender_avg_expr %>% t() %>% data.frame()
    # normalize columnwise
    sender_avg_expr = sender_avg_expr %>% mutate_at(make.names(ligand_pearson$test_ligand), scale)
    # restore rownames
    rownames(sender_avg_expr) <- sender_ct
    ## transpose again to get ligand x celltype data.frame
    sender_avg_expr = sender_avg_expr %>% t() %>% data.frame()
    colnames(sender_avg_expr) <- sender_ct
    ### status message ###
    print("Computing contribution of each sender celltype to effect of each ligand")

    ## step 3: sum all ligand activity scores for each cell type
    # assertion check: check whether ligands in 'ligand_pearson' and 'avg_expr_table' have the same order
    stopifnot(make.names(ligand_pearson$test_ligand) == rownames(sender_avg_expr))
    # multiply pearson of each ligand with average expression of the corresponding cell type 
    # to quantify contribution of each sender to effect of a given ligand
    weight_celltypes <- function(x)(x*ligand_pearson$pearson)
    sender_avg_expr = sender_avg_expr %>% mutate_at(colnames(sender_avg_expr), weight_celltypes)
    # sum all contributions per cell to get one score per cell
    sender_avg_expr = sender_avg_expr %>% summarise_at(colnames(sender_avg_expr), sum)
    ### status message ###
    print("Summing contributions of each sender celltype")
    
    ## clean data
    # create results data.frame
    q_interactions <- data.frame(cell_type = colnames(sender_avg_expr),
                      receiver = rep(receiver_ct, length(colnames(sender_avg_expr))),
                      nr_intact = as.numeric(sender_avg_expr[1,]))
    # scale results between 0 and 10 
    q_interactions$nr_intact <- scale(q_interactions$nr_intact)*100 + 1# to make differences more apparent
    # round results
    q_interactions$nr_intact <- round(q_interactions$nr_intact, 1)
    # order output data.frame
    q_interactions = q_interactions %>% arrange(cell_type)
    
    ## return results
    return(q_interactions)
}

## function to create CONNECTOME graph
create_connectome_graph <- function(interactions, color_pal, min_intact = 0,
                            node_size = 20, width_factor = 1, text_size = 1.0){
    
    ## check if color_pal is one of the 19 material design color palettes
    color_pal_avail <- c("red", "pink", "purple", "deep-purple", "indigo", "blue", "light-blue", "cyan",
                         "teal", "green", "light-green", "lime", "yellow", "amber", "orange", "deep-orange",
                         "brown", "grey", "blue-grey")
    if(!(color_pal %in% color_pal_avail)){
        print(paste0("ERROR: Selected color palette \'", color_pal, "\' is not available, please choose one of the following:\n", color_pal_avail))
        return
    }
    # create graph from interaction results
    G <- graph_from_edgelist(as.matrix(interactions[, 1:2]), directed = F)
    E(G)$width <- interactions$nr_intact * width_factor # edge width
    E(G)$color <- color.scale(interactions$nr_intact, extremes = c(pal_material(color_pal)(10)[1], pal_material(color_pal)(10)[10])) # edge color gradient
    V(G)$label.cex <- text_size # text size
    
    return(G)
}

## function to plot CONNECTOME graph
plot_connectome <- function(G, coordinates, color_pal, title, title_size = 2, n_connectome_graph, save_fig = TRUE){
    
    ## check if color_pal is one of the 19 material design color palettes
    color_pal_avail <- c("red", "pink", "purple", "deep-purple", "indigo", "blue", "light-blue", "cyan",
                         "teal", "green", "light-green", "lime", "yellow", "amber", "orange", "deep-orange",
                         "brown", "grey", "blue-grey")
    if(!(color_pal %in% color_pal_avail)){
        print(paste0("ERROR: Selected color palette \'", color_pal, 
                     "\' is not available, please choose one of the following:\n", color_pal_avail))
        return
    }
    ### status message ###
    print("Plotting connectome plot")
    ## create connectome plot
    plot(G, edge.width = E(G)$width, vertex.size = 10, vertex.label.family = "Helvetica",
                           vertex.label.color= "black", vertex.color = pal_material(color_pal)(10)[1],
                           layout = coordinates) + title(paste(title), cex.main = title_size)
    connectome_plot <- recordPlot()
    if(save_fig == TRUE){
       png(paste0("connnectome_plot_", n_connectome_graph , ".png"), width=1000, height=1000, units="px")
       print(connectome_plot)
       dev.off() 
    }
    return(connectome_plot)
}