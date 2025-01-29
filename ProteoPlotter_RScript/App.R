#Esther Olabisi-Adeniyi
#Jennifer Geddes-McAlister Lab - University of Guelph
#August 8, 2023

#All packages used
# my.packages <- c("ggplot2", "ggpubr", "ggrepel", "tidyverse", "DT", "reshape2", "scales", "RColorBrewer", "colourpicker", "ggnewscale", "ggfortify", "shiny", "shinyjs", "shinyWidgets", "eulerr", "webshot", "shinydashboard")
#Check if packages are installed
# installed.pkgs <- my.packages %in% rownames(installed.packages())
#Install packages that aren't already installed
# if(any(installed.pkgs == FALSE)) {
# install.packages(my.packages[!installed.pkgs])
# }
#Load all packages
#invisible(lapply(my.packages, library, character.only = T))

#Option to run packages individually when running app on Shinyapps.io
library(ggplot2); library(ggpubr); library(ggrepel); library(tidyverse); library(DT); library(reshape2); library(scales); library(RColorBrewer); library(colourpicker); library(ggnewscale); library(ggfortify); library(shiny); library(shinyjs); library(shinyWidgets); library(ggvenn); library(shinyBS); library(eulerr); library(shinydashboard); library(UpSetR); library(purrr); library(spsComps); library(dplyr); library(randomcoloR); library(tools)


#### EXTRACT IDENTIFIERS ----
extract_identifiers <- function(df){
  identifiers <- df[,1] %>% 
    grep(pattern = "^#", ., value = T) %>%
    sub(pattern = ".*\\{(.*)\\}.*", replacement = "\\1", ., perl = TRUE) %>% 
    sub(pattern = ".*:", replacement = "", ., perl = TRUE)
  identifiers <- identifiers[-(grep("Typ", identifiers))] #remove "Type" option
  return(identifiers)
} 


####1. PROCESS DATA MATRIX ----
#Read protein.txt matrix processed in Perseus for a volcano plot. Remove all comment rows. 
#Change column names to shorter generic names that are easier to work with. 
#prot.dataset: Imported from Perseus following quality control, two samples t-test, and GO enrichment
#gid = group identifier for recognizing the "intensity" columns

process_df <- function(prot.dataset, fdr = 0.05, s0 = 1, gid = NULL) {
  col.names <- names(prot.dataset)
  gid <- ifelse(is.null(gid), "group", gid)
  gid.index <- grep(gid, prot.dataset[ ,1], ignore.case = T) #get index of the row containing group-names
  intensity.cols <- prot.dataset[gid.index, ] %>% grep("^$", ., invert = T) #intensity columns indexes
  grp.names <- prot.dataset[gid.index, intensity.cols] #obtain group names
  grp.names <- make.names(sub(".*}", "", grp.names, perl = T), unique = T)
  colnames(prot.dataset)[intensity.cols] <- grp.names #replace the intensity column headers with the group names
  major.prot.id <- grep("majority", col.names, ignore.case = T, value = T)
  comment_rows <- grep("^#", prot.dataset[ ,1])
  prot.dataset <- prot.dataset[-(comment_rows), ] #remove comment rows
  
  #if T-test results are present, process relevant columns
  if(isTruthy(grepl("Student", col.names, ignore.case = T))) {
    grp.name.col <- grep(".+significant", col.names, value = T) 
    prot.dataset <- prot.dataset %>% 
      setNames(case_when(
        grepl(".*p.value.*", names(.)) ~ "Minus.log.pval",
        grepl(".*Difference.*", names(.)) ~ "Difference",
        grepl(".*q.value.*", names(.)) ~ "q_val",
        grepl(".*Significant.*", names(.)) ~ "Significant",
        grepl(".*(biological|GOBP).*", names(.), ignore.case = TRUE) ~ "Biological Process",
        grepl(".*(cellular|GOCP).*", names(.), ignore.case = TRUE) ~ "Cellular Component",
        grepl(".*(molecular|GOMF).*", names(.), ignore.case = TRUE) ~ "Molecular Function",
        grepl(".*keywords.*", names(.), ignore.case = TRUE) ~ "Keywords",
        TRUE ~ names(.)
      )) %>%
      mutate_at(vars(Difference, Minus.log.pval, q_val, intensity.cols), as.numeric) %>%
      dplyr::select(all_of(intensity.cols), 
                    all_of(major.prot.id), 
                    all_of(grp.name.col),
                    any_of(c("Minus.log.pval", "q_val", "Difference", "Significant", 
                             "Biological Process", "Cellular Component", "Molecular Function", "Keywords")))
    
    
    #Extract group labels from the "significant" column for use in legend
    sig.df <- prot.dataset[prot.dataset$Significant == "+",]  
    grp.names <- unique(sig.df[, grp.name.col]) #extract label 
    grp.names <- unlist(str_split(grp.names, "_", n = 2)) #separate labels
    new.intensity.names <- names(prot.dataset)[intensity.cols]
    first.group <- startsWith(new.intensity.names, grp.names[1]) %>% new.intensity.names[.]
    second.group <- startsWith(new.intensity.names, grp.names[2]) %>% new.intensity.names[.]
    
    
    #Add the group labels based on positive or negative x axis then calculate mean across groups
    prot.dataset <- prot.dataset %>%
      mutate(Expression = case_when(Difference>=s0 & q_val<=fdr ~ paste0("Higher abundance in ", grp.names[1]),
                                    Difference<=-s0 & q_val<=fdr ~ paste0("Higher abundance in ", grp.names[2]),
                                    TRUE ~ "Not significant")) 
    # if (purrr::is_empty(major.prot.id) == F) {
    #   prot.dataset <- prot.dataset %>% column_to_rownames(var = major.prot.id) 
    # } else { prot.dataset }
    
    #calculate row means as matrix since it's faster
    prot.dataset$`First-Group` <- rowMeans(as.matrix(prot.dataset[,first.group]))
    prot.dataset$`Second-Group` <- rowMeans(as.matrix(prot.dataset[,second.group]))
  } else {  
    prot.dataset
  }
  return(prot.dataset)
}



####2. 1D HEATMAPS ----
#This function accepts the dataframe produced after 1D annotation enrichment in Perseus. The title argument is optional
onedheatmap <- function(oned.df, plot.title = "", fdr = 0.05, score = 1.0) {
  
  oned.df <- oned.df %>%
    mutate_at(c("Score", "Benj. Hoch. FDR"), as.numeric) %>%
    filter(`Benj. Hoch. FDR` <= fdr,
           Score >= -score & Score <= score)
  annotations <- oned.df$Name
  samples <- sort(unique(oned.df$Column))
  
  #generate 1D score matrix
  M.score <- matrix(ncol=length(samples), nrow=length(annotations))
  for(i in 1:length(samples)){
    for(j in 1:length(annotations)){
      score.value <- oned.df$Score[oned.df$Column==samples[i] & oned.df$Name==annotations[j]]
      if(length(score.value)==0){score.value <- NA}
      M.score[j,i] <- score.value
    }
  }
  
  #print there are duplicates if length of score.value > 1
  
  m <- M.score
  rownames(m) <- annotations
  colnames(m) <- samples
  #Assign zero to NA values.
  m[is.na(m)] <- 0
  #collapse matrix for use in ggplot
  m <- melt(m)
  colnames(m) <- c("Annotations", "T-test differences", "value")
  #n.columns = length(unique(m$`T-test differences`))
  #n.rows = length(m$`Annotations`)
  
  return(
    ggplot(data = m, aes(x = `T-test differences`, 
                         y = reorder(Annotations, value), 
                         fill = value))+
      geom_tile(colour = "grey", linewidth = 1)+
      scale_fill_gradientn(colors = c(low = "blue", 
                                      mid = "white",
                                      high = "red"),
                           na.value = "grey",
                           limits = c(-score, score), 
                           breaks = c(-(score/2), 0, score/2))+
      guides(fill = guide_colourbar(title = "Score",
                                    ticks.colour = "black"))+
      theme(axis.title = element_text(face = "bold", size = 11),
            axis.text.x = element_text(angle = 0),
            axis.text = element_text(colour = "black"),
            plot.title = element_text(size = 11),
            text = element_text(size = 12),
            legend.key.size = unit(3, "mm"),
            legend.text = element_text(size=9),
            legend.title = element_text(size = 10),
            legend.position = "right")+
      coord_fixed(ratio = 0.5)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      labs(title = paste0(plot.title, "\n"), x = 'T-test comparison',
           y = "Annotation")
  )
  
}



####3. VOLCANO Plot ----
#The function requires a dataframe that has been processed by process_df() (step 1). It also requires the curves.df exported from Perseus to create the lines on the volcano plot. The other arguments are optional.
#plot.title = custom plot title top left
#group.colors = colours for the default plot in the following order: higher-left, higher-right colours. 
#s0 = s0 value to be shown on the plot. Default = 0.1
#fdr = fdr value to be shown on the plot. Default = 0.05
#fdr.lines = 'TRUE' to show lines and 'FALSE' to remove lines
#num.breaks = number of breaks for the legend on the "average intensity values" plot 

volcano_plot <- function(df, curves.df = NULL, go.category = "Keywords", plot.title = "", 
                         s0 = 1, fdr = 0.05, fdr.lines = "No", group.colors = c("#00B1B2", "#FF6666"),
                         avg.intensity = "No", num.breaks = 3,
                         select.ids = c(), select.go.terms = c(), sig_or_all = "All") {
  
  #Subset significant and non-significant(ns) proteins into their respective dfs. 
  #Extract the number of unique GO terms in both sig and ns dfs. The numbers will be used for color-coding the GO terms downstream. 
  #Subset the selected proteins so that we can use their protein ids as point labels on the plot
  #Also, extract the x-axis label (in "Group1_Group2" format) from "significant" column
  sig.df <- df %>% 
    filter(q_val<=fdr & Difference >= s0 | q_val<=fdr & Difference <= -s0) 
  ns.df <- df %>%
    filter(!(q_val<=fdr & Difference >= s0 | q_val<=fdr & Difference <= -s0))
  col.names <- names(df)
  grp.name.col <- grep(".+significant", col.names, value = T) #column index 
  xlab <- unique(df[df$Significant == "+", grp.name.col]) 
  xlab <- sub("_", " - ", xlab) #for x axis label 
  grp.names <- str_split_1(xlab, " - ") #for legend label
  
  
  
  #Custom settings for ggplots
  a = 0.4 #alpha
  s = 21  #shape
  ns.sz = 2 #non-sig shape size  
  ns.color = "grey" #non-sig colour
  s.sz = 4  #sig shape size
  text.sz = 20
  
  #Base volcano plot to plot points, fdr curves, axes boundaries, and other custom arguments.
  base_vplot <- ggplot(data = ns.df, 
                       aes(x = Difference,
                           
                           y = Minus.log.pval, fill = ns.color))+
    geom_point(size = ns.sz,
               alpha = a, shape = s)+
    labs(title = plot.title,
         caption = paste("s0 =", s0, "  ", "FDR =", fdr),
         x = paste0("Difference (", xlab, ")"), 
         y = expression(-log[10]("p-value")))+
    coord_cartesian(xlim = c(min(df$Difference),
                             max(df$Difference)),
                    ylim = c(min(df$Minus.log.pval),
                             max(df$Minus.log.pval)),
                    expand = T)+
    theme_minimal()+
    theme(plot.title = element_text(face="bold", size = 22),
          plot.caption = element_text(colour = "darkviolet", size = 13),
          axis.line = element_line(colour = "black"),
          text = element_text(size = text.sz))
  
  ###Adding layers to the base plot
  #Layer for protein ID labeling
  if(is.null(select.ids)){
    label.pt.plot <- NULL
  } else {
    major.prot.id <- grep("majority", col.names, ignore.case = T, value = T)
    select.ids <- df[rownames(df) %in% select.ids, ]
    label.pt.plot <- geom_label_repel(data = select.ids,
                                      aes(label = !!sym(major.prot.id)),
                                      nudge_y = 0.2, nudge_x = 0.2, fill = "lightgrey", min.segment.length = 0.05, 
                                      size = 3.5, max.overlaps = getOption("ggrepel.max.overlaps", nrow(select.ids)))
  }
  
  #Layer for labeling GO terms
  if(is.null(select.go.terms)){
    go.terms.plot <- NULL
  } else{
    if (sig_or_all == "All"){
      df <- df
    } else if(sig_or_all == "Significant"){
      df <- sig.df
    }
    hilight.go.terms <- gsub("([][.*+?^$(){}|\\\\])", "\\\\\\1", select.go.terms, perl = TRUE) %>% 
      paste0(., collapse = "|") #escape regex characters in GO terms
    hilight.go.terms <- df %>% 
      filter(grepl(pattern = hilight.go.terms,
                   x = !!sym(go.category))) %>% #select rows containing chosen GO terms 
      mutate(!!sym(go.category) := lapply(str_split(!!sym(go.category), ";"), 
                                          \(x) grep(pattern = hilight.go.terms, x, value = T))) #keep only the chosen go terms per row
    go.terms.plot <- geom_label_repel(data = hilight.go.terms,
                                      aes(label = !!sym(go.category)),
                                      nudge_x = 0.3, nudge_y = 0.3, fill = alpha("lightgreen", 0.7), 
                                      size = 3.5, 
                                      max.overlaps = getOption("ggrepel.max.overlaps", 
                                                                           nrow(hilight.go.terms)
                                      )
    )
  }
  
  #Create default plot
  default.plot <- base_vplot+
    scale_fill_identity(name = "", labels = "Not significant", guide = "legend")+
    new_scale_fill()+
    geom_point(data = sig.df, aes(fill = Expression),
               size = s.sz, shape = s)+
    scale_fill_manual(values = group.colors)+
    guides(fill = guide_legend(order = 2)) +  label.pt.plot + go.terms.plot+
    theme(legend.spacing = unit(-1, "cm")) 
  
  #Create avg intensity plot 
  all.avg.intensities <- c(sig.df$`First-Group`, sig.df$`Second-Group`)
  min.intensity <- floor(min(all.avg.intensities)) #round minimum avg intensity of both groups
  max.intensity <- ceiling(max(all.avg.intensities)) #round maximum avg intensity
  avg.intensity.pl <- base_vplot +
    scale_fill_identity()+
    geom_point(data = sig.df[sig.df$Difference > 0, ],
               aes(size = `First-Group`, fill = alpha(group.colors[2], a)),
               shape = s)+
    geom_point(data = sig.df[sig.df$Difference < 0, ],
               aes(size = `Second-Group`, fill = alpha(group.colors[1], a)),
               shape = s)+
    scale_size_binned(range = c(1,8),
                      name = "Avg LFQ Intensity",
                      n.breaks = num.breaks, limits = c(min.intensity, max.intensity))+
    guides(size = guide_bins(override.aes = list(fill = "darkgray"), show.limits = T))+
    label.pt.plot + go.terms.plot
  
  
  if (avg.intensity == "Yes") {
    #Sizes
    final.plot <- avg.intensity.pl
  } else {
    final.plot <- default.plot
  }
  
  if(fdr.lines == "Yes" && is.null(curves.df)){
    final.plot <- stop("Please upload the matrix containing x and y coordinates for the FDR curves.")
  } else if (fdr.lines == "Yes" && is.null(curves.df) == FALSE) {
    #Fdr curves layer
    curve.plot <- geom_line(data = curves.df, aes(x, y), linetype=2)
    final.plot <- final.plot + curve.plot
    final.plot
  } else if (fdr.lines == "No"){
    final.plot
  }
  
  return(list(final.plot))
}





####4. Venn Diagram ----
#This function takes the matrix from Perseus and replaces the column names with group names. The function returns a narrowed down dataframe with intensity and protein ID columns only
#gid = group identifier
groupnames_to_colnames <- function(venn.df, gid = NULL) {
  col.names <- names(venn.df)
  gid <- ifelse(is.null(gid), "group", gid)
  
  #extract group names then use them as new column names
  gid.index <- grep(gid, venn.df[ ,1], ignore.case = T) #index of the row containing group-names
  intensity.cols <- venn.df[gid.index, ] %>% grep("^$", ., invert = T) #intensity columns indexes
  grp.names <- venn.df[gid.index, intensity.cols] #obtain group names
  
  if (is_empty(intensity.cols) == TRUE)  {
    stop("Cannot detect intensity columns.")
  }
  
  colnames(venn.df)[intensity.cols] <- sub(".*}", "", grp.names, perl = T) #replace the intensity column headers with the group names
  comment_rows <- grep("^#", venn.df[, 1]) #index of commented rows
  venn.df <- venn.df[-(comment_rows), ] #remove commented rows
  venn.df[intensity.cols] <- as.data.frame(lapply(venn.df[intensity.cols], as.numeric),
                                           check.names = F) #ensure that intensities are "numeric"
  id.col <- grep("IDs", names(venn.df), ignore.case = T) #id column
  venn.df <- venn.df[, c(intensity.cols, id.col)]   #keep intensity + protein ID columns
}

#This functions takes the dataframe from the above step and evaluates the presence/absence of proteins based on minimum valid values by group (unique group names as input)
evaluate_valid_vals <- function(venn.df, group.names) {
  min.percent = 75/100
  df_list <- lapply(group.names, function(x) 
    select(venn.df, starts_with(x))
  ) #split the main df into the different groups
  binary.res <- lapply(df_list, \(x) ifelse(
    rowSums(!is.na(x))/ncol(x) >= min.percent,
    1, 0))  #binary results: If % valid values >/= minimum % in each df i.e group, the protein is present and evaluates to 1. Else, it evaluates to 0 for absent proteins
  names(binary.res) <- group.names  #add labels to the binary results 
  #Insert binary results in venn.df
  merged.df <- data.frame(df_list, binary.res, check.names = F)
  
  #Return combined df or binary df only
  output <- list(binary.result = binary.res,
                 merged.df = merged.df)
  return(output)
}


####5. PCA plot ---- 
#This function also accepts the main protein matrix and creates a PCA of all intensity columns.
pca_plot <- function(df, gid = NULL, ellipse = "No", .colours =  NULL) {
  gid <- ifelse(is.null(gid), "group", gid)
  counts <- groupnames_to_colnames(df, gid) %>%
    dplyr::select(!contains("IDs")) %>%
    mutate_all(as.numeric) %>% 
    na.omit() %>% t() %>% as.data.frame() #transpose df of intensity values
  
  pca <- prcomp(counts)
  counts$Group.names <- sub("\\.\\d+", "", rownames(counts)) #remove serial numbers to extract group names
  text.sz = 19
  point.sz = 4.8
  scl = 0
  if (!is.null(.colours)) {
    colour.layer <- scale_colour_manual(values = .colours)
  } else { colour.layer <- NULL}
  layers <- list(
    colour.layer,
    guides(colour=guide_legend("SampleZ", override.aes =list(size=4)), fill = "none", size = "none"),
    theme_classic(),
    theme(axis.line = element_line(colour = "black"),
          text = element_text(size = text.sz),
          panel.border = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "grey"))
  )
  pca_ellipse <- autoplot(pca, data = counts, size = point.sz, 
                          scale = scl, colour = "Group.names", 
                          frame = T, frame.type = "t") + layers
  
  
  pca_plain <- autoplot(pca, data = counts, size = point.sz, 
                        scale = scl, colour = "Group.names") + 
    guides(colour=guide_legend("Sample Type"), fill = "none", size = "none")+ layers
  #extract wider axes limits from pca_ellipse
  pca_ellipse_build <- ggplot_build(pca_ellipse)
  pca_ellipse_x <- pca_ellipse_build$layout$panel_params[[1]]$x.range
  pca_ellipse_y <- pca_ellipse_build$layout$panel_params[[1]]$y.range
  #apply to pca_plain
  pca_plain$coordinates$limits$x <- pca_ellipse_x
  pca_plain$coordinates$limits$y <- pca_ellipse_y
  
  
  if (ellipse == "Yes") {
    pca_ellipse
    
  } else {
    pca_plain
  }
}



####6. S-curve plot ----
#The function accepts the main protein matrix to produce S-curves 
#The function create one s-curve plot for each group and puts them on the same image with ggarrange

scurve_plot <- function(df, gid = NULL, y.increment = 1, .colours = NULL) {
  gid <- ifelse(is.null(gid), "group", gid)
  df <- groupnames_to_colnames(df, gid) %>%
    dplyr::select(!contains("IDs")) %>%
    mutate_if(is.numeric, ~2^(.)) %>% mutate_if(is.numeric, ~log10(.)) 
  all.groups <- names(df)
  all.groups <- gsub("\\.\\d+", "", all.groups) #remove serial numbers to extract group names
  uniq.groups <- unique(all.groups)
  groupmeans <- mapply(function(each.grp)
    rowMeans(as.matrix(df[each.grp])),
    uniq.groups)
  colnames(groupmeans) <- paste0(colnames(groupmeans), "_means")
  df <- as.data.frame(cbind(df, groupmeans), check.names = F)
  col.names <- names(df)
  text.sz = 13
  
  if (!is.null(.colours)) {
    .colours <- .colours
  } else { 
    .colours <- hcl.colors(length(uniq.groups))
  }
  ranked.df <- df %>%
    mutate(across(ends_with("means"), ~ rank(dplyr::desc(.)), .names = "{.col}_ranked"))
  
  
  
  #function to create plot
  plot.s.range <- function(grp.name, my.col) {
    y.name <- paste0(grp.name, "_means")
    y.max <- max(ranked.df[y.name])
    ggplot(data = ranked.df,
           aes(x = !!sym(paste0(grp.name, "_means_ranked")),
               y = !!sym(paste0(grp.name, "_means"))))+
      geom_point(size = 2, colour = my.col, shape = 1)+
      labs(x = paste0('Protein Rank ', '(',grp.name, ')'),
           y = bquote(log[10] * ' ' * 'LFQ Intensity (' * .(grp.name) * ')'))+
      scale_y_continuous(breaks = seq(0, y.max, by = y.increment))+
      theme_minimal()+
      theme(axis.line = element_line(colour = "black"),
            text = element_text(size = text.sz))
  }
  plot.list <- map2(uniq.groups, .colours, plot.s.range) #create all plots
  
  return(
    #plot s.range plots
    ggarrange(plotlist = plot.list)
  )
}  







####7. SHINY ----
#It accepts the protein matrix from Perseus and processes it using step 1 for the remaining program. It also accepts the 1D matrix from Perseus to create 1D heatmaps and to filter out less significant GO terms on the volcano plot. Thirdly, it accepts the curve matrix from Perseus for the lines on the volcano plot 
#A Shiny app's function has two sections as seen below: the UI and the Server

dl.unit = "px"

#Shiny UI ---- 
ui <- shinyUI(fluidPage(
  titlePanel(h1(div(img(src = "ProteoPlotter.png")), align = "center"), windowTitle = 'ProteoPlotter'),
  tags$style(HTML(".shiny-output-error-validation {color: red; font-size:18px;}")), #change error messages to colour red
  
  fluidRow(
    column(2, offset = 8,
           tags$style(HTML(".navbar {left:80px;}")),
           navbarPage(title = h6(icon("images"), "Dimensions"), 
                      navbarMenu("",
                                 tabPanel(
                                   numericInput("width.val", label = "Width (px)", value = 1920, width = '70%')
                                 ),
                                 tabPanel(
                                   numericInput("height.val", label = "Height (px)", value = 1080, width = '70%')
                                 )
                      ),
           )
    ),
    column(2,
           navbarPage(title = h6(icon("images"), "DPI"),
                      navbarMenu("",
                                 tabPanel(
                                   numericInput("dpi.val", label = "DPI", value = 300, width = '70%')
                                 )
                      )
           )
    )
  ),
  tabsetPanel(
    tags$style(HTML(".homepg { margin-top: -2px; }")),
    tabPanel(h2(class= "homepg", 
                icon("home"), style = "font-weight: bold; color: black; background-color: #EDEBDAFC;"),
             fluidRow(
               column(4),
               column(4, style = "background-color: cornsilk;",
                      h4("ProteoPlotter empowers users to visualize proteomics data after quality control and data analysis in Perseus.", br(), br(),
                         div(img(src = "ProteoPlotter_workflow.png"), style="text-align: center;"),
                         br(), br()),
                      p("The source code and sample files are available on", 
                        a("Github.", href = "https://https://github.com/JGM-Lab-UoG/ProteoPlotter", target = "_blank"), style = "color: darkred;"),
                      p("Contact us: jgeddesm@uoguelph.ca", style = "color: darkred;"),
                      br(), br(), br()),
               column(4)
             )),
    #Panel for 1D Heatmap
    tabPanel("1D Annotation Heatmap",
             sidebarLayout(
               sidebarPanel(width = 4,
                            fluidRow(
                              column(width = 12,
                                     fileInput("onedfile", label = "1D Annotation Enrichment Matrix")),
                              column(width = 12, selectizeInput("one.d.order", 
                                                                label = "Order Of Groups:", 
                                                                choices = c(), multiple = T)),
                              column(width = 6, numericInput("change.fdr",
                                                             label = "FDR", 
                                                             value = 0.05, width = '50%', step = 0.01)),
                              column(width = 6, numericInput("change.score",
                                                             label = "Min/Max 1D Score", 
                                                             value = 1.0)),
                              p("Angle of x-axis label:", 
                                style="text-align: left; margin-left: 12px;"),
                              column(width = 3, numericInput("x.angle.bp",
                                                             label = "GOBP", 
                                                             value = 0, min = 0, max = 360)),
                              column(width = 3, numericInput("x.angle.cc",
                                                             label = "GOCC", 
                                                             value = 0, min = 0, max = 360)),
                              column(width = 3, numericInput("x.angle.mf",
                                                             label = "GOMF", 
                                                             value = 0, min = 0, max = 360)),
                              column(width = 3, numericInput("x.angle.kw",
                                                             label = "KW", 
                                                             value = 0, min = 0, max = 360))
                            )
               ),
               mainPanel(
                 br(),
                 br(),
                 ##1D plots output
                 h4("1D Annotation Heatmaps", align = "center", ),
                 
                 plotOutput("onedhm1.bio.proc"),
                 downloadBttn("dlhm1", "Heatmap 1", size = "xs"),
                 br(),
                 br(),
                 plotOutput("onedhm2.cell.comp"),
                 downloadBttn("dlhm2", "Heatmap 2", size = "xs"),
                 br(),
                 br(),
                 plotOutput("onedhm3.mol.func"),
                 downloadBttn("dlhm3", "Heatmap 3", size = "xs"),
                 br(),
                 br(),
                 plotOutput("onedhm4.keywords"),
                 downloadBttn("dlhm4", "Heatmap 4", size = "xs"), 
               )
             )
    ),
    
    #Main Tab for Volcano plot, Scurve, and PCA
    tabPanel("Volcano Plot", 
             sidebarLayout(
               sidebarPanel(width = 4,
                            fluidRow(
                              column(width = 6,
                                     fileInput(inputId = "proteinfile", label = "Matrix With T-test *") ),
                              column(width = 6, 
                                     fileInput("curvesfile", label = "FDR Curves Matrix")),
                              column(width = 6, actionButton("gen.vplot", label = "Generate Plot")),
                              column(width = 12, br())
                            ),
                            
                            useShinyjs(),
                            fluidRow(
                              div(id = "reset.grp",
                                  h4("Plot Options"), 
                                  column(width = 3, numericInput("xmin",
                                                                 label = "x min", 
                                                                 value = NA) ),
                                  column(width = 3, numericInput("xmax", 
                                                                 label = "x max", 
                                                                 value = NA) ),
                                  column(width = 3, numericInput("ymin", 
                                                                 label = "y min", 
                                                                 value = NA) ),
                                  column(width = 3, numericInput("ymax", 
                                                                 label = "y max", 
                                                                 value = NA) ),
                                  column(width = 6, colourpicker::colourInput("left.col",
                                                                              label = "Left Group",
                                                                              value = "#00B1B2") ),
                                  column(width = 6, colourpicker::colourInput("right.col",
                                                                              label = "Right Group",
                                                                              value = "#FF6666") ),
                                  column(width = 12, br()),
                                  column(width = 6, numericInput("s.knot",
                                                                 label = "s0 Value",
                                                                 value = 1, step = 0.5)),
                                  column(width = 6, numericInput("fdr.val",
                                                                 label = "FDR",
                                                                 value = 0.05, step = 0.01)),
                                  column(width = 12, br()),
                                  column(width = 4, radioButtons("fdrlines",
                                                                 label = "Display Curves",
                                                                 choices = list("Yes", "No"),
                                                                 selected = "No", inline = T) ),
                                  column(width = 5, radioButtons("show.avg",
                                                                 label = "Average Protein Intensities",
                                                                 choices = list("Yes", "No"),
                                                                 selected = "No", inline = T) ),
                                  column(width = 3, numericInput("legend.brks",
                                                                 label = "# Legend Breaks",
                                                                 value = 3, step = 1)),
                                  column(width = 12, h5("Gene Ontology Terms", style = "color: darkred;")),
                                  column(width = 12, selectInput("go.categories",
                                                                 label = "Hover over the plot to display:",
                                                                 choices = list("Keywords", "Cellular Component", "Molecular Function", "Biological Process"),
                                                                 selected = "Keywords",
                                                                 width = '100%')),
                                  column(width = 12, selectInput("go.terms",
                                                                 label = "Annotation terms:",
                                                                 choices = list(), multiple = T)),
                                  column(width = 12, radioButtons("which.proteins",
                                                                  label = "All vs. significant proteins",
                                                                  choices = list("All", "Significant"),
                                                                  selected = "Significant", inline = T, width = '100%')),
                                  column(width = 12, checkboxInput("filter.terms.by.1d",
                                                                   label = HTML("Keep only GO terms from the 1D enrichment data<br/>(Requires the 1D matrix*)"))),
                                  column(width = 12, br()),
                                  column(width = 12, actionButton("reset", label = "Reset Volcano Plot")),
                              )
                            )
               ),
               mainPanel(
                 br(), br(), br(),
                 
                 h4("Volcano Plot",  align = "center"),
                 ##Volcano plot output with brush selection of points
                 uiOutput("curves.err"),
                 plotOutput("volcanoplot",
                            brush = brushOpts("vplot.brush"), 
                            hover = hoverOpts("vplot.hover")),
                 h6("Click and drag cursor to highlight protein IDs.", style = "color: darkred;"),
                 column(width=4,
                        actionButton("deselect.ids", label = "Deselect IDs")),
                 br(), br(),
                 verbatimTextOutput("hover.info"),
                 br(), br(), br(),
                 column(width=4, 
                        textInput("vplot.title",
                                  label = "Volcano Plot Title (Optional)",
                                  value = NULL,
                                  placeholder = "Enter text...")),
                 column(width = 4, 
                        downloadBttn("dlvplot", "Volcano plot", size = "xs")),
                 br(), br(), br(),br(),
               )
             )
    ),
    
    #PCA Tab
    tabPanel("PCA Plot",
             sidebarLayout(
               sidebarPanel(
                 # br(),
                 column(width = 8, p("Upload matrix under the Volcano Plot tab")),
                 
                 column(width = 4, actionButton("gen.pca", label = "Generate Plot")),
                 column(12, br()),
                 column(width = 12, h6("Group identifier: based on grouping levels in the uploaded matrix.")),
                 column(width = 6, selectizeInput("grp.identifier", 
                                                  label = "Group Identifier",
                                                  choices = c())),
                 column(width = 4, radioButtons("ellipse", label = "Ellipse",
                                                choices = list("Yes", "No"), selected = "No",
                                                inline = T), offset = 2),
                 column(width = 12, uiOutput("pca.colour.picker")),
                 
                 br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br()
               ),
               mainPanel(
                 #PCA plot output
                 h4("PCA Plot", align = "center"),
                 plotOutput("pcaplot", height = "500px"),
                 downloadBttn("dlpca", "PCA", size = "xs")
               )
             )),
    
    tabPanel("Dynamic Range Plot",
             sidebarLayout(
               sidebarPanel(
                 column(width = 12, p("Upload matrix under the Volcano Plot tab")),
                 column(width = 12, actionButton("gen.scurve", label = "Generate Plot")), 
                 column(12, br()), br(), br(), br(),
                 column(width = 5, numericInput("y.incr", label = "Increment (y-axis)", value = 1)),
                 column(width = 12, uiOutput("scurve.colour.picker")),
                 br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br()
                 
                 
               ),
               mainPanel(
                 ##S-curve plot output
                 h4("Dynamic Range Plot",  align = "center"),
                 plotOutput("scurveplot"),
                 downloadBttn("dlscurve", "Dynamic Range Plot", size = "xs"),
                 br(), br()
               )
             )),
    #Venn diagram tab
    tabPanel("Venn Diagram", 
             sidebarLayout(
               sidebarPanel(width = 4,
                            fluidRow(
                              column(width = 12,
                                     fileInput("vennfile", label = "Venn Diagram / UpSet Plot Matrix")),
                              column(width = 12, h6("Group identifier: based on grouping levels in the uploaded matrix.")),
                              column(width = 6, selectizeInput("grp.identifier.v", 
                                                               label = "Group Identifier", 
                                                               choices = c(), width = "80%")),
                              column(width = 6, actionButton("gen.venn", label = "Generate Venn", 
                                                             style="color: blue;")),
                              column(width = 12, selectInput("venngroups", 
                                                             label = "Select Groups:", 
                                                             choices = c(), multiple = T)),
                              shinyBS::bsTooltip(id = "venngroups", title = "Select groups to be compared.", "focus"),
                              column(width = 12, br()),
                              column(width = 12, uiOutput("col.title.venn")),
                              column(width = 12, uiOutput("venn.colour.picker")),
                              column(width = 12, br()),
                              column(width = 12, checkboxGroupInput("count_percent", 
                                                                    label = "Counts / Percent",
                                                                    choices = c("counts", "percent"), 
                                                                    selected = "counts",
                                                                    inline = T )),
                              column(width = 3, selectInput("linetype", 
                                                            label = "Stroke (type)",
                                                            choices = c("blank", "solid", "dashed", "dotted", "longdash"),
                                                            selected = "solid")),
                              column(width = 3, numericInput("linewidth", 
                                                             label = "Stroke (width)",
                                                             min = 0, max = 10, value = 1)),
                              column(width = 3, numericInput("labelsize", 
                                                             label = "Size (labels)",
                                                             min = 4, value = 12)),
                              column(width = 3, numericInput("numsize", 
                                                             label = "Size (values)",
                                                             min = 3, value = 10))
                            )
                            
               ),
               mainPanel(
                 tags$head(
                   tags$style(
                     HTML("#venndf\\.preview 
                                     {font-size: 15px; color: red; font-family: 'Courier New';}")
                   )),
                 htmlOutput("venndf.preview", align="center"),
                 br(),
                 DTOutput("venntable"),
                 br(), br(), br(), br(), br(),
                 plotOutput("vennplot"),
                 downloadBttn("dlvennplot", "Venn Plot", size = "xs")
               )
             )
    ),
    #UpSet plot tab
    tabPanel("UpSet Plot", 
             sidebarLayout(
               sidebarPanel(width = 4,
                            fluidRow(
                              column(width = 8, p("Upload matrix under the Venn Diagram tab")),
                              column(width = 12, selectInput("upsetgroups", 
                                                             label = "Select Groups:", 
                                                             choices = c(), multiple = T)),
                              column(width = 6, actionButton("gen.upset", label = "Generate Plot")),
                              column(width = 12, br()),
                              column(width = 4, numericInput("upset.text.sz", label = 'Text size', 
                                                             value = 1.5, step = 0.5, max = 10)),
                              column(width = 4, numericInput("upset.pt.sz", label = 'Dot size', 
                                                             value = 2.5, step = 0.2, max = 10)),
                              column(width = 3, numericInput("upset.line.sz", label = 'Line size', 
                                                             value = 1, step = 0.2, max = 10)),
                              column(width = 12, uiOutput("col.title.upSet")),
                              column(width = 12, uiOutput(outputId = "upSet.colour.picker"))
                            )
               ),
               mainPanel(
                 plotOutput("upset.plot"),
                 downloadBttn("dlupsetplot", "UpSet Plot", size = "xs")
               ))
    ),
    tabPanel(p(icon("circle-info"), "Guide", style = "font-weight: bold; color: #A10000FE;"),
             fluidRow(style = "background-color: #DEDEDEFD;",
                      column(12, 
                             p(br(), 
                               "Note: All matrices are to be exported from Perseus in the .txt format. The 'Dimensions' and DPI dropdown menu on the upper-right enable the user to modify the download width, height, and resolution.")),
                      column(4,
                             p(h5("1D Annotation Heatmap",  style="font-weight:bold; color:#000045FF;"),
                               h5("- Visualizes the functional enrichment of annotation terms within samples, based on 1D annotation enrichment scores ranging from -1 to 1.", style="color:#000045FF;"),
                               tags$b("Input:"), "Matrix containing 1D annotation scores computed with the '1D annotation enrichment' feature in Perseus - exported in the .txt format.", br(),
                               tags$b("Create:"), "Upload the matrix to the file-input field within the '1D Annotation Heatmap' tab to generate 1D heatmaps for Gene Ontology (GO) categories and Keywords. The customization options allow you to reorder the heatmap columns. You can also filter the 1D annotation results to retain those under a certain FDR threshold or within a specific range of the 1D score i.e. -1D Score : +1D Score.", br(), br(),
                               
                               h5("Volcano Plot (only reads Student's T-test input data)",  style="font-weight:bold; color:#000045FF;"),
                               h5("- Illustrates proteins that significantly differ in abundance between compared proteomes.", style="color:#000045FF;"),
                               tags$b("Input 1:"), "Imputed matrix containing T-test result columns for only one pairwise comparison and the other columns for protein intensities, Protein IDs, etc. When computing the T-test in Perseus, please ensure that q-values are reported in the matrix.", br(),
                               "GO annotation is not required for the volcano plot matrix.", br(),
                               tags$b("Input 2:"), "Optional. You can add upload a second matrix which contains x-y coordinates for the curves indicating the FDR threshold on the volcano plot. This matrix can be exported after generating a volcano plot for the group pair within Perseus.", br(),
                               tags$b("Create:"), "Upload the required matrix to the file-input field within the 'Volcano Plot' tab. Select the group identifier based on the commented rows in the matrix, then generate the volcano plot.", br(),
                               "To display GO annotation for proteins, the matrix must include GO annotation columns. You can change colours and filter proteins by FDR and s0 values using the customization options.", br(), br()
                             )),
                      
                      column(4, 
                             p(h5("PCA and Dynamic Range Plots",  style="font-weight:bold; color:#000045FF;"),
                               h5("- The PCA plot reduces the dimensionality of proteins into their principal components, allowing the user to identify patterns and degree of variation among replicates or groups. The dynamic range plot illustrates the distribution of protein abundance across proteomes.",
                                  style="color:#000045FF;"),
                               tags$b("Input:"), "The volcano plot matrix is used to generate both the PCA plot and Dynamic range plot.", br(),
                               "T-test result columns are not required for these plots so the imputed matrix can be exported without calculating T-test results.", br(),
                               tags$b("Create:"), "Simply upload the matrix within the 'Volcano Plot' tab, select the appropriate group identifier, then return to the 'PCA Plot' or 'Dynamic Range Plot' tab to generate the figures.", br(), 
                               "For the PCA plot, you can add ellipses to show the probability distribution of observations within a 95% CI.", br(),
                               "*The PCA plot allows ellipses only when there are 4 or more groups.", br(), br(),
                               
                               h5("Venn Diagram",  style="font-weight:bold; color:#000045FF;"),
                               h5("- Compares the number of proteins across up to 5 sample groups, revealing their degree of similarity and difference.",  style="color:#000045FF;"),
                               tags$b("Input:"), "Accepts the matrix containing average protein intensity values computed using the 'average groups' tab in Perseus.", br(),
                               tags$b("Create:"), "Upload the matrix within the 'Venn Diagram' tab, select the appropriate identifier, then generate the Venn diagram.", br(),
                               "Using the 'Set Group' customization option, you can choose or reorder the groups to compare within the figure.", br(),
                               "*See the UpSet plot feature to visualize more than 5 groups."
                             )),
                      
                      column(4,
                             p(h5("UpSet Plot",  style="font-weight:bold; color:#000045FF;"),
                               h5("- Similar to the Venn diagram, this compares the number of proteins across up to 25 sample groups.",  style="color:#000045FF;"),
                               tags$b("Input:"), "Utilizes the same matrix as the Venn diagram.", br(),
                               tags$b("Create:"), "Upload the matrix within the 'Venn Diagram' tab, select the identifier, then generate the UpSet plot within the 'UpSet Plot' tab, choosing or reordering groups as needed.", br(), br(), br() 
                             ))
             )
    )
  )
)  
)




#Shiny Server ----
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2) #increase input size limit to 30MB
  #Read files -----
  #Protein input file
  df.prot <- eventReactive(input$proteinfile, {
    f.path  <- input$proteinfile$datapath
    validate(need(file_ext(f.path) == "txt", "Cannot detect a .txt file. Please check."))
    read.delim(f.path, na.strings = c("NA", "NaN"), check.names = F, quote = "")
  })
  #Populate grp.identifier
  observe({req(df.prot())
    tryCatch({
      identifiers <- extract_identifiers(df.prot())
      updateSelectizeInput(session, "grp.identifier",
                           choices = identifiers,
                           server = T)
    })
    # print(is.null(identifiers))
    # },
    # warning = function(w) {
    #   showNotification(w$message, '', duration = NULL, type = "error")
    #   return()
    # }, 
    # error = function(e) {
    #   showNotification(e$message, '', duration = NULL, type = "error")
    #   return()
    # }, 
    # silent=TRUE )
  })
  #Curves file
  df.curves <- eventReactive(input$curvesfile, {
    f.path <- input$curvesfile$datapath
    validate(need(file_ext(f.path) == "txt", "Cannot detect a .txt file for 'FDR Curves Matrix'. Please check.")) #Ensure that input is a txt file
    df <- read.delim(f.path, header = T, nrows = 1) #If txt file, read the df but only the first row
    headers <- c("x", "y") #define headers to look for below
    validate(need(headers %in% colnames(df), "Please upload the matrix containing x and y coordinates for the FDR curves.")) #Ensure that column names in the input file match the defined headers (x and y)
    df <- read.table(f.path, header = T) #If column names are x and y, read the full df
  })
  #1D annotation file
  df.1d <- eventReactive(input$onedfile, {
    f.path <- input$onedfile$datapath
    validate(need(file_ext(f.path) == "txt", "Cannot detect a .txt file. Please check."))
    df <- read.delim(f.path, comment.char = "#", check.names = F, quote = "", nrows = 1) 
    headers <- c("Column", "Type", "Name", "Score")
    validate(need(headers %in% colnames(df), "At least one of the following columns is missing within the 1D matrix:
                  'Column', 'Type', 'Name', 'Score'"))
    df <- read.delim(f.path, comment.char = "#", check.names = F, quote = "") %>%
      mutate(Column = sub(".*Difference ", "", Column)) #shorten the content of the column called 'Column' to retain the names of groups only
  })
  #update UI with 1d group names
  observe({df.1d()
    tryCatch({
      updateSelectizeInput(session,
                           "one.d.order",
                           choices = unique(df.1d()$Column),
                           selected = unique(df.1d()$Column),
                           server = T)
    })
  })
  #Reorder 1d groups based on UI selection
  df.1d2 <- reactive({req(df.1d(), input$one.d.order)
    tryCatch({
      df <- df.1d()
      df %>%
        filter(Column %in% input$one.d.order) %>% #keep selected groups
        mutate(Column = factor(Column, levels = input$one.d.order)) #order 
    })
  })
  #Process the protein matrix 
  df.prot2 <- reactive({
    req(df.prot())
    tryCatch({
      if (length(input$grp.identifier) > 1){ 
        grp.identifier = input$grp.identifier
      } else {
        grp.identifier = NULL
      }
      df <- process_df(prot.dataset = df.prot(), fdr = input$fdr.val, s0 = input$s.knot, gid = grp.identifier)
    })
  })

  
  
  ####Volcano plot -----  
  #hover to show all GO Terms per point
  hovered.pts <- reactive({req(input$go.categories, input$vplot.hover)
    pts <- nearPoints(df.prot2(), input$vplot.hover, threshold = 10)
    pts[[input$go.categories]] 
    # current.pts <- brushed.pts() #access reactive val
    # current.pts <- union(current.pts, ids) #return every selected point to the reactive value
    # brushed.pts(current.pts)
    # brushed.pts()
  })
  #display hover
  output$hover.info <- renderPrint({ 
    tryCatch({
      pts <- hovered.pts()
      if (!is_empty(pts)){
        pts <- gsub(";", ", ", pts)
        pts <- str_wrap(pts, width = 60)
        cat(input$go.categories, ":\n", pts, sep = "") 
      } 
    })
  })
  #brush to highlight protein IDs
  brushed.pts <- reactiveVal()
  observeEvent(input$vplot.brush, {
    tryCatch({
      pts <- brushedPoints(df.prot2(), input$vplot.brush) #select brushed points as rows of df.prot
      # major.prot.id <- grep("majority", names(df.prot2()), ignore.case = T, value = T)
      ids <- rownames(pts) #extract protein row indices
      current.pts <- brushed.pts() #access reactive val
      current.pts <- union(current.pts, ids) #return every selected point to the reactive value
      brushed.pts(current.pts)
      brushed.pts()
    })
  })
  observeEvent(input$deselect.ids, {
    brushed.pts(NULL) #clear IDs highlighted on the plot
    session$resetBrush("vplot.brush") #clear frame brushed onto the plot
  })
  #Populate list of go terms to UI
  observe({req(input$go.categories, df.prot2())
    tryCatch({
      go.category <- input$go.categories
      df.prot2 <- df.prot2()
      if(go.category %in% names(df.prot2)){
        all.go.terms <- df.prot2[, go.category] %>% str_split(., ";") %>% ifelse(. == "", NA, .) %>% 
          unlist() %>% unique() %>% sort(na.last = T) #extract unique GO terms & sort so NA is last
        
        
        if (input$filter.terms.by.1d == T) {
          all.go.terms <- all.go.terms[all.go.terms %in% df.1d()$Name]
        } 
        
        updateSelectInput(session,
                          "go.terms",
                          choices = c("Select from dropdown"='', all.go.terms))
      }
      
    })
  })
  
  
  vplot <- reactive({
    req(df.prot2())
    req.pattern <- grepl("Student", x = names(df.prot2()))
    validate(
      need(isTruthy(req.pattern), 
           "Cannot detect Student's T-test columns for the volcano plot within the uploaded matrix. Please check."))
    req(input$gen.vplot, input$left.col, input$right.col)
    # curves.df <- df.curves()
    # print(curves.df)
    tryCatch({
      vplot <- volcano_plot(df = df.prot2(),
                            curves.df = df.curves(),
                            plot.title = input$vplot.title,
                            go.category = input$go.categories,
                            s0 = input$s.knot,
                            fdr = input$fdr.val,
                            group.colors = c(input$left.col, input$right.col),
                            fdr.lines = input$fdrlines,
                            avg.intensity = input$show.avg, 
                            num.breaks = input$legend.brks, 
                            select.ids = c(brushed.pts()), 
                            select.go.terms = input$go.terms,
                            sig_or_all = input$which.proteins)
    })
  })
  #send vplot to the UI
  output$volcanoplot <- renderPlot({
    tryCatch({
      if (is.na(input$xmin) && is.na(input$xmax) && is.na(input$ymin) && is.na(input$ymax)) {
        vplot()[[1]]
      } else {
        #Axes limits adjustment
        vplot()[[1]]+
          coord_cartesian(xlim = c(input$xmin, input$xmax),
                          ylim = c(input$ymin, input$ymax),
                          expand = T)
      }
    })
  })
  #Display warning for overlaps if annotating many long terms
  observe({req(input$go.terms, vplot())
    pl.b <- ggplot_build(vplot()[[1]])
    #Extract annotation data from ggplot
    ann.terms <- pl.b$data[[3]]
    n.terms = nrow(ann.terms) #how many annotated terms 
    n.char = max(sapply(input$go.terms, nchar)) #num characters per selected term
    if(n.terms >= 25 && n.char >= 25) {
      showNotification({ 
        h4("Warning: When a lengthy annotation term appears across many proteins, some annotations may be obscured due to an increase in overlap.",  style = "color: red; font-size: 13px;") }, duration = NULL)
    }
  })
  #Add error message for when FDR curve file is missing.
  observe({
    fdr.opts <- input$fdrlines
    curve.file <- input$curvesfile
    if(fdr.opts == "Yes" && is.null(curve.file)){
      output$curves.err <- renderUI({ h4("Please upload the matrix containing x and y coordinates for the FDR curves.",  style = "color: red;") })
    } else if(fdr.opts == "No" | !is.null(curve.file)){  
      output$curves.err  <- renderUI({ NULL }
      )}
  })
  #Reset volcano to default plot
  observeEvent(input$reset, {
    shinyjs::reset("reset.grp")
  })
  #Download volcano plot
  output$dlvplot <- downloadHandler( 
    filename = "Volcano_plot.png",
    content = function(file) {
      ggsave(file, plot = vplot()[[1]],
             width = input$width.val, height = input$height.val, dpi = input$dpi.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  
  
  
  
  ####PCA plot -----
  pca.groups <- reactive({req(input$grp.identifier, df.prot()) 
    tryCatch({
      df <- groupnames_to_colnames(df.prot(), gid = input$grp.identifier) %>% dplyr::select(!contains("IDs")) #add group names as new header
      #Extract group names from df. Also, Shiny does not like non-alphanumeric symbols (\\W) in input ids so remove if any
      grp.names <- gsub("\\.\\d+", "", names(df)) %>% unique(.) %>% 
        gsub("\\W", "_", ., perl = T) 
    })
  })
  observeEvent(pca.groups(), {
    pca.groups <- pca.groups()
    output$pca.colour.picker <- renderUI({
      lapply(pca.groups, function(p) {
        column(width = 4, colourpicker::colourInput(inputId = p, label = paste0(p),
                                                    value = "grey21", allowTransparent = T))
      })
    })
    output$scurve.colour.picker <- renderUI({
      lapply(pca.groups, function(s) {
        formatted.s <- paste0("s", s) #add s in front of scurve plot colour ids to differentiate from pca colour ids.
        column(width = 4, colourpicker::colourInput(inputId = formatted.s, label = paste0(s),
                                                    value = "grey21", allowTransparent = T))
      })
    })
  })
  pca.colour.ids <- reactive({req(pca.groups()) #put all colour ids together
    sapply(pca.groups(), function(p) { 
      input[[p]]
    })
  })
  scurve.colour.ids <- reactive({req(pca.groups()) 
    sapply(pca.groups(), function(s) { 
      formatted.s <- paste0("s", s)
      input[[formatted.s]]
    })
  })
  pca.pl <- reactive({req(input$gen.pca, input$grp.identifier, df.prot()) 
    tryCatch({
      pca_plot(df.prot(), gid = input$grp.identifier, ellipse = input$ellipse, .colours = pca.colour.ids())
    })
  })
  output$pcaplot <- renderPlot({
    pca.pl()
  })
  output$dlpca <- downloadHandler( 
    filename = "PCA.png",
    content = function(file) {
      ggsave(file, plot = pca.pl(),
             width = input$width.val, height = input$height.val, dpi = input$dpi.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  ####S-curve plot -----
  scurve.pl <- reactive({
    req(input$gen.scurve, input$grp.identifier)
    tryCatch({
      scurve_plot(df.prot(), gid = input$grp.identifier, y.increment = input$y.incr, .colours = scurve.colour.ids()) 
    })
  })
  output$scurveplot <- renderPlot({
    scurve.pl() 
  })
  #download
  output$dlscurve <- downloadHandler( 
    filename = function() { paste("scurve.png") },
    content = function(file) {
      ggsave(file, plot = scurve.pl(),
             width = input$width.val, height = input$height.val, dpi = input$dpi.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  
  
  
  
  ####Heatmap plots -----
  #create heatmap list
  plotlist.1d <- reactive({
    tryCatch({
      unique.annot.types <- df.1d2() %>%
        group_by(Type) %>%
        group_split()
      lapply(unique.annot.types, 
             function(x) onedheatmap(x, 
                                     plot.title = paste(unique(x$Type)), 
                                     fdr = input$change.fdr,
                                     score = input$change.score
             ))
    })
  })
  #plot heatmaps individually
  bio.proc <- reactive(plotlist.1d()[[1]]
                       + theme(axis.text.x = element_text(angle = input$x.angle.bp, 
                                                          hjust = ifelse(input$x.angle.bp > 0, 1, 0.5))
                       ))
  cell.comp <- reactive(plotlist.1d()[[2]]
                        + theme(axis.text.x = element_text(angle = input$x.angle.cc, 
                                                           hjust = ifelse(input$x.angle.cc > 0, 1, 0.5))
                        ))
  mol.func <- reactive(plotlist.1d()[[3]]
                       + theme(axis.text.x = element_text(angle = input$x.angle.mf, 
                                                          hjust = ifelse(input$x.angle.mf > 0, 1, 0.5))
                       ))
  keywords <- reactive(plotlist.1d()[[4]]
                       + theme(axis.text.x = element_text(angle = input$x.angle.kw, 
                                                          hjust = ifelse(input$x.angle.kw > 0, 1, 0.5))
                       ))
  output$onedhm1.bio.proc <- renderPlot({ 
    bio.proc()
  }) 
  output$onedhm2.cell.comp <- renderPlot({ 
    cell.comp()
  }) 
  output$onedhm3.mol.func <- renderPlot({ 
    mol.func()
  }) 
  output$onedhm4.keywords <- renderPlot({ 
    keywords()
  })
  #downloads
  output$dlhm1 <- downloadHandler( 
    filename = "Bio_process.png",
    content = function(file) {
      ggsave(file, plot = bio.proc(), 
             width = input$width.val, height = input$height.val, dpi = input$dpi.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  output$dlhm2 <- downloadHandler( 
    filename = "Cell_component.png",
    content = function(file) {
      ggsave(file, plot = cell.comp(),
             width = input$width.val, height = input$height.val, dpi = input$dpi.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  output$dlhm3 <- downloadHandler( 
    filename = "Mol_function.png",
    content = function(file) {
      ggsave(file, plot = mol.func(),
             width = input$width.val, height = input$height.val, dpi = input$dpi.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  output$dlhm4 <- downloadHandler( 
    filename = "Keywords.png",
    content = function(file) {
      ggsave(file, plot = keywords(), 
             weight = input$width.val, height = input$height.val, dpi = input$dpi.val,
             device = "png", bg = "white", units = dl.unit)
    }
  )
  
  
  
  
  ####Venn Diagram and UpSet Plot-----
  #Read venn matrix 
  df.venn <- eventReactive(input$vennfile, {
    f.path = input$vennfile$datapath #file path
    validate(need(file_ext(f.path) == "txt", "Cannot detect a .txt file. Please check."))
    read.delim(f.path, check.names = F, quote = "")
  })
  observe({req(df.venn())
    tryCatch({
      identifiers.v <- extract_identifiers(df.venn())
      updateSelectizeInput(session, "grp.identifier.v",
                           choices = identifiers.v,
                           server = T)
    })
  })
  
  observe({req(input$grp.identifier.v)
    tryCatch({
      if(!is.null(input$grp.identifier.v)){
        df <- groupnames_to_colnames(df.venn(), gid = input$grp.identifier.v)
        df <- select(df, !contains("IDs"))
        venn.groups <- names(df)
        #Show group names on UI
        updateSelectInput(session,
                          inputId = "venngroups",
                          choices = venn.groups,
                          selected = venn.groups)
        updateSelectInput(session,
                          inputId = "upsetgroups",
                          choices = venn.groups,
                          selected = venn.groups)
      }
    })
  })
  
  #Convert data to binary
  venn.input <- reactive({req(df.venn(), input$venngroups)
    tryCatch({
      df <- groupnames_to_colnames(df.venn(), gid = input$grp.identifier.v)
      evaluate_valid_vals(venn.df = df, 
                          group.names = input$venngroups)
    })
  })
  upset.input <- reactive({req(df.venn(), input$upsetgroups)
    tryCatch({
      df <- groupnames_to_colnames(df.venn(), gid = input$grp.identifier.v)
      evaluate_valid_vals(venn.df = df, 
                          group.names = input$upsetgroups)
    })
  })
 
  observe({req(input$venngroups)
    venn.groups <- input$venngroups
    output$col.title.venn <- renderUI({ 
      p("Colours:", style = "font-weight: bold; color: #8b0000; background-color: lightgrey;")
    })
    output$venn.colour.picker <- renderUI({
      lapply(venn.groups, function(v) {
        formatted.v <- gsub("\\W", "_", v, perl = T)
        column(width = 4, colourpicker::colourInput(inputId = formatted.v, label = paste0(v),
                                                    value = "white", allowTransparent = T))
      })
    })
  })
  observe({req(input$upsetgroups)
    upset.groups <- input$upsetgroups
    output$col.title.upSet <- renderUI({ 
      p("Colours:", style = "font-weight: bold; color: #8b0000; background-color: lightgrey;")
    })
    output$upSet.colour.picker <- renderUI({
      lapply(upset.groups, function(u) {
        formatted.u <- paste0("u", u) %>% #add u in front of upset plot colour ids to differentiate from venn colour ids. 
          gsub("\\W", "_", ., perl = T)
        column(width = 4, colourpicker::colourInput(inputId = formatted.u,
                                                    label = paste0(u), value = "grey21", 
                                                    allowTransparent = T))
      })
    })
    
  })
  #Put all colour input ids in vectors for customizing the venn and upset plots
  venn.colour.ids <- reactive({req(input$venngroups)
    grps <- gsub("\\W", "_", input$venngroups, perl = T) #remove non-alphanumeric symbols (\\W) in input ids 
    sapply(grps, function(v) { 
      input[[v]]
    })
  })
  upset.colour.ids <- reactive({req(input$upsetgroups)
    grps <- gsub("\\W", "_", input$upsetgroups, perl = T)
    sapply(grps, function(u) {
      formatted.u <- paste0("u", u) 
      input[[formatted.u]]
    })
  })
  
  
  #Create venn diagram with binary values and colourinput variables
  venn.pl <- reactive({req(venn.input(), input$gen.venn)
    tryCatch({
      venn.fit <- lapply(venn.input()[["binary.result"]],
                         \(x) which(x==1))
      venn.fit <- eulerr::venn(venn.fit)
      pl <- plot(venn.fit,
                 quantities = list(type = input$count_percent,
                                   fontsize = input$numsize),
                 fills= list(fill = venn.colour.ids()),
                 edges = list(lwd = input$linewidth, lty = input$linetype),
                 labels = list(fontsize = input$labelsize))
      return(pl)
    })
  })
  #Generate Venn diagram and add download handler
  output$vennplot <- renderPlot({ venn.pl() })
  output$dlvennplot <- downloadHandler(
    filename = "Venn.png",
    content = function(file) {
      ggsave(file, plot = venn.pl(),
             device = "png", bg = 'white', 
             width = input$width.val, height = input$height.val, dpi = input$dpi.val,
             units = dl.unit)
    }
  )
  
  
  #UpSet Plot
  upset.pl <- reactive({req(df.venn(), input$gen.upset, upset.input())
    upset.colour.ids <- upset.colour.ids()
    tryCatch({
      upset.fit <- as.data.frame(upset.input()[["binary.result"]], check.names = F)
      sn = names(upset.fit)
      pl <- upset(upset.fit, order.by = c("freq"), 
                  text.scale = input$upset.text.sz, 
                  sets.bar.color = upset.colour.ids,
                  keep.order = T,
                  nsets = length(upset.fit), nintersects = NA, 
                  main.bar.color = 'grey21',
                  point.size = input$upset.pt.sz,
                  line.size = input$upset.line.sz,
                  sets = sn)
      return(pl)
    })
  })
  output$upset.plot <- renderPlot({upset.pl()})
  output$dlupsetplot <- downloadHandler(
    filename = "UpSet.png",
    content = function(file) {
      ggsave(file, plot = print(upset.pl()),
             width = input$width.val, height = input$height.val, dpi = input$dpi.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  
  
  
  
  
  # if (!interactive()) {
  #   session$onSessionEnded(function() {
  #     stopApp()
  #     q("no")
  #   })
  # }
  
}

####RUN PROTEOMICS APP -----
shinyApp(ui, server)
