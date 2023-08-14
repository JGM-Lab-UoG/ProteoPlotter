#Esther Olabisi-Adeniyi
#Bioinformatics Master's Research Project - BINF6999
#August 8, 2023

require(ggplot2)
require(ggpubr)
require(gplots)
require(ggrepel)
require(tidyverse)
require(reshape2)
require(shinyWidgets)
require(scales)



####1. PROCESS DATA MATRIX ----
require(RColorBrewer)
require(ggnewscale)

#Read protein.txt matrix processed in Perseus for a volcano plot. Remove all comment rows. 
#Change column names to shorter generic names that are easier to work with. 
#prot.dataset: Imported from Perseus following quality control, two samples t-test, and GO enrichment
#left.range: Intensity columns belonging to the left group of the volcano plot / two samples test in Perseus
#right.range: The same as left.range but for the right group

process.df <- function(prot.dataset, left.range, right.range) {
  comment_char <- grep("^#", prot.dataset[,1])
  prot.dataset <- prot.dataset[-(comment_char),]
  col.names <- colnames(prot.dataset)
  intensity.cols <- grep("intensity.+", col.names, ignore.case = T, value = T)
  major.prot.col <- grep("majority", col.names, ignore.case = T, value = T)
  grp.name.col <- grep("significant", col.names, value = T) 
  
  #select relevant columns only
  prot.dataset <- prot.dataset %>% 
    setNames(sub(".*p.value.*", "minus.log10.pval", names(.))) %>% #-log10pval
    setNames(sub(".*Difference.*", "Difference", names(.))) %>% #Difference
    setNames(sub(".*Significant.*", "Significant", names(.))) %>%
    setNames(sub(".*biological.*", "Biological Process", ignore.case = T, names(.))) %>%
    setNames(sub(".*cellular.*", "Cellular Component", ignore.case = T, names(.))) %>%
    setNames(sub(".*molecular.*", "Molecular Function", ignore.case = T, names(.))) %>%
    setNames(sub(".*keywords.*", "Keywords", ignore.case = T, names(.))) %>%
    select(all_of(intensity.cols), 
           all_of(major.prot.col), 
           all_of(grp.name.col),
           minus.log10.pval, Difference, Significant, `Biological Process`, `Cellular Component`, `Molecular Function`, Keywords)
  
  
  #Extract group labels from the "significant" column for use in legend
  sig.df <- prot.dataset[prot.dataset$Significant == "+",]  
  grp.names <- unique(sig.df[, grp.name.col]) #extract label 
  grp.names <- str_split_1(grp.names, "_") #separate labels
  
  #Add the group labels based on positive or negative x axis then calculate mean across groups
  prot.dataset <- prot.dataset %>%
    mutate(Expression = case_when(Difference >= 0 & Significant == "+" ~ paste0("Higher abundance in ", grp.names[1]),
                                  Difference <= 0 & Significant == "+" ~ paste0("Higher abundance in ", grp.names[2]),
                                  TRUE ~ "Not significant")) %>%
    mutate(Expression = factor(Expression, levels = c(
      paste0("Higher abundance in ", grp.names[2]),
      paste0("Higher abundance in ", grp.names[1]),
      "Not significant"))) %>%  #reorder
    mutate_at(c('Difference', 'minus.log10.pval', intensity.cols), as.numeric) %>%
    rowwise() %>%
    mutate(Left.group = mean(
      c_across(all_of(left.range)))) %>%
    mutate(Right.group = mean(
      c_across(all_of(right.range)))) %>%
    remove_rownames() %>%
    column_to_rownames(var = major.prot.col)
  
  
  return(prot.dataset)
}





####2. 1D HEATMAPS ----
onedheatmap <- function(oned.df, plot.title = "") {
  
  oned.df <- oned.df %>%
    mutate_at(c("Score"), as.numeric)
  
  #some annotations have the same name but belong to a different database
  #add additional column combining those entries
  oned.df$unique.annotation <- paste(oned.df$Name, " (", oned.df$Type , ")", sep="")
  
  #unique annotations and celltypes
  annotations <- sort(unique(oned.df$unique.annotation))
  annotations <- annotations[!grepl("^\\+", annotations)]
  celltypes <- sort(unique(oned.df$Column))
  
  #generate 1D score matrix
  M.score <- matrix(ncol=length(celltypes), nrow=length(annotations))
  
  for(i in 1:length(celltypes)){
    for(j in 1:length(annotations)){
      score.value <- oned.df$Score[oned.df$Column==celltypes[i] & oned.df$unique.annotation==annotations[j]]
      if(length(score.value)==0){score.value <- NA}
      M.score[j,i] <- score.value
    }
  }
  
  m <- M.score
  rownames(m) <- annotations
  colnames(m) <- celltypes
  
  #Sequence of samples ('cell types' as currently above) for columns of heat map
  m <- m[,c(1:8)]
  m[is.na(m)] <- 0
  
  #collapse matrix for use in ggplot
  m <- melt(m)
  colnames(m) <- c("Annotations", "T-test differences", "value")
  
  return(
    ggplot(data = m, aes(x = `T-test differences`, 
                         y = reorder(Annotations, value), 
                         fill = value))+
      geom_tile(colour = "grey", linewidth = 1)+
      scale_fill_gradientn(colors = c(low = "blue", mid = "white", high = "red"),
                           na.value = "grey")+
      scale_y_discrete(position = "right")+
      guides(fill = guide_colourbar(title = "Colour Key"))+
      theme(axis.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 90),
            axis.text = element_text(colour = "black"),
            legend.position = "left")+
      labs(title = plot.title, 
           y = "Annotations")
  )
  
}





####3. VOLCANO Plot ----
#The function requires a dataframe that has been processed by process.df() (step 1). It also required the curves.df exported from Perseus to create the lines on the volcano plot. The other arguments are optional.
#go.terms = display GO terms for significant or non-significant proteins
#plot.title = custom plot title top left
#group.cols = colours for the default plot in the following order: higher-left, higher-right colours. 
#s0 = s0 value to be shown on the plot. Default = 0.1
#fdr = fdr value to be shown on the plot. Default = 0.05
#fdr.lines = 'yes' to show lines and 'no' to remove lines
#palette.col = based on hcl.color() palettes. "Viridis" is the default one

volcano_plot <- function(df, curves.df, 
                         go.terms = "", which.go = "Keywords", plot.title = "", 
                         s0 = 1, fdr = 0.05, fdr.lines = "yes", 
                         palette.col = "Viridis", group.cols = c("#00B1B2", "#FF6666"),
                         vary.sizes = "no",
                         select.pts = c()) {
  
  #Subset significant and non-significant(ns) proteins into their respective dfs. 
  #Extract the number of unique GO terms in both sig and ns dfs. The numbers will be used for color-coding the GO terms downstream. 
  #Subset the selected proteins so that we can use their protein ids(rownames) as point labels on the plot
  #Also, extract the x-axis label (in "Group1_Group2" format) from "significant" column
  sig.df <- df[df$Significant == "+",]
  ns.df <- df[df$Significant == "",] 
  num.sig.go <- length(unique(sig.df[, which.go])) 
  num.ns.go <- length(unique(ns.df[, which.go])) 
  select.pts <- df[select.pts, ]
  grp.name.col <- grep("significant", colnames(sig.df)) #column index 
  xlab <- unique(sig.df[, grp.name.col]) 
  xlab <- sub("_", " - ", xlab) #for x axis label 
  grp.names <- str_split_1(xlab, " - ") #for legend label
  
  
  #Custom settings for ggplots
  a = 0.2 #alpha
  s = 21  #shape
  ns.sz = 2 #non-sig shape size  
  ns.col = "grey" #non-sig colour
  s.sz = 4  #sig shape size
  
  #Base volcano plot to plot points, fdr curves, axes boundaries, and other custom arguments.
  base_vplot <- ggplot(data = ns.df, 
                       aes(x = Difference,
                           y = minus.log10.pval))+
    geom_point(aes(fill = ns.col), size = ns.sz,
               alpha = a, shape = s)+
    labs(title = plot.title,
         caption = paste("s0 =", s0, "  ", "FDR =", fdr),
         x = paste0("Difference (", xlab, ")"), 
         y = expression(-log[10]("p-value")))+
    coord_cartesian(xlim = c(min(df$Difference),
                             max(df$Difference)),
                    ylim = c(min(df$minus.log10.pval),
                             max(df$minus.log10.pval)),
                    expand = T)+
    theme_minimal()+
    theme(plot.title = element_text(face="bold", size = 20),
          plot.caption = element_text(
            colour = "darkviolet", size = 11))
  
  
  #Fdr curves layer 
  curve.plot <- geom_line(data = curves.df, aes(x, y), linetype=2)
  
  #Layer for labelling points 
  label.pt.plot <- geom_label_repel(data = select.pts, 
                                    aes(label = rownames(select.pts)),
                                    nudge_y = 0.4,  fill = "grey")
  
  ###Adding layers to the base plot
  #Add GO terms for sig proteins
  sig.go <- base_vplot +
    scale_fill_identity()+  #hides grey legend that comes from the base plot
    new_scale_fill()+
    geom_point(data = sig.df,
               aes(fill = !!sym(which.go)), 
               size = s.sz,
               shape = s)+
    scale_fill_manual(values = hcl.colors(
      num.sig.go, 
      palette.col), 
      na.value = alpha(ns.col, a)) + 
    label.pt.plot
  
  #Add GO terms for non-sig proteins. Grey out sig proteins. 
  nonsig.go <- base_vplot +
    scale_fill_identity()+
    new_scale_fill()+
    geom_point(aes(fill = !!sym(which.go)),
               size = ns.sz, 
               shape = s)+
    scale_fill_manual(values = hcl.colors(
      num.ns.go,
      palette.col),
      na.value = alpha(ns.col, a)) + 
    geom_point(data = sig.df, fill = "grey",
               size = s.sz, alpha = a,
               shape = s) + 
    label.pt.plot
  
  ###Add layer showing average protein abundance across samples 
  vary.size.pl <- base_vplot +
    scale_fill_identity()+
    new_scale_fill()+
    geom_point(data = sig.df[sig.df$Difference > 0, ],
               aes(size = Right.group),
               fill = alpha("#FF6666", a),
               shape = s)+
    scale_size_binned(range = c(1,8),
                      name = paste0("Average Intensity ", grp.names[1]),
                      n.breaks = 4)+
    new_scale(new_aes = "size")+
    geom_point(data = sig.df[sig.df$Difference < 0, ],
               aes(size = Left.group), 
               fill = alpha("#00B1B2", a),
               shape = s)+
    scale_size_binned(range = c(1,8),
                      name = paste0("Average Intensity ", grp.names[2]),
                      n.breaks = 4) + 
    label.pt.plot
  
  #Options: GO term plot for sig proteins, with fdr lines and without fdr lines
  if (go.terms == "significant" && fdr.lines == "yes") {
    return (sig.go + curve.plot)
    
  } else if (go.terms == "significant" && fdr.lines == "no") {
    return (sig.go)
    
    #Options: GO terms for non-significant proteins with fdr lines, without fdr lines, and default plot with color-coding of high and low expression proteins. 
  } else if (go.terms == "non-significant" && fdr.lines == "yes") {
    return (nonsig.go + curve.plot)
    
  } else if (go.terms == "non-significant" && fdr.lines == "no") {
    return (nonsig.go)
    
    #Sizes
  } else if(vary.sizes == "yes" && fdr.lines == "yes") {
    return (vary.size.pl + curve.plot)
    
  } else if (vary.sizes == "yes" && fdr.lines == "no") {
    return (vary.size.pl)
    
  } else if (fdr.lines == "yes") {
    return(
      base_vplot +
        scale_fill_identity(name = NULL,
                            labels = "Not significant",
                            guide = "legend")+
        new_scale_fill()+
        geom_point(data = sig.df, aes(fill = Expression),
                   size = s.sz, shape = s) +
        scale_fill_manual(values = group.cols)+
        guides(fill = guide_legend(order = 1))+
        label.pt.plot + curve.plot + 
        theme(legend.spacing = unit(-0.5, "cm"))
    )
    
  } else if (fdr.lines == "no") {
    return(
      base_vplot +
        scale_fill_identity(name = NULL,
                            labels = "Not significant",
                            guide = "legend")+
        new_scale_fill()+
        geom_point(data = sig.df, aes(fill = Expression),
                   size = s.sz, shape = s) +
        scale_fill_manual(values = group.cols)+
        guides(fill = guide_legend(order = 1))+
        label.pt.plot +
        theme(legend.spacing = unit(-0.5, "cm"))
    )
  }
  
}





####4. PCA plot ---- 
library(ggfortify)
library(ggrepel)

#This function also accepts the processed dataframe from step 1 and creates a PCA of all intensity columns.
#Serial ID are removed from group names so that eclipses are drawn by group
pca_plot <- function(df) {
  counts <- df[, grep("intensity", colnames(df), value = T)]
  counts <- na.omit(counts)
  counts <- as.data.frame(t(counts)) 
  pca <- prcomp(counts) 
  counts$Group.names <- sub("\\.\\d+$", "", rownames(counts)) #remove serial numbers to extract group names
  
  autoplot(pca, data = counts, scale = 0, colour = "Group.names", frame = T, frame.type = "norm")+
    guides(colour=guide_legend("Sample Type"), fill = "none")+
    coord_fixed()+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          panel.border = element_blank())
}




####5. S-curve plot ----
#The function accepts the same dataframe produced after step 1 but only produces S-curves for the group pair (e.g WT vs MT) in the volcano plot
#The function create one s-curve plot for each group and puts them on the same image with ggarrange

scurve <- function(df) {
  ranked.df <- df %>%
    mutate(Left.rank = rank(desc(Left.group)))%>%
    mutate(Right.rank = rank(desc(Right.group)))
  
  sig.df <- df[df$Significant == "+", ]
  grp.names <- unique(
    sig.df[, grep("significant", colnames(sig.df))] ) 
  grp.names <- str_split_1(grp.names, "_")  
  
  left.gr.pl <- ggplot(data = ranked.df,
                       aes(x = Left.rank, y = Left.group))+
    geom_point(size = 2, colour = "royalblue")+
    labs(x = paste(grp.names[2], "Proteins ranked by iBAQ Intensity"),
         y = bquote("iBAQ Intensity for " ~ .(grp.names[2]) * " " * (Log[10])))+
    theme_minimal()
  
  right.gr.pl <- ggplot(data = ranked.df,
                        aes(x = Right.rank, y = Right.group))+
    geom_point(size = 2, colour = "royalblue")+
    labs(x = paste0(grp.names[1], "proteins ranked by iBAQ Intensity"),
         y = bquote("iBAQ Intensity for " ~ .(grp.names[1]) * " " * (Log[10])))+
    theme_minimal()

  
  return(
    #plot with a tiny column in between to serve as a larger gap
    ggarrange(left.gr.pl, NULL, right.gr.pl, 
              nrow = 1, 
              widths = c(1, 0.05, 1))
  )
}  





####6. SHINY ----
require(shiny)
require(shinyjs)
require(colourpicker)

#It accepts the protein matrix from Perseus and processes it using step 1 for the remaining program. It also accepts the 1D matrix from Perseus to create 1D heatmaps and to filter out less significant GO terms on the volcano plot. Thirdly, it accepts the curve matrix from Perseus for the lines on the volcano plot 
#A Shiny app's function has two sections as seen below: the UI and the Server

ProteomicsApp <- shinyApp(
  #Shiny UI ---- 
  ui <- fluidPage(
    titlePanel("Visualization for Proteomics"),
    sidebarLayout(
      
      sidebarPanel(width = 4,
                   
                   fluidRow(
                     column(width = 4,
                            fileInput(inputId = "proteinfile", label = "Input Protein Groups file") ),
                     column(width = 4, 
                            fileInput("curvesfile", label = "Input file for volcano curves") ),
                     column(width = 4,
                            fileInput("onedfile", label = "1D annotation data input") ),
                     column(width = 6, 
                            actionButton("gen.pca", label = "Generate PCA plot")),
                     column(width = 6, 
                            actionButton("gen.scurve", label = "Generate S-curve plot")),
                     column(width = 12, br()),
                     uiOutput("enter.range"),
                     column(width = 12, selectInput("left.gr", 
                                                    label = "Column range for Left group",
                                                    choices = c(),
                                                    multiple = T) ),
                     column(width = 12, selectInput("right.gr", 
                                                    label = "Column range for Right group",
                                                    choices = c(),
                                                    multiple = T) ),
                     column(width = 12, br()),
                     column(width = 12, br()),
                   ),
                   
                   useShinyjs(),
                   fluidRow(
                     div(id = "reset.grp",
                         h4("ggplot Options"), 
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
                         
                         column(width = 12, actionButton("reset", label = "Reset Volcano Plot")),
                         column(width = 12, br()),
                         column(width = 6, colourpicker::colourInput("left.col",
                                                                     label = "Left group colour",
                                                                     value = "#00B1B2") ),
                         column(width = 6, colourpicker::colourInput("right.col",
                                                                     label = "Right group colour",
                                                                     value = "#FF6666") ),
                         column(width = 12, br()),
                         column(width = 6, numericInput("s.knot",
                                                        label = "s0 value",
                                                        value = 1) ),
                         column(width = 6, numericInput("fdr.val",
                                                        label = "FDR value",
                                                        value = 0.05) ),
                         column(width = 12, 
                                helpText("Note: Changing the s0 and FDR values does not change the plot. The caption only indicates the parameters used for creating the plot.") ),
                         column(width = 12, br()),
                         column(width = 6, radioButtons("fdr.lines",
                                                        label = "Display curves",
                                                        choices = list("yes", "no"),
                                                        selected = "yes", inline = T) ),
                         column(width = 6, radioButtons("pt.sizes",
                                                        label = "Average protein intensities",
                                                        choices = list("yes", "no"),
                                                        selected = "no", inline = T) ),
                         column(width = 12, radioButtons("go.term",
                                                         label = "Gene Ontology terms",
                                                         choices = list("significant", "non-significant", "reset"),
                                                         selected = "reset", inline = T) ),
                         column(width = 7, selectInput("go.types",
                                                       label = "Gene Ontology annotation category for volcano plot",
                                                       choices = list("Keywords", "Cellular Component", "Molecular Function", "Biological Process"),
                                                       selected = "Keywords", 
                                                       width = '100%') ),
                         column(width = 5, selectInput("palette",
                                                       label = "GO terms color palette",
                                                       choices = list("Viridis", "Plasma", "Inferno", "Rocket", "Lajolla", "Turku", "Hawaii", "Batlow", "Spectral", "Blue-Red", "Green-Orange", "RdYlBu", "Zissou 1", "Roma"),
                                                       selected = "Viridis") )
                     )
                   )           
      ),
      
    
        
      mainPanel(
        textOutput("text"),
        uiOutput("df.preview"),
        dataTableOutput("dataframe"),
        br(),
        
        h4("Volcano Plot",  align = "center"),
        ##Volcano plot output with brush selection of points
        plotOutput("volcanoplot",
                   brush = "vplot.brush"),
        textInput("vplot.title",
                  label = "Volcano Plot Title (Optional)",
                  value = NULL,
                  placeholder = "Enter text..."),
        fluidRow(
          ##initialize drop down list of protein labels
          column(width = 6, selectizeInput("protein.labs",
                                           label = "Manually select protein IDs:",
                                           choices = list(),
                                           multiple = T)),
          column(width = 6, verbatimTextOutput("brushed.ids")),
          column(width = 6, downloadButton("dlvplot", "Download Volcano plot"))
        ),
        br(),
        br(),
        
        ##S-curve plot output
        h4("S-curve Plots",  align = "center"),
        plotOutput("scurveplot"),
        downloadButton("dlscurve", "Download S-curves"),
        br(),
        br(),
        
        #PCA plot output
        h4("PCA Plot", align = "center", ),
        plotOutput("pcaplot"),
        downloadButton("dlpca", "Download PCA"),
        
        br(),
        br(),
        ##1D plots output
        h4("1D Annotation Heatmaps", align = "center", ),
        
        plotOutput("onedhm1.bio.proc"),
        downloadButton("dlhm1", "Download Heatmap 1"),
        br(),
        plotOutput("onedhm2.cell.comp"),
        downloadButton("dlhm2", "Download Heatmap 2"),
        br(),
        plotOutput("onedhm3.mol.func"),
        downloadButton("dlhm3", "Download Heatmap 3"),
        br(),
        plotOutput("onedhm4.keywords"),
        downloadButton("dlhm4", "Download Heatmap 4"),
      )
    )
  ),
  
  
  
  
  #Shiny Server ----
  server <- function(input, output, session) {
    
    #Curves file
    df.curves <- eventReactive(input$curvesfile, {
      read.table(input$curvesfile$datapath, header = T)
    })
    
    #Protein input file
    df.prot <- eventReactive(input$proteinfile, {
      read.delim(input$proteinfile$datapath, na.strings = c("NA", "NaN"))
    })
    
    #1D annotation file
    df.1d <- eventReactive(input$onedfile, {
      df <- read.delim("../1d annot.txt", stringsAsFactors = F, header=T)
      df <- df[-c(1),]
    })
    
    #Update column selection on the UI
    observe({
      req(df.prot())
      updateSelectInput(session,
                        inputId = "left.gr",
                        choices = colnames(df.prot()))
      updateSelectInput(session,
                        inputId = "right.gr",
                        choices = colnames(df.prot()))
      output$enter.range <- renderUI({
        helpText("To calculate average protein intensities across groups, select columns for each group:")
      })
    })
    
    #Process the protein dataframe and calculate means across group columns
    df.prot2 <- reactive({
      req(df.prot(), input$left.gr, input$right.gr)
      df <- process.df(prot.dataset = df.prot(),
                       left.range = input$left.gr,
                       right.range = input$right.gr )
      
      if (is.null(input$onedfile)) {
        df
        
      } else {
        # Retain only those GO terms found in the 1D Annotation file
        df$Keywords[-which(df$Keywords %in% df.1d()$Name)] <- NA
        df$`Cellular Component`[-which(df$`Cellular Component` %in% df.1d()$Name)] <- NA
        df$`Molecular Function`[-which(df$`Molecular Function` %in% df.1d()$Name)] <- NA
        df$`Biological Process`[-which(df$`Biological Process` %in% df.1d()$Name)] <- NA
        
        df
      }    
    })

    #Dataframe preview on UI
    observe({
      req(df.prot2())
      output$df.preview <- renderUI({
        helpText(h4("Dataframe preview")) })
      output$dataframe <- renderDataTable( df.prot2() %>% head(3) )
      
      #fill dropdown list with protein names (labels) from which to choose
      updateSelectizeInput(session,
                           inputId = "protein.labs",
                           choices = rownames(df.prot2()), 
                           server = T )
    })
    
    #activate brush selection    
    brushed.pts <- reactive({
      brushedPoints(df.prot2(), input$vplot.brush)
    })
    #display selected ids as text on the UI
    output$brushed.ids <- renderPrint({ 
      cat("Selected proteins: ", rownames(brushed.pts()) ) 
    }) 
    

        
    ####create volcano plot using UI input variables  
    vplot <- reactive({
      vplot <- volcano_plot(df = df.prot2(),
                            curves.df = df.curves(),
                            go.terms =  input$go.term,
                            palette.col = input$palette,
                            plot.title = input$vplot.title,
                            which.go = input$go.types,
                            s0 = input$s.knot,
                            fdr = input$fdr.val,
                            group.cols = c(input$left.col, input$right.col),
                            fdr.lines = input$fdr.lines,
                            vary.sizes = input$pt.sizes,
                            #show protein labels either by quick or fixed selection 
                            select.pts = c(rownames(brushed.pts()), 
                                           input$protein.labs)
      )
      
      if (is.na(input$xmin) && is.na(input$xmax) && is.na(input$ymin) && is.na(input$ymax)) {
        vplot
        
      } else {
        #Axes limits adjustment
        vplot+
          coord_cartesian(xlim = c(input$xmin, input$xmax),
                          ylim = c(input$ymin, input$ymax),
                          expand = T)
      }
    })
    #send vplot to the UI
    output$volcanoplot <- renderPlot({
      vplot()
    })
    #Reset to default plot
    observeEvent(input$reset, {
      shinyjs::reset("reset.grp")
    })
    output$dlvplot <- downloadHandler( 
      filename = "Volcano_plot.png",
      content = function(file) {
        ggsave(file, plot = vplot(), width = 13, height = 8, dpi = 600, bg = 'white')
      }
    )
    
    
    
    ####PCA plot
    pca.pl <- reactive({req(input$gen.pca) 
      pca_plot( df.prot2() )
    })
    output$pcaplot <- renderPlot({
      pca.pl()
    })
    output$dlpca <- downloadHandler( 
      filename = "PCA.png",
      content = function(file) {
        ggsave(file, plot = pca.pl(), width = 13, height = 8, dpi = 600, bg = 'white')
      }
    )
    
    
    
    
    ####S-curve plot 
    scurve.pl <- reactive({
      req(input$gen.scurve)
      scurve(df.prot2()) 
    })
    output$scurveplot <- renderPlot({
      scurve.pl() 
    })
    #download
    output$dlscurve <- downloadHandler( 
      filename = function() { paste("scurve.png") },
      content = function(file) {
        ggsave(file, plot = scurve.pl(), width = 13, height = 7, dpi = 600, bg = 'white')
      }
    )
    
    
    
    
    ####Heatmap plots
    #create heatmap list
    plotlist.1d <- reactive({
      unique.annot.types <- df.1d() %>%
        group_by(Type) %>%
        group_split()
      lapply(unique.annot.types, 
             function(x) onedheatmap(x, 
                                     plot.title = paste(unique( x$Type ))
             ))
    })
    #plot heatmaps individually
    bio.proc <- reactive(plotlist.1d()[[1]])
    cell.comp <- reactive(plotlist.1d()[[2]])
    mol.func <- reactive(plotlist.1d()[[3]])
    keywords <- reactive(plotlist.1d()[[4]])
    output$onedhm1.bio.proc <- renderPlot({ bio.proc() }) 
    output$onedhm2.cell.comp <- renderPlot({ cell.comp() }) 
    output$onedhm3.mol.func <- renderPlot({ mol.func() }) 
    output$onedhm4.keywords <- renderPlot({ keywords() })
    #downloads
    output$dlhm1 <- downloadHandler( 
      filename = "Bio_process.png",
      content = function(file) {
        ggsave(file, plot = bio.proc(), width = 13, height = 8, dpi = 600, bg = 'white')
      }
    )
    output$dlhm2 <- downloadHandler( 
      filename = "Cell_component",
      content = function(file) {
        ggsave(file, plot = cell.comp(), width = 13, height = 8, dpi = 600, bg = 'white')
      }
    )
    output$dlhm3 <- downloadHandler( 
      filename = "Mol_function.png",
      content = function(file) {
        ggsave(file, plot = mol.func(), width = 13, height = 8, dpi = 600, bg = 'white')
      }
    )
    output$dlhm4 <- downloadHandler( 
      filename = "Keywords.png",
      content = function(file) {
        ggsave(file, plot = keywords(), width = 13, height = 8, dpi = 600, bg = 'white')
      }
    )
    
    
  }
  
)

####RUN PROTEOMICS APP -----
runApp(ProteomicsApp)
