# Libraries
library(shiny, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(stringr)
library(DT)
library(gridExtra)
library(ggvis)
library(leaflet, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")
library(shinydashboard, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")
library(stargazer, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")
# Functions
source("readData.R")
source("summ.R")

# Don't use factors unless explicitly needed
options(stringsAsFactors = F)

createLink <- function(val) {
  sprintf(
    '<a href="https://www.ebi.ac.uk/gwas/search?query=%s" target="_blank">%s</a>',
    val,
    val
  )
}


# metadata files
sampleInfo = read.csv("/projects/grs-data/explore-app/data/20130606_sample_info.csv",
                      fileEncoding = "latin1")
sampleInfo = sampleInfo[, c(1, 3, 5)]
colnames(sampleInfo) = c("IND", "POP", "GEN")
supInfo = read.csv("/projects/grs-data/explore-app/data/super_population_info.csv",
                   fileEncoding = "latin1")
supInfo = supInfo[, c(1, 3, 7)]
colnames(supInfo) = c("POP", "SUP", "FULL")
sampleInfo = merge(sampleInfo, supInfo, by = "POP")
points = read.csv("/projects/grs-data/explore-app/data/geo.csv")

# color settings
supNames = c("AFR", "EUR", "AMR", "SAS", "EAS")
c5 = c("#484496", "#F4A500", "#328A4F", "#944116", "#D92414")
names(c5) = supNames
fullNames = c("African",
              "European",
              "Admixed American",
              "Southeast Asian",
              "East Asian")
cF = c("#484496", "#F4A500", "#328A4F", "#944116", "#D92414")
names(cF) = fullNames

# path
# path="/projects/grs-data/explore-app/data/snpSets/"
path = "/projects/grs-data/20180403-explore/"

shinyServer(function(input, output, session) {
  session$allowReconnect(TRUE)
  session$onSessionEnded(stopApp)
  observeEvent(input$legend, {
    isolate({
      proxy <- leafletProxy("mymap")
      proxy %>% clearControls()
      if (input$legend) {
        palF <- colorFactor(cF, levels = names(cF))
        proxy %>% addLegend(
          "bottomright",
          pal = palF,
          values = points$FULL,
          title = "Populations"
        )
      }
    })
  })
  
  
  # SNP set categories
  categories = reactive({
    cats = read.delim(
      "/projects/grs-data/explore-app/data/EFO_categories.txt",
      sep = "\t",
      fileEncoding = "latin1"
    )
    cats$Category = ifelse(nchar(cats$Category) == 0, NA, cats$Category)
    cats$Subcategory = ifelse(nchar(cats$Subcategory) == 0, NA, cats$Subcategory)
    cats$Trait = paste0(cats$Trait, ".txt")
    cats
  })
  
  catsNames = reactive({
    cats = categories()
    sort(unique(cats$Category))
  })
  
  output$catsControls <- renderUI({
    cnames = catsNames()
    checkboxGroupInput("selSNPcat", "", cnames, inline = T)
  })
  
  # pre-computed ANOVA
  snpset = reactive({
    cats = categories()
    cnames = catsNames()
    preCmpAnova = read.csv("/projects/grs-data/explore-app/data/preComputedANOVA2.csv",
                           fileEncoding = "latin1")
    # preCmpAnova=preCmpAnova[which(preCmpAnova$p.value < 1), ]
    preCmpAnova = merge(
      preCmpAnova,
      cats,
      by.x = "Names",
      by.y = "Trait",
      all.x = T
    )
    preCmpAnova = preCmpAnova[order(-preCmpAnova$Fstat),]
    if (!is.null(input$selSNPcat)) {
      preCmpAnova[which(preCmpAnova$Category %in% input$selSNPcat),]
    }
    else{
      preCmpAnova
    }
  })
  
  fileList = reactive({
    snpset = snpset()
    fileList = snpset$Names
    featNames = gsub("-unweighted.txt$", "", fileList)
    sel = which(nchar(featNames) > 1)
    fileList = fileList[sel]
    featNames = featNames[sel]
    fileList = data.frame(fileList, featNames)
    colnames(fileList) = c("fileName", "SNPset")
    fileList
  })
  
  # Pre-computed ANOVA
  output$x1 = DT::renderDataTable(
    DT::datatable({
      snpset = snpset()
      snpset[, 1] = gsub(".txt$", "", snpset[, 1])
      snpset[, 1:4]
    }, escape = FALSE, selection = "single", rownames = F, options = list(lengthMenu = c(5, 30, 50, 100), pageLength = 5)) %>%
      formatStyle(columns = 1:4,
                  color = "black")
  )
  output$SNPtable = DT::renderDataTable(
    DT::datatable({
      counts = read.table(
        paste0(path, gsub(
          ".txt", "-unweighted.snpcount", selSnpSet()
        )),
        sep = "\t",
        header = T,
        row.names = NULL
      )
      counts$Calls = counts$Calls / 2504
      counts$SNP <- createLink(counts$SNP)
      counts
    },
    escape = FALSE, selection = "single", rownames = F, options = list(lengthMenu = c(5, 30, 50, 100), pageLength = 5)) %>%
      formatStyle(columns = c("SNP"),
                  color = "black") %>%
      formatPercentage(columns = c('Calls'), 2) %>%
      formatStyle(
        'Calls',
        background = styleColorBar(counts$Calls * 0.5, 'steelblue'),
        backgroundSize = '100% 88%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  )
  selSnpSet <- reactive({
    if (is.null(input$x1_row_last_clicked)) {
      selected = as.character(snpset()[1, 1])
    } else{
      selected = as.character(snpset()[input$x1_row_last_clicked, 1])
    }
    selected
  })
  
  grs <- reactive({
    grs = readGrs(paste0(path, gsub(
      ".txt", "-unweighted.txt", selSnpSet()
    )), sampleInfo)
  })
  
  anovaStat <- reactive({
    grs = grs()
    # ANOVA
    grsMod <<- lm(GRS ~ SUP, data = grs)
    grsAno <<- anova(grsMod)
    fstat = grsAno$`F value`[1]
    pval = grsAno$`Pr(>F)`[1]
    list(fstat = fstat, pval = pval)
  })
  
  
  medians <- reactive({
    grs = grs()
    # Medians
    grsdt = as.data.table(grs)
    grsdt$GRS = as.numeric(grsdt$GRS)
    grsdt[, SUPMED := median(GRS), by = SUP]
    grsdt[, POPMED := median(GRS), by = POP]
    mediansSUP = as.data.frame(grsdt[which(!duplicated(grsdt$SUP)), 5:7, with = F])
    mediansPOP = as.data.frame(grsdt[which(!duplicated(grsdt$POP)), c(3, 5, 8), with =
                                       F])
    list(mediansSUP = mediansSUP, mediansPOP = mediansPOP)
  })
  
  output$mytable = DT::renderDataTable({
    summData = grs()
    summData = summarySE2(summData,
                          measurevar = "GRS",
                          groupvars = c("FULL"))
  },
  options = list(
    paging = FALSE,
    autoWidth = TRUE,
    searching = FALSE
  ),
  rownames = FALSE,)
  
  output$aovSummary = renderText({
    anovaStat()
    stargazer(
      grsAno,
      grsMod,
      type = "html",
      title = c("Summary of ANOVA", "Summar of Linear Regression")
    )
  })
  output$boxplot <- renderPlot({
    grs = grs()
    title = gsub(".txt$", "", selSnpSet())
    fstat = anovaStat()[[1]]
    pval = anovaStat()[[2]]
    medians = medians()[[1]]
    orderedFull = medians$FULL[order(medians$SUPMED)]
    
    # Order SUP by medians
    grs$FULL = ordered(grs$FULL, levels = orderedFull)
    
    g = ggplot(grs, aes(x = FULL, y = GRS)) +
      geom_point(
        position = position_jitter(width = 0.6, height = 0.02),
        aes(
          col = FULL,
          size = .01,
          fill = FULL
        ),
        colour = "black",
        pch = 21,
        show.legend = F
      ) +
      scale_fill_manual(values = cF, guide = F) +
      geom_boxplot(outlier.shape = NA, alpha = 0.3) +
      theme_bw() +
      theme(panel.grid = element_blank()) + theme(text = element_text(size =
                                                                        15)) +
      labs(
        title = title,
        x = paste(
          "F value=",
          format(fstat, digit = 5),
          ", P-value",
          format(pval, digit = 3)
        ),
        y = "Polygenic Trait Score"
      )
    print(g)
  })
  
  output$stripplot <- renderPlot({
    grsdf = grs()
    title = gsub(".txt$", "", selSnpSet())
    medians = medians()[[2]]
    orderedPop = medians$POP[order(medians$POPMED)]
    medians = medians[order(medians$POPMED),]
    
    grsdf$POP = ordered(grsdf$POP, levels = orderedPop)
    grsdf = grsdf[order(grsdf$POP, grsdf$GRS),]
    grslist = split(grsdf, f = grsdf$POP)
    names(grslist) = orderedPop
    
    grslist = lapply(grslist, function(x) {
      x$x = rep(1, nrow(x))
      x$x = 100 * x$x / sum(x$x)
      x
    })
    
    
    grsdf = do.call(rbind, grslist)
    spacer_size = 10
    grsdf$x = cumsum(grsdf$x)
    ymax = max(grsdf$GRS) * 1.2
    ymin = min(grsdf$GRS) * 0.8
    medians$x = seq(1, nrow(medians) * 100, by = 100)
    medians$xend = seq(100, nrow(medians) * 100, by = 100)
    medians$midx = medians$x + (medians$xend - medians$x + 1) / 2
    if (nrow(medians) %% 2 == 0) {
      medians$col = factor(rep(1:2, nrow(medians) / 2))
    } else{
      medians$col = factor(c(rep(1:2, (
        nrow(medians) - 1
      ) / 2), 1))
    }
    medians$midy = (ymax - ymin) / 2 + ymin
    medians$h = ymax - ymin
    medians$ymin = ymin
    #
    ggplot() +
      geom_tile(data = medians, aes(
        x = midx,
        y = midy,
        fill = col,
        height = h
      )) +
      scale_fill_manual(values = c("lightgrey", "white"),
                        guide = F) +
      geom_point(
        data = grsdf,
        aes(x = x, y = GRS, col = SUP),
        alpha = 0.5,
        size = 1
      ) +
      scale_color_manual(values = c5, guide = F) +
      geom_segment(
        data = medians,
        aes(
          x = x,
          y = POPMED,
          xend = xend,
          yend = POPMED
        ),
        size = 1,
        col = "red"
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      ) + theme(text = element_text(size = 15)) +
      labs(y = "Polygenic Trait Score", x = "Population", title = title) +
      geom_text(
        data = medians,
        aes(
          x = midx,
          y = ymin,
          label = POP,
          color = SUP
        ),
        angle = -45,
        vjust = 1,
        hjust = -0.2,
        size = 4,
        fontface = "bold"
      ) +
      coord_cartesian(ylim = c(ymin - (ymax - ymin) * 0.1, ymax))
    
  })
  output$mymap <- renderLeaflet({
    title = gsub(".txt$", "", selSnpSet())
    points = read.csv("/projects/grs-data/explore-app/data/geo.csv")
    mediansPOP = medians()[[2]][, c(1, 3)]
    points = merge(points, mediansPOP, by = "POP")
    points$POPMEDScaled = points$POPMED / mean(points$POPMED)
    points$description = paste0(
      "<p style=\"color:black\"><b>Sample information: </b>",
      points$description,
      " (",
      points$POP,
      ")",
      ". <b>Population size:</b> ",
      points$nSample,
      ". <b>Median GRS for this population:</b> ",
      round(points$POPMED, 4),
      "</p>"
    )
    # Colors
    pal <- colorFactor(c5, levels = names(c5))
    palF <- colorFactor(cF, levels = names(cF))
    if (input$legend) {
      leaflet(options = leafletOptions(dragging = FALSE, zoomControl = FALSE)) %>%
        addProviderTiles("CartoDB", options = providerTileOptions(noWrap = TRUE)) %>%
        fitBounds(
          lng1 = max(points$lng),
          lat1 = 60,
          lng2 = min(points$lng),
          lat2 = min(points$lat)
        ) %>%
        addCircleMarkers(
          lng = ~ lng,
          lat = ~ lat,
          weight = 1,
          radius = ~ (10 * 4 ^ POPMEDScaled / max(POPMEDScaled)),
          popup = ~ description,
          color = ~ pal(SUP),
          stroke = F,
          fillOpacity = 0.5,
          data = points
        ) %>%
        addLegend(
          "bottomright",
          pal = palF,
          values = points$FULL,
          title = "Populations"
        )
    }
    else{
      leaflet(options = leafletOptions(dragging = FALSE, zoomControl = FALSE)) %>%
        addProviderTiles("CartoDB", options = providerTileOptions(noWrap = TRUE)) %>%
        #addTiles(options = providerTileOptions(noWrap = TRUE)) %>%
        fitBounds(
          lng1 = max(points$lng),
          lat1 = 60,
          lng2 = min(points$lng),
          lat2 = min(points$lat)
        ) %>%
        addCircleMarkers(
          lng = ~ lng,
          lat = ~ lat,
          weight = 1,
          radius = ~ (20 * 2 ^ POPMEDScaled / max(POPMEDScaled)),
          popup = ~ description,
          color = ~ pal(SUP),
          stroke = F,
          fillOpacity = 0.5,
          data = points
        )
    }
    
  })
  output$snpFromSet <- renderText({
    ## Change path before deploy on server
    path = "/projects/grs-data/explore-app/data/grsSnp/EFO_new/"
    allfiles = list.files(path)
    filename = selSnpSet()
    if (filename %in% allfiles) {
      tb = as.data.table(read.table(
        paste0(path, filename),
        sep = "\t",
        row.names = NULL
      ))
      if (nrow(tb) == 0) {
        res = "SNP set file is empty"
      } else{
        nSnp = paste0("Number of SNPs: ", length(tb$V3))
        snpL = paste(tb$V3, tb$V4, sep = '-')
        res = snpL
      }
    } else{
      res = "No SNP set information available"
    }
    res
  })
})
