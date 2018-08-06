library(DT)
library(ggvis)
library(leaflet, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")
library(shinydashboard, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")
library(shiny, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")
library(shinycssloaders, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")
library(glue, lib.loc = "/projects/grs-data/R/x86_64-redhat-linux-gnu-library/3.4/")


body <- dashboardBody(

    includeHTML('header.html'),
    tags$head(
	  	tags$title("GADGET Explorer"),
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
        includeScript("www/ga.js"),
	tags$style("#mymap { background-color: #D4DADC;}")
    ),
    fluidRow(
      column(width = 12,
             box(title = "Browse by Category", width = NULL, status = "primary", solidHeader = T,
                 uiOutput("catsControls")
             )
      ),
      column(width = 6,
             box(width = NULL, status = "primary",
                 DT::dataTableOutput('x1')
             ),
             box(title = "Summary statistics", width = NULL, status = "warning", solidHeader = T, collapsible = T, collapsed = T,
              DT::dataTableOutput("mytable"),
              DT::dataTableOutput("SNPtable")%>% withSpinner(color="#F39C12"),
              htmlOutput('aovSummary') %>% withSpinner(color="#F39C12"),
              br(),
              # plotOutput("barplot"),
              br()
            
             ),
             box(title = "Variants used in calculation", width = NULL, status = "warning", solidHeader = T, collapsible = T, collapsed = F,
                 p(class = "text-muted",
                   br(),
                   "SNPs and risk alleles:"
                 ),
                 textOutput("snpFromSet")
             )
      ),
      column(width = 6,
             box(title = "Distribution of PTS Worldwide", width = NULL, status = "primary", collapsible = T, collapsed = F,
                 leafletOutput("mymap", height=300) %>% withSpinner(color="#3C8DBC") ,
    checkboxInput("legend", "Show legend", TRUE)
             ),
             box(title= "PTS Distribution by Population", width = NULL, status = "warning", 
                 plotOutput("stripplot") %>% withSpinner(color="#F39C12") #,
                 #textOutput("test")
             ),
             box(title = "Human population abbreviation key", width = NULL, status = "warning", solidHeader = T, collapsible = T, collapsed = T,
                 p(class = "text-muted",
                   br(),
                   # Labels
                   includeMarkdown("labels.md")
                 )
             ),
             box(title="PTS Distribution by Continental Group", width = NULL, status = "warning",  collapsible = T, collapsed = F,
                 plotOutput("boxplot") %>% withSpinner(color="#F39C12")
             )
             
             
             
      )
    ),
    fluidRow(column(width = 12,
             box(title =NULL, width = NULL, status = "danger", solidHeader = F,
                 p("GADGET is intended as a tool for researchers to explore population-specific distributions of genetic variants that have been associated with a wide variety of human traits.  
      Users of this site should treat the results with caution, as the interpretation PTS across populations can be complicated by a number of factors.  We provide more detail on these issues on the
       Learn page and in our manuscript describing the server.")
             )
             )
    )

)

dashboardPage(
	title = "GADGET Explore",
  dashboardHeader(disable = TRUE),
  dashboardSidebar(disable = TRUE),
  body
)
