library("shinydashboard")
library("shinyjs")
library("shinyHeatmaply")
source("functions.R")

dashboardPage(
  dashboardHeader(
  	title = "CIPHOLD"
  ),
  dashboardSidebar(
  	shinyjs::useShinyjs(),
  	sidebarMenu(id="tab",
  		menuItem("Multifile Groups", tabName="multiGroups", icon = icon("cog", lib = "glyphicon")),
      menuItem("Run Clustering", tabName = "clustering", icon = icon("dashboard")),
      menuItem("SCAFFOLD Analysis", tabName = "analysis", icon = icon("th-list")),
      menuItem("Map Exploration", tabName = "mapping", icon = icon("braille")),
      menuItem("Map Dataset", tabName = "mapdata", icon = icon("signal")),
      menuItem("Export Data", tabName = "export", icon=icon("download"))
    )
  ),
  dashboardBody(
  	shinyjs::useShinyjs(),
  	tags$head(tags$script(src = "jquery-ui.min.js")),
		singleton(tags$head(
			tags$link(rel = 'stylesheet', type = 'text/css', href = 'custom.css')
		)),
		tags$head(tags$script(src = "d3.min.js")),
		tags$head(tags$script(src = "graph.js")),
		tags$head(tags$script(src = "rect_select.js")),
		singleton(tags$head(
			tags$link(rel = 'stylesheet', type = 'text/css', href = 'rect_select.css')
		)),
		singleton(tags$head(
			tags$link(rel = 'stylesheet', type = 'text/css', href = 'graph.css')
		)),
		singleton(tags$head(
			tags$link(rel = 'stylesheet', type = 'text/css', href = 'style.css')
		)),
		fluidRow(
	  	column(9,
		  	tabItems(
		  		tabItem(tabName="multiGroups",
		  			uiOutput("boxDelBut"),
		  			uiOutput("boxMultiGroup")
		  		),
		      tabItem(tabName = "clustering",
		      	uiOutput("boxPreProcess"),
		      	uiOutput("boxRunClustering"),
		      	uiOutput("boxOverClustering"),
		      	uiOutput("boxHeatmap")
					),
		      tabItem(tabName = "analysis",
		        column(9,
		          uiOutput("boxPreProcessGated"),
	            uiOutput("boxMapMarkers"),
		          uiOutput("boxParameters")
		        ),
		      	column(3,
			      	uiOutput("boxGated"),
			      	uiOutput("boxGatedTable")
			      )
		      ),
		      tabItem(tabName = "mapping",
		      	column(12,
		      		uiOutput("boxMap"),
		      		uiOutput("boxTable")
		      	)
		      ),
		      tabItem(tabName = "mapdata",
		      	column(9,
		          uiOutput("boxMapDataMarkers"),
		          uiOutput("boxMapDataParams")
		         ),
		        column(3,
		          uiOutput("boxScaffoldMapData"),
		          uiOutput("boxScaffoldMapView")
		        )
		      ),
		      tabItem(tabName = "export",
		      	column(9,
		      		uiOutput("boxExportNames"),
		      		uiOutput("boxSelectEnrichment")
		      	),
		      	column(3,
		      		uiOutput("boxScaffoldMapExport"),
		      		uiOutput("boxScaffoldMapPops")
		      	)
		      )
		    )
		  ),
	  	column(3,
	  		uiOutput("boxInput"),
	  		uiOutput("boxInputTable"),
	  		uiOutput("boxOutput"),
	  		uiOutput("warning"),
	  		uiOutput("boxScaffoldMap"),
	  		uiOutput("boxScaffoldMapUI"),
	  		uiOutput("boxScaffoldMapDownload")
	  	)
	  )
  )
)