# CIPHoLD
<h1> Install dependencies </h1>
To use this tool you need a windows computer with :
<ul>
  <li>R version 3.5.1 or later </li>
  <li>R studio </li>
  <li>R tools 3.5 (administrateor installation and check "Path environment" box during installation)</li>
  <li>Google Chrome browser (to download zip output correctly)</li>
</ul>
<p>Now run your Rstudio and you can copy paste this commande line to install all apcajges requierd by CIPHoLD tools: </p>

<code>
  install.packages(c("devtools","shiny","shinydashboard","shinyjs","gtools","shinyHeatmaply"))
  install.packages("lme4")
  install.packages("BiocManager")
  BiocManager::install("flowCore", version = "3.8")
  BiocManager::install("ggcyto", version = "3.8")
  library("devtools")
  install_github("cipheLab/CIPHoLD")
</code>

<h1> Run Tool </h1>

<code>
  library("CIPHoLD")
  CIPHoLD.run()
</code>

CIPHoLD Tools is a rewrite of SCAFFOLD tools (https://github.com/nolanlab/scaffold). 
<ul> Difference with scaffold 
  <li> File manager </li>
  <li> Add other clustering method </li>
  <li> Write annotation in FCS files </li>
  <li> All interface update with dashboard </li>
</ul>
We decide to change the input and output system. 
