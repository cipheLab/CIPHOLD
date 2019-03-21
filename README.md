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

<h3> Run Tool </h3>

<code>
  library("CIPHoLD")
  CIPHoLD.run()
</code>

<h3> Presentation </h3>
CIPHoLD Tools is a rewrite of SCAFFOLD tools (https://github.com/nolanlab/scaffold). 
<ul> <h4>Difference with scaffold</h4> 
  <li> File manager </li>
  <li> Upload CSV file (clustering from other tools) </li>
  <li> Add other clustering method </li>
  <li> Write annotation in FCS files </li>
  <li> All interface update with dashboard </li>
  <li> Different output format and information </li>
</ul>
We decide to change the input and output system. 

<h3> Upload Data </h3>
<p>You can upload FCS in first tab and create group. Group creation is just a concatenation of FCS with new parameters used to save the files id (sample)</p>
<p>You can upload csv table, with this files you don't need to clustering. Its result of clustering parts </p>
  
 <h3> Clustering </h3>
 
 
 <h3> Gated Files and Map creation </h3>
 
 
 <h3> Map new dataset </h3>
 
 
 <h3> Explore Mapping </h3>
 
 
 <h3> Download Annotation </h3>
 
 
 
