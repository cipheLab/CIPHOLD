# CIPHoLD
<h2> Install dependencies </h2>
To use this tool you need a windows computer with :
<ul>
  <li>R version 3.5.1 or later </li>
  <li>R studio </li>
  <li>R tools 3.5 (administrator installation and check "Path environment" box during installation)</li>
  <li>Google Chrome browser (to download zip output correctly)</li>
</ul>
<p>Now run your Rstudio and you can copy paste this commande line to install all apcajges requierd by CIPHoLD tools: </p>

```
  install.packages(c("devtools","shiny","shinydashboard","shinyjs","gtools","shinyHeatmaply"))
  install.packages("lme4")
  install.packages("BiocManager")
  BiocManager::install("flowCore", version = "3.8")
  BiocManager::install("ggcyto", version = "3.8")
  library("devtools")
  install_github("cipheLab/CIPHoLD")
```

<h2> Run Tool </h2>

```
  library("CIPHoLD")
  CIPHoLD.run()
```

<h2> Presentation </h2>
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

<h2> Upload Data </h2>
<p>You can upload FCS in first tab and create group. Group creation is just a concatenation of FCS with new parameters used to save the files id (sample)</p>
![alt text](https://github.com/cipheLab/CIPHoLD/blob/master/doc/img/01.png)

<p>You can upload csv table, with this files you don't need to clustering. Its result of clustering parts </p>
  
 <h2 Clustering </h2>
 
 
 <h2> Gated Files and Map creation </h2>
 
 
 <h2> Map new dataset </h2>
 
 
 <h2> Explore Mapping </h2>
 
 
 <h2> Download Annotation </h2>
 
 
 
