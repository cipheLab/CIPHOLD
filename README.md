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

install.packages(c("devtools","shiny","shinydashboard","shinyjs","gtools","shinyHeatmaply"))<br/>
install.packages("lme4")<br/>
install.packages("BiocManager")<br/>
BiocManager::install("flowCore", version = "3.8")<br/><br/>
BiocManager::install("ggcyto", version = "3.8")<br/>
install_github("cipheLab/CIPHoLD")<br/>
library("CIPHoLD")<br/>
CIPHoLD.run()<br/>


<h1> Run Tool </h1>
