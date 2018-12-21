# CIPHoLD
<h1> Install dependencies </h1>
To use this tool you need a windows computer with :
<ul>
  <li>R version 3.5.1 or later </li>
  <li>R studio </li>
  <li>R tools 3.5 (administrateor installation and check "Path environment" box during installation)</li>
  <li>Google Chrome browser (to download zip output correctly)</li>
</ul>
<code>
install.packages(c("devtools","shiny","shinydashboard","shinyjs","gtools","shinyHeatmaply"))
install.packages("lme4")
</code>

<code>source("https://bioconductor.org/biocLite.R")</code>

<code>biocLite("ggcyto")</code>

<code>library("devtools") </code>

<code>install_github("cipheLab/CIPHoLD")</code>

<code>library("CIPHoLD")</code>

<code>CIPHoLD.run()</code>
