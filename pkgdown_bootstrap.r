library(pkgdown)

# this doesn't work. I have no idea why.
#setwd("./BayesianTools")

# sink("_pkgdown.yml")
# template_articles("./BayesianTools")
# template_navbar("./BayesianTools")
# template_reference("./BayesianTools")
# 
# sink()
# 
# dir = getwd()
# if (!dir.exists("./docs")) dir.create("docs")
# build_home(pkg = file.path(dir, "BayesianTools"), path = file.path(dir, "docs"))
#setwd("..")

setwd("BayesianTools")

sink("../_pkgdown.yml")
template_navbar(".")
template_reference(".")
template_articles(".")
sink()

if (!dir.exists("../docs")) dir.create("docs")
build_site(path = "../docs")

setwd("..")
