using Documenter
using CSplineUPM

makedocs(
    sitename = "CSplineUPM",
    format = Documenter.HTML(),
    modules = [CSplineUPM]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
