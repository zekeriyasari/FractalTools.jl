using Documenter, FractalTools

makedocs(;
    modules=[FractalTools],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/zekeriyasari/FractalTools.jl/blob/{commit}{path}#L{line}",
    sitename="FractalTools.jl",
    authors="Zekeriya Sarı, Gizem Kalender",
    assets=String[],
)

deploydocs(;
    repo="github.com/zekeriyasari/FractalTools.jl",
)
