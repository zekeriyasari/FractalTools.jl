using Documenter, FractalTools

DocMeta.setdocmeta!(FractalTools, :DocTestSetup, :(using FractalTools); recursive=true)

makedocs(;
    modules=[FractalTools],
    format=Documenter.HTML(),
    pages=[
        
        "Notes" => [
            "1D Fractal Interpolation" => "manual/one_dimensional_interpolation.md",
            ],
        "API" => "index.md",
        ],
    repo="https://github.com/zekeriyasari/FractalTools.jl/blob/{commit}{path}#L{line}",
    sitename="FractalTools.jl",
    authors="Zekeriya Sarı, Gizem Kalender"
)

deploydocs(;
    repo="github.com/zekeriyasari/FractalTools.jl",
)
