using Documenter, FractalTools

DocMeta.setdocmeta!(FractalTools, :DocTestSetup, :(using FractalTools); recursive=true)

makedocs(;
    modules=[FractalTools],
    format=Documenter.HTML(),
    pages=[
        "API" => "index.md",
        # "Notes" => [
        #     "Interpolation" => [
        #         "1D Fractal Interpolation" => "manual/interpolation/one_dimensional_interpolation.md",
        #         "2D Fractal Interpolation" => "manual/interpolation/two_dimensional_interpolation.md",
        #         "1D Hidden Fractal Interpolation" => "manual/interpolation/one_dimensional_hidden_interpolation.md",
        #         "2D Hidden Fractal Interpolation" => "manual/interpolation/two_dimensional_hidden_interpolation.md",
        #         ],
        #     "Integration" => [
        #         "1D Fractal Integration" => "manual/integration/one_dimensional_integration.md",
        #         "2D Fractal Integration" => "manual/integration/two_dimensional_integration.md",
        #         "1D Hidden Fractal Integration" => "manual/integration/one_dimensional_hidden_integration.md",
        #         "2D Hidden Fractal Integration" => "manual/integration/two_dimensional_hidden_integration.md",
        #         ],
        # ],
    ],
    repo="https://github.com/zekeriyasari/FractalTools.jl/blob/{commit}{path}#L{line}",
    sitename="FractalTools.jl",
    authors="Zekeriya Sarı, Gizem Kalender"
)

deploydocs(;
    repo="github.com/zekeriyasari/FractalTools.jl",
)
