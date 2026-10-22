using Documenter, StoMPS

makedocs(;
    modules=[StoMPS],
    sitename="StoMPS.jl",
    authors="Jianxin Gao, Qiaoyi Li, Kunyang Du",
    format=Documenter.HTML(edit_link="main"),
    pages=[
        "Home" => "index.md",
        "Getting started" => "guide/quickstart.md",
        "Applications" => "guide/applications.md",
        "Parameters" => "guide/parameters.md",
        "Reading results" => "guide/results.md",
        "Convergence and uncertainty" => "guide/validation.md",
        "Troubleshooting" => "guide/troubleshooting.md",
        "Design" => [
            "Paper and notation" => "design/paper.md",
            "Physical conventions" => "design/physical_conventions.md",
            "Tensor conventions" => "design/tensor_conventions.md",
            "Born sampling" => "design/born_sampling.md",
            "Thermal evolution" => "design/finite_temperature.md",
            "Resampling" => "design/resampling.md",
            "Execution and output" => "design/runtime.md",
        ],
        "Library tutorials" => "guide/library.md",
        "API reference" => "api.md",
    ],
)
