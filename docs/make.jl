using Documenter, StoMPS

makedocs(;
    modules=[StoMPS],
    sitename="StoMPS.jl",
    authors="Jianxin Gao, Qiaoyi Li, Kunyang Du",
    format=Documenter.HTML(edit_link="main"),
    pages=[
        "Home" => "index.md",
        "Execution and output" => "design/runtime.md",
        "Annealed MPI resampling" => "design/resampling.md",
    ],
)
