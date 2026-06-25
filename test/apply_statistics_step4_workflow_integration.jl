#!/usr/bin/env julia

# Applies the Statistical Step 4 integration to the legacy workflow.
# The operation is idempotent and aborts if the expected source blocks differ.

path = joinpath(@__DIR__, "..", "src", "workflow.jl")
text = read(path, String)

old_return = """    return (data=output_data, cell=cell, log=log_data, dir=ldir)
end
"""

new_return = """    statistics_report = analyze_md_cycle_statistics(cycle, ldir)

    return (
        data=output_data,
        cell=cell,
        log=log_data,
        dir=ldir,
        statistics=statistics_report.statistics,
        statistics_report=statistics_report,
    )
end
"""

if occursin(old_return, text)
    text = replace(text, old_return => new_return; count=1)
elseif !occursin("statistics_report = analyze_md_cycle_statistics(cycle, ldir)", text)
    error("Could not locate the step_npt! return block in src/workflow.jl")
end

old_summary = """        summary = step_analyze!(wp, cycle, npt.cell, gcmc.output, ref_cell, npt.data)
        push!(history, summary)
"""

new_summary = """        summary = step_analyze!(wp, cycle, npt.cell, gcmc.output, ref_cell, npt.data)
        summary[\"md_statistics_csv\"] = npt.statistics_report.csv_path
        summary[\"md_statistics_json\"] = npt.statistics_report.json_path
        summary[\"md_statistics_valid\"] = npt.statistics.valid
        push!(history, summary)
"""

if occursin(old_summary, text)
    text = replace(text, old_summary => new_summary; count=1)
elseif !occursin("summary[\"md_statistics_csv\"]", text)
    error("Could not locate the legacy workflow summary block in src/workflow.jl")
end

write(path, text)
println("Integrated MD statistics reporting into src/workflow.jl")
