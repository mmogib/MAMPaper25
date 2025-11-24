	function basenameparts(filename::String)
	parts = split(basename(filename), ".")
	if length(parts) == 1
		return parts[1], "UNKOWN"
	else
		return join(parts[1:end-1], "."), parts[end]
	end
end


function outputfilename(file_name::String; dated::Bool = true)
	fname = basename(file_name)
	root = dirname(file_name)
	root_dir = mkpath(root)
	filename = if dated
		fparts = basenameparts(fname)
		zdt = now(tz"Asia/Riyadh")
		dayfolder = Dates.format(zdt, "yyyy_mm_dd")
		hourfolder = Dates.format(zdt, "HH_MM")
		"$root_dir/$(fparts[1])_$(dayfolder)_$(hourfolder).$(fparts[2])"
	else
		"$root_dir/$fname"
	end

	filename
end

function save_solutions(filename::String, sols)
	header = ["n" "Problem" "Best value" "Search Direction" "Step Size" "Nonmonotone Parameter" "Iterations" "Backtracking Iterations" "Function Evaluations" "Gradient Evaluations" "Message"]
	parts = basenameparts(filename)
	ext = parts[end]
	filename = ext == "UNKOWN" ? "$filename.csv" : filename
	filename = outputfilename(filename)

	if uppercase(ext) in ["CSV", "UNKOWN"]
		for (n, sol_arr) in sols
			open(filename, "w") do io
				writedlm(io, header, ',')
				for sol in sol_arr
					data =
						[
							n sol.problem.name    sol.problem.f(sol.x)    sol.search_direction    sol.step_size    sol.nonmonotone_parameter    sol.stats.iterations    sol.stats.backtracking_iterations     sol.stats.function_evaluations     sol.stats.grad_evaluations sol.message.message
						]
					writedlm(io, data, ',')
				end
			end
		end
	else
		for (n, sol_arr_str) in sols
			open(filename, "w") do io
				for sol in sol_arr_str
					write(io, "n=$n\n" * sol, '\n')
				end
			end
		end
	end
end

function solutions_to_dataframe(solutions::Vector{Solution})
	df = DataFrame(
		dim = [sol.problem.n for sol in solutions],
		problem = [sol.problem.name for sol in solutions],
		best_value = [isnothing(sol.x) ? "" : sol.problem.f(sol.x) for sol in solutions],
		min_norm_value = [sol.norm_x_d for sol in solutions],
		search_direction = [sol.search_direction for sol in solutions],
		step_size = [sol.step_size for sol in solutions],
		nonmonotone_parameter = [sol.nonmonotone_parameter for sol in solutions],
		iterations = [sol.stats.iterations for sol in solutions],
		backtracking_iterations = [sol.stats.backtracking_iterations for sol in solutions],
		function_evaluations = [sol.stats.function_evaluations for sol in solutions],
		grad_evaluations = [sol.stats.grad_evaluations for sol in solutions],
		time = [sol.stats.time for sol in solutions],
		message = [sol.message.message for sol in solutions],
	)
	return df
end

function savedf(filename::String, df::DataFrame; dated::Bool = true)
	fname = outputfilename(filename; dated)
	fparts = basenameparts(fname)
	ext = isnothing(fparts[2]) ? "UNKOWN" : lowercase(fparts[2])
	if ext in ["txt", "csv"]
		@info "Saving data in $(fname)."
		CSV.write(fname, df)
	elseif ext == "xlsx"
		@info "Saving data in $(fname)."
		XLSX.writetable(fname, df; overwrite = true)
	else
		@warn "The file $fname is not supported."
	end
end




function dataFrames2Matrix(data::Vector{DataFrame}, field::Symbol, np::Int64, ns::Int64)

	T = map(data) do df
		values = df[!, field]
		values .|> x -> isempty(x) ? NaN : Float64(x)
	end |> d -> vcat(d...) |> d -> reshape(d, np, ns)
	T
end

function dataFrames2Matrix(data::Vector{DataFrame}, field::Symbol)
	ns = length(data)
	np, = size(data[1])
	dataFrames2Matrix(data, field, np, ns)
end

function readPaperData(dfs::DataFrame, names::Vector{String})
	dfs = transform(dfs, :algorithm => ByRow(strip) => :algorithm)
	data = map(enumerate(names)) do (i, name)
		filter(row -> row[:algorithm] == name, dfs)
	end
	data

end
function readPaperData(input_file::String, names::Vector{String})
	dfs = CSV.File(input_file) |> DataFrame
	readPaperData(dfs, names)
end

function plotPerformanceProfile(T::Matrix{<:Number}, title::String, labels::Vector{String}; colors::Union{ColorPalette, Vector{Symbol}} = palette(:tab10), logscale::Bool = true)
	(w, h) = Plots._plot_defaults[:size]
	(_, _, max_ratio) = performance_profile_data(T)
	p = performance_profile(PlotsBackend(),
		T, labels;
		size = (1.2w, h),
		logscale,
		title = title,
		xlabel = L"\tau",
		ylabel = L"\rho(\tau)",
		legendfontsize = 8,
		linestyle = :dash,
		palette = colors,
		linewidth = 2.5,
		minorgrid = true, leg = :bottomright,
	)
	p = if (logscale)
		xs = 1:ceil(max_ratio + 0.5)
		ys = 0:0.1:1.01
		plot(p,
			xticks = (xs, map(x -> "$(Int(x))", xs)),
			yticks = (ys, map(x -> "$x", ys)),
			# framestyle=:origin
		)
	else
		p
	end
	p

end
function produceProfileImages(input_file::String)
	colors = [:red, :black, :green, :blue, :purple]

	output_folder = "./results/profiles/imgs"


	names = String.(Symbol.([MOPCGM, GMOPCGM, CGPM, GCGPM, STTDFPM, Framework]))
	printstyled("reading stored data in "; color = :blue)
	printstyled("$input_file\n"; color = :blue, bold = true)
	data = readPaperData(input_file, names)
	printstyled("creating plots ...\n"; color = :green)

	ns = length(data)
	np, = size(data[1])
	# 
	plts = map([
		(:iters, "Iterations"),
		(:fun_evaluations, "Function Evaluations"),
		(:time, "CPU Time")]) do (item, title)
		T = dataFrames2Matrix(data, item, np, ns)
		# pretty_table(T)
		printstyled("creating plots for $title \n"; color = :green)
		p = plotPerformanceProfile(T, title, names; colors)
		p, title
	end

	map(plts) do (p, title)
		file_name = replace(title, " " => "_") |> lowercase
		printstyled("saving plots for $title to file \n"; color = :reverse)
		png_file = outputfilename("$output_folder/$(file_name).png"; dated = false)
		savefig(p, png_file)
		# svg_file = outputfilename(file_name, "svg"; root=output_folder, suffix="performance_profiles")
		# savefig(p, svg_file)
		# pdf_file = outputfilename(file_name, "pdf"; root=output_folder, suffix="performance_profiles")
		# savefig(p, pdf_file)
		# eps_file = outputfilename(file_name, "eps"; root=output_folder, suffix="performance_profiles")
		# savefig(p, eps_file)
	end
	printstyled("Finished saving in $output_folder\n"; color = :green)
	plts
end

function save_plots(df::DataFrame; algorithms::Vector{Symbol}, label::Symbol, folder_name::String, file_prefix::String = "")
	colors = [:red,          # 1
		:blue,         # 2
		:green,        # 3
		:darkorange1,  # 4
		:purple,       # 5
		:yellow,       # 6
		:cyan,         # 7
		:magenta,      # 8
		:black,        # 9
		:lime,
	]
	(w, h) = Plots._plot_defaults[:size]
	ys = 0:0.1:1.01

	plts = map([(:iterations, "Iterations"), (:time, "Time"),
		(:function_evaluations, "Function Evaluations"), (:grad_evaluations, "Gradient Evaluations")]) do (data_point, data_header)
		T = Matrix{Float64}(undef, Int(DataFrames.nrow(df) / length(algorithms)), length(algorithms))
		foreach(enumerate(algorithms)) do (i, name)
			d = select(filter(r -> r[label] == String(name), df),
				[data_point, :message] => ByRow((itr, f) -> occursin("Converged", f) ? itr : NaN) => Symbol(data_header))
			T[:, i] = d[!, Symbol(data_header)]
		end
		pdata = performance_profile_data(T)
		max_ratio = pdata[3]
		xs = 1:ceil(max_ratio + 0.5)
		p = performance_profile(PlotsBackend(), T, String.(algorithms);
			title = "$(data_header)",
			logscale = true,
			size = (1.2w, h),
			xlabel = "Performance ratio",
			ylabel = "Solved problems (%)",
			legendfontsize = 8,
			linestyle = :dash,
			palette = colors,
			linewidth = 2.5,
			minorgrid = true, leg = :bottomright,
		)
		plot(p, xticks = (xs, map(x -> "$(Int(x))", xs)),
			yticks = (ys, map(x -> "$x", ys))), data_header
	end

	for (p, name) in plts
		# output_file_svg = outputfilename("$(folder_name)_$name", "svg"; suffix)
		output_file_png = outputfilename("$(folder_name)/$(file_prefix)$(name).png")

		# Plots.svg(p, output_file_svg)
		Plots.png(p, output_file_png)
	end
end

