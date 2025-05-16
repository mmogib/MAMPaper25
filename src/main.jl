include("depts.jl")
include("types.jl")
include("utils.jl")
include("functions.jl")
include("problems.jl")
# Example 1
begin
	rng = Xoshiro(2025)
	backtracking_iters = 1000
	max_iters = 1000
	ϵ = 1e-4
	ν0 = 0.5
	trycatch = true
	szs = [20, 50, 100, 150, 200]
	# szs = [20]
	SearchDirections = [SD1, SD2]
	# SD2,St2,MP2
	StepSizes = [St1, St2, St3, St4]
	# StepSizes = [St2]
	MPs = [MP0, MP1, MP2, MP3, MP4, MP5]
	# MPs = [MP2]
	problems = [spd_problem1, spd_problem2, po_problem1, po_problem2]
	# problems = [spd_problem2]
	all_sols = []
	total_number_of_experiments = length(SearchDirections) * length(StepSizes) * length(MPs) * length(szs) * length(problems)
	counter = 0
	for problem_fn in problems
		for n in szs
			sols = Vector{Solution}(undef, length(SearchDirections) * length(StepSizes) * length(MPs))
			sols_str = Vector{String}(undef, length(SearchDirections) * length(StepSizes) * length(MPs))



			problem = problem_fn(n; rng = rng)


			# x0 = rand(rng, problem.M) #
			x0 = problem.x0
			i = 0
			for sd in SearchDirections
				for st in StepSizes
					for mp in MPs
						i += 1
						counter += 1
						sol = GMP1(problem, x0;
							SD = sd,
							LS = st,
							MP = mp,
							backtracking_iters = backtracking_iters,
							max_iters = max_iters,
							ϵ = ϵ,
							ν0 = ν0,
							rng = rng,
							trycatch = trycatch,
						)
						# printstyled(join(String.(Symbol.([sd, st, mp])), " :: "), " value ", problem.f(sol_mp5[1]), " Iterations: ", sol_mp5[2], "\n", color = :green)
						color = (isa(sol.message, SuccessMessage)) ? :green : :red
						@printf "=========================================\n"
						@printf "            Solving (%s)     \n" String(Symbol(problem_fn))
						@printf "            Solution Log (n==%d)     \n" n
						@printf "            Experiment (%d / %d)     \n" counter total_number_of_experiments
						@printf "=========================================\n"
						printstyled(sol, "\n", color = color)
						sols[i] = sol
						sols_str[i] = "$sol"
					end
				end
			end
			push!(all_sols, ("$n", sols, sols_str))
		end
	end
	df = solutions_to_dataframe(vcat(map(x -> x[2], all_sols)...))
	save_plots(df, algorithms = Symbol.(StepSizes), label = :step_size, folder_name = "results/profiles", file_prefix = "step_size_")
	save_plots(df, algorithms = Symbol.(MPs), label = :nonmonotone_parameter, folder_name = "results/profiles", file_prefix = "nonmonoton_params_")
	save_plots(df, algorithms = Symbol.(SearchDirections), label = :search_direction, folder_name = "results/profiles", file_prefix = "search_direction_")

	filename = "results/solutions.csv"
	savedf(filename, df)
	filename = "results/solutions.txt"
	save_solutions(filename, map(x -> (x[1], x[3]), all_sols))
	filename = "results/solutions.xlsx"
	savedf(filename, df)
end
