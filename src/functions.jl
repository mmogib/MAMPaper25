
function SD1(∇f, arg...; kw...)
	function gs(_x0, x, args...; kwargs...)
		-∇f(x)
	end
	return gs
end

function SD2(∇f, M, arg...; αmin = 1e-3, αmax = 1e3, kw...)
	UD = Uniform(αmin, αmax)
	α0 = rand(UD)
	function spectral(x0, x, d, args...; α = 1, iter = 1, kwargs...)
		if iter == 1
			return -α0 * ∇f(x0)
		else
			gradf = ∇f(x)
			sk = parallel_transport_to(M, x0, d, x)

			yk = gradf + (1 / α) * sk
			norm_sk_2 = inner(M, x, sk, sk)
			inner_sk_yk = inner(M, x, yk, sk)

			αk1 = max(αmin, min(norm_sk_2 / inner_sk_yk, αmax))
			d = -αk1 * gradf
			return d
		end
	end
	return spectral
end



function St1(M, arg...; λL = x -> x, kw...)
	function line_search(x, args...; kwargs...)
		t = λL(x)
		# UD = Uniform(1e-9, t)
		# rand([0.001, 0.01, 0.0001]) * rand(UD), 0, 0
		rand([0.001, 0.01, 0.0001]) * t, 0, 0
	end
	return line_search
end

function St2(M, f, arg...; kw...)
	# α = 3.0,
	function ArmijoLS(x, d, args...; ν = 0.0, ρ = 0.5, β = 0.8, max_iter = 1000, kwargs...)
		grad_inner = inner(M, x, -d, d)
		iter = 1
		fx = f(x)
		f_evals = 1
		while true
			if iter > max_iter
				println("Inner Iters Reached")
				return ρ, max_iter, f_evals
			end
			τ = β^iter
			xtrial = exp(M, x, τ * d)
			LHS = f(xtrial)
			f_evals += 1
			RHS = fx + ρ * τ * grad_inner + ν
			if LHS <= RHS
				return τ, iter, f_evals
			end
			iter += 1
		end
	end
	return ArmijoLS
end

function St3(M, f, arg...; kw...)
	# lambda_max value is  3.9 or 1.5 .
	# α = 3.0,
	function AdaptativeStepsize(x, d, args...; ν = 0.0, ρ = 0.5, β = 0.8, max_iter = 1000, kwargs...)
		UD = Uniform(1e-4, 3.9)
		λ0 = rand(UD)
		grad_inner = inner(M, x, -d, d)
		iter = 1
		fx = f(x)
		f_evals = 1
		while true
			if iter > max_iter
				println("Inner Iters Reached")
				return ρ, max_iter
			end
			τ = λ0 * β^iter
			xtrial = exp(M, x, τ * d)
			LHS = f(xtrial)
			f_evals += 1
			RHS = fx + ρ * τ * grad_inner + ν
			if LHS <= RHS
				return τ, iter, f_evals
			end
			λ0 = λ0 * β^(iter - 1)
			iter += 1
		end
	end
	return AdaptativeStepsize
end

function St4(M, f, arg...; kw...)
	# lambda_max value is  3.9 or 1.5 .
	# α = 3.0,
	function AdaptativeStepsize(x, d, args...; ν = 0.0, ρ = 0.5, β = 0.8, max_iter = 1000, kwargs...)
		UD = Uniform(1e-4, 1.5)
		λ0 = rand(UD)
		grad_inner = inner(M, x, -d, d)
		iter = 1
		fx = f(x)
		f_evals = 1
		while true
			if iter > max_iter
				println("Inner Iters Reached")
				return ρ, max_iter
			end
			τ = λ0 * β^iter
			xtrial = exp(M, x, τ * d)
			LHS = f(xtrial)
			f_evals += 1
			RHS = fx + ρ * τ * grad_inner + ν
			if LHS <= RHS
				return τ, iter, f_evals
			end
			# λ0 = λ0 * β^(iter - 1)
			iter += 1
		end
	end
	return AdaptativeStepsize
end
function MP0(arg...; kwargs...)
	function MP(args...; kwargs2...)
		0.0, 2
	end
	MP
end

function MP1(f; rng = Xoshiro(2025), kw...)
	function MP(x_0, x_1, args...; ν0 = 1, δ0 = 0.89, kwargs...)
		UD = Uniform(δ0, 1)
		δk = rand(rng, UD)
		νk = (1 - δk) * (f(x_0) - f(x_1) + ν0)
		νk, 2
	end
	MP
end

function MP2(arg...; max_iters = 1000, rng = Xoshiro(2025), kw...)
	seq = [1 / t^2 for t in 1:max_iters]
	function MP(args...; kwargs...)
		rand(rng, seq), 0
	end
	return MP
end
function MP3(arg...; max_iters = 1000, rng = Xoshiro(2025), kw...)
	seq = [1 / (t + 1) for t in 1:max_iters]
	function MP(args...; kwargs...)
		rand(rng, seq), 0
	end
	MP
end
function MP4(arg...; max_iters = 1000, rng = Xoshiro(2025), kw...)
	seq = [1 / t for t in 1:max_iters]
	function MP(args...; kwargs...)
		rand(rng, seq), 0
	end
	MP
end
function MP5(arg...; max_iters = 1000, rng = Xoshiro(2025), kw...)
	seq = [1 / (1 + t^2) for t in 1:max_iters]
	function MP(args...; gradnorm = 1, kwargs...)
		0.5 * gradnorm * rand(rng, seq), 0
	end
	MP
end

function GMP1(problem::ManifoldProblem, x0;
	SD::Function = SD2,
	LS::Function = St2,
	MP::Function = MP1,
	ν0 = 0.5,
	max_iters = 1000,
	backtracking_iters = 2000,
	ϵ = 1e-5,
	rng = Xoshiro(2025),
	trycatch = True,
)
	M, f, ∇f, λL = problem.M, problem.f, problem.∇f, problem.λL
	if trycatch
		try

			sol, time = @timed _GMP1(M, x0, SD(∇f, M), LS(M, f; λL = λL), MP(f; max_iters = 1000, rng = rng), ν0, max_iters, backtracking_iters, ϵ)
			x, stat, msg, min_norm = sol
			timed_stat = Stats(stat, time)
			return Solution(
				x, problem, String(Symbol(SD)), String(Symbol(LS)), String(Symbol(MP)), timed_stat, msg, min_norm,
			)
		catch
			return Solution(
				nothing, problem, String(Symbol(SD)), String(Symbol(LS)), String(Symbol(MP)), Stats(), FailureMessage("Err"), Inf64,
			)
		end
	else
		sol, time = @timed _GMP1(M, x0, SD(∇f, M), LS(M, f; λL = λL), MP(f; max_iters = 1000, rng = rng), ν0, max_iters, backtracking_iters, ϵ)
		x, stat, msg, min_norm = sol
		timed_stat = Stats(stat, time)
		return Solution(
			x, problem, String(Symbol(SD)), String(Symbol(LS)), String(Symbol(MP)), timed_stat, msg, min_norm,
		)
	end
end
function _GMP1(M, x0, SD::Function, LS::Function, MP::Function, ν0, max_iters, backtracking_iters, ϵ)
	iter = 1
	x = copy(x0)
	bk_iter_total = 0
	ν = ν0
	d = rand(M)
	α = 1
	grad_evals = 0
	f_evals = 0
	# ProgressBar = Progress(max_iters, dt = 1.0, desc = "Working ...", color = :magenta)
	ProgressBar = ProgressUnknown(desc = "Working ...", color = :magenta, spinner = true)
	while true
		d = SD(x0, x, d; α = α, iter = iter)
		grad_evals += 1
		norm_x_d = norm(M, x, d)
		# iter % 10 == 0 && printstyled("$(String(Symbol(MP))) Itreation $iter, norm value = $norm_x_d\n", color = :green)
		if norm_x_d <= ϵ
			finish!(ProgressBar)
			stat = Stats(iterations = iter, backtracking_iterations = bk_iter_total, function_evaluations = f_evals, grad_evaluations = grad_evals)
			msg = SuccessMessage("Converged")
			return x, stat, msg, norm_x_d
		end
		if iter >= max_iters
			finish!(ProgressBar, spinner = '❌')
			stat = Stats(iterations = iter, backtracking_iterations = bk_iter_total, function_evaluations = f_evals, grad_evaluations = grad_evals)
			msg = FailureMessage("Max Iterations Reached")
			return x, stat, msg, norm_x_d
		end
		α, bk_iter, ls_f_evals = LS(x, d; ν = ν, max_iter = backtracking_iters)
		# printstyled("$(String(Symbol(MP))) Itreation $iter, alpha = $α\n", color = :magenta)
		bk_iter_total += bk_iter
		f_evals += ls_f_evals

		x, x0 = exp(M, x, α * d), x
		ν, mp_f_evals = MP(x0, x; gradnorm = norm_x_d, iter = iter)
		f_evals += mp_f_evals
		# printstyled("$(String(Symbol(MP))) Itreation $iter, ν = $ν\n", color = :magenta)
		update!(ProgressBar, iter, showvalues = [("Iterations", iter), ("||x-d||", norm_x_d), ("α", α)], spinner = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
		iter += 1
	end
end

# "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
# "🌑🌒🌓🌔🌕🌖🌗🌘"


# function MyRiemanniagradient(M, f; ∇f = nothing)
# 	# n, = representation_size(M)
# 	x0 = rand(M)
# 	arm_stepsize = Manopt.ArmijoLinesearch(M)
# 	# arm_stepsize = Armijo2Linesearch(M)
# 	opt = OptimizationManopt.GradientDescentOptimizer()
# 	prob = if isnothing(∇f)
# 		optf = OptimizationFunction(f, Optimization.AutoZygote())
# 		OptimizationProblem(
# 			optf, x0, []; manifold = M, stepsize = arm_stepsize)
# 	else
# 		optf = OptimizationFunction(f, grad = ∇f)
# 		OptimizationProblem(optf, x0; manifold = M)
# 	end
# 	sol = solve(prob, opt)
# 	return sol
# end
