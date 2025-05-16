function generate_positive_definite_matrix(n::Int; min_eigen = 1e-5, max_eigen = 20.0, rng = Xoshiro(2025))
	# Step 1: Generate random eigenvalues in (min_eigen, max_eigen)
	UD = Uniform(min_eigen, max_eigen)
	eigenvalues = rand(UD, n)

	# Step 2: Generate a random orthogonal matrix Q
	Q, _ = qr(randn(rng, n, n))  # QR decomposition gives an orthonormal matrix Q
	# Q = BigFloat.(Q, RoundUp, precision = 64)
	# Step 3: Construct A = QΛQ^T
	A = Q * Diagonal(eigenvalues) * Q'
	A = (A + A') / 2
	return A
end

function generate_positive_definite_matrix(n::Int, m::Int; kwargs...)
	return map(_ -> generate_positive_definite_matrix(n; kwargs...), 1:m)
end
# SymmetricPositiveDefinite (Example 4.4.)
# problem == 3 
function spd_problem1(n, ars...; rng = Xoshiro(2025), α = 0.0, β = 10.0)
	UD = Uniform(α, β)
	a = rand(rng, UD)
	b = rand(rng, UD)

	f(X) = begin
		det_x = det(X)
		log_det_x = log(det_x)
		a .* log_det_x^2 - b * log_det_x
	end
	M = SymmetricPositiveDefinite(n)
	∇f(x) = (2 * a * log(det(x)) - b) * x
	λL = (args...; kwargs...) -> inv(2 * sqrt(n) * a)
	x0 = rand(rng, M)
	return ManifoldProblem(name = "Problem 1 a=$a, b = $b", M = M, f = f, ∇f = ∇f, n = n, λL = λL, x0 = x0)
end

# SymmetricPositiveDefinite (Example 4.5.)
# problem == 4
function spd_problem2(n, ars...; rng = Xoshiro(2025), α = 0.0, β = 10.0)
	UD = Uniform(α, β)
	a = rand(rng, UD)
	b = rand(rng, UD)
	d = rand(rng, UD)
	c = rand(rng) * a * d
	M = SymmetricPositiveDefinite(n)
	f(x) = begin
		x = real(x)
		detx = det(x)
		f = a * log(detx^d + b) - c * log(detx)
		f
	end
	∇f(x) = begin
		detx = det(x)
		g = (a * d * detx^d / (detx^d + b) - c) * inv(x)
		g = real(g)
		g = (g + g') / 2
		g
	end
	λL = (args...; kwargs...) -> inv(n * a * d^2)
	# x0 = generate_positive_definite_matrix(n; rng = rng)
	x0 = rand(rng, M)
	return ManifoldProblem(name = "Problem 2 a=$a, b = $b, c=$c, d=$d", M = M, f = f, ∇f = ∇f, n = n, λL = λL, x0 = x0)
end

# SymmetricPositiveDefinite
# 5.2.1. The center of mass on the SPD matrices cone.
# problem == 6
function spd_problem3(n::Int; m::Int = 5, rng = Xoshiro(2025))
	As = generate_positive_definite_matrix(n, m; rng = rng)
	M = SymmetricPositiveDefinite(n)
	f(x) = begin
		invsqrtx = sqrt(x) / I
		sum(log(sum((invsqrtx * A * invsqrtx) .^ 2)) for A in As)
	end
	∇f(x) = begin
		sqrtx = sqrt(x)
		As_inv = map(A -> A / I, As)
		sum(sqrtx * log(sqrtx * Ainv * sqrtx) * sqrtx for Ainv in As_inv)
	end
	x0 = exp(mean(map(log, As)))
	λL = (args...; kwargs...) -> begin
		# invsqrtx = inv(sqrt(x))
		# radius = maximum(log(invsqrtx * norm(A) * invsqrtx) for A in As)
		radius = maximum(distance(M, x0, A) for A in As)
		radius4 = 4 * radius
		t = inv(radius4 * coth(radius4))
		return 1.99 * t

	end
	return ManifoldProblem(name = "spd_problem3 $(length(As)) As", M = M, f = f, ∇f = ∇f, n = n, λL = λL, x0 = x0)
end


# Positive Orthant (Example 4.1.)
# problem == 1
function po_problem1(n, ars...; rng = Xoshiro(2025), α = 0.0, β = 10.0, l = 0.0, u = 20.0)
	UD = Uniform(α, β)
	a = rand(rng, UD, n)
	b = rand(rng, UD, n)
	d = rand(rng, UD, n)
	c = 1.1 * a + 3.9 * a
	f(x) = sum(-a .* exp.(-b .* x) .+ c .* log.(x) .^ 2 + d .* log.(x))
	∇f(x) = a .* b .* exp.(-b .* x) + 2 * c .* log.(x) ./ x + d ./ x
	M = PositiveVectors(n)
	λL = (args...; kwargs...) -> inv(sum((a + 2c) .^ 2))
	x0 = rand(rng, Uniform(l, u), n)
	ManifoldProblem(name = "PO Problem 1", M = M, f = f, ∇f = ∇f, n = n, λL = λL, x0 = x0)
end

# Positive Orthant (Example 4.2.)
# problem == 2
function po_problem2(n, ars...; rng = Xoshiro(2025), α = 0.0, β = 10.0, l = 0.0, u = 20.0)
	UD = Uniform(α, β)
	a = rand(rng, UD, n)
	b = rand(rng, UD, n)
	UD = Uniform(α + 2, β)
	d = rand(rng, UD, n)
	c = rand(rng) * a .* d
	f(x) = sum(a .* log.(x .^ d .+ b) .- c .* log.(x))
	∇f(x) = (a .* d .* (x .^ (d .- 1))) ./ (x .^ d .+ b) - c ./ x
	M = PositiveVectors(n)
	λL = (args...; kwargs...) -> inv(sum(a .^ 2 .* (d .^ 4)))
	x0 = rand(rng, Uniform(l, u), n)
	ManifoldProblem(name = "PO Problem 2", M = M, f = f, ∇f = ∇f, n = n, λL = λL, x0 = x0)
end

# Positive Orthant 5.2.2. The center of mass on the positive orthant.
# problem == 5 
function po_problem3(n, ars...; m::Int = 5, rng = Xoshiro(2025), l = 1e-6, u = 100.0)
	M = PositiveVectors(n)
	UD = Uniform(l, u)
	ws = map(_ -> rand(rng, UD, n), 1:m)
	f(x) = begin
		0.5 * sum(sum((log.(w ./ x)) .^ 2) for w in ws)
	end
	∇f(x) = begin
		map(i -> x[i] * sum(log(x[i] / ws[j][i]) for j in 1:m), 1:n)
	end
	x0 = rand(rng, UD, n)
	λL = (args...; kwargs...) -> begin
		# invsqrtx = inv(sqrt(x))
		# radius = maximum(log(invsqrtx * norm(A) * invsqrtx) for A in As)
		radius = maximum(distance(M, x0, w) for w in ws)
		radius4 = 4 * radius
		t = inv(radius4 * coth(radius4))
		return 1.99 * t

	end

	ManifoldProblem(name = "PO Problem 2", M = M, f = f, ∇f = ∇f, n = n, λL = λL, x0 = x0)
end




