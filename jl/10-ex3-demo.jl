# 拡張ユークリッドアルゴリズムのデモ（練習問題3）
using GaloisFields
using Polynomials
using CodingTheoryUtils

const PRIMITIVE_POLY = "x^4 + x + 1"
const T = 3
const GF, α = GaloisField(F2, Polynomial(string2coefvec(PRIMITIVE_POLY), "α"))

function pretty_coeffs(poly)
    return join(gf_pretty.(poly.coeffs, α), ", ")
end

function show_poly(label, poly)
    deg_val = degree(poly)
    deg_str = deg_val == -Inf ? "-Inf" : string(deg_val)
    println("$label (deg $deg_str): [$(pretty_coeffs(poly))]")
end

function syndrome_polynomial(rx, t)
    coeffs = zeros(typeof(α), 2t)
    for i in 1:2t
        coeffs[i] = rx(α^i)
    end
    return Polynomial(coeffs, "x")
end

function euclid_with_trace(Sx, t)
    base_type = typeof(α)
    head_coeffs = zeros(base_type, 2t + 1)
    head_coeffs[end] = base_type(1)
    r_prev = Polynomial(head_coeffs, "x")
    r_curr = Sx
    w_prev = Polynomial([base_type(0)], "x")
    w_curr = Polynomial([base_type(1)], "x")
    println("=== 初期化 ===")
    show_poly("r₀(x) = x^{2t}", r_prev)
    show_poly("r₁(x) = S(x)", r_curr)
    show_poly("w₀(x)", w_prev)
    show_poly("w₁(x)", w_curr)
    step = 1
    while degree(r_curr) ≥ t && !iszero(r_curr)
        q, r = divrem(r_prev, r_curr)
        println("\n=== Step $step: r_{step-1}(x) ÷ r_{step}(x) ===")
        show_poly("r_{step-1}(x)", r_prev)
        show_poly("r_{step}(x)", r_curr)
        show_poly("q_$step(x)", q)
        show_poly("r_{step+1}(x)", r)
        r_prev, r_curr = r_curr, r
        w_prev, w_curr = w_curr, w_prev - q * w_curr
        show_poly("w_{step}(x)", w_prev)
        show_poly("w_{step+1}(x)", w_curr)
        step += 1
    end
    σ = w_curr / w_curr.coeffs[1]
    println("\n=== 正規化 ===")
    show_poly("σ(x) before norm", w_curr)
    show_poly("σ(x)", σ)
    Ω = r_curr
    show_poly("Ω(x)", Ω)
    return σ, Ω
end

function find_error_positions(σ)
    positions = Int[]
    println("\n=== 誤り位置探索 ===")
    for exp in 0:order(α)-1
        val = σ(α^exp)
        if val == zero(val)
            pos = mod(-exp, order(α))
            push!(positions, pos)
            println("σ(α^$exp) = 0 → X = α^$(pos) → ビット位置 $(pos + 1)")
        end
    end
    return positions
end

function main()
    println("=== 練習問題3 デモ ===")
    println("GF(2^4) 生成多項式: $PRIMITIVE_POLY, t = $T")
    rx = string2F2poly("x + x^3 + x^5")
    println("\n受信多項式 r(x): $rx")
    Sx = syndrome_polynomial(rx, T)
    println("\n=== シンドローム計算 ===")
    show_poly("S(x)", Sx)
    σ, Ω = euclid_with_trace(Sx, T)
    positions = find_error_positions(σ)
    println("\n推定誤り位置指数 (0-index): $positions")
end

main()
