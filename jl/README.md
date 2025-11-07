# 2025年度千葉大学集中講義 Julia Scripts

## 10-ex3-demo.jl の実行環境

### 必須ライブラリ
- [GaloisFields.jl](https://github.com/tkluck/GaloisFields.jl)
- [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl)
- [CodingTheoryUtils.jl](https://github.com/uchkw/CodingTheoryUtils.jl)

### パッケージ導入手順
- Julia 1.10 以降を起動し，Pkg モードに入る（`]` キー）
- 上記３つのライブラリを通常の General レジストリから取得する
  ```julia
  (@v1.10) pkg> add GaloisFields Polynomials CodingTheoryUtils
  ```

### スクリプトの実行
- 本ディレクトリで `julia 10-ex3-demo.jl` を走らせるだけで，ガロア体 GF(2^4) 上の練習問題 3（拡張ユークリッド法）のトレース出力が得られる
