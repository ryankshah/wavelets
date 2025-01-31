# lib/wavelets/custom_wavelet.ex
defmodule Wavelets.CustomWavelet do
  @moduledoc """
  Tools for designing custom wavelets and verifying their properties
  """

  alias Wavelets.Filter

  @doc """
  Creates a custom wavelet filter from given coefficients.
  Verifies orthogonality and other wavelet properties.
  """
  @spec create(
          list(number),
          list(number),
          list(number),
          list(number),
          keyword()
        ) ::
          {:ok, Filter.t()} | {:error, String.t()}
  def create(decomp_low, decomp_high, recon_low, recon_high, opts \\ []) do
    with :ok <-
           verify_length_match(decomp_low, decomp_high, recon_low, recon_high),
         :ok <-
           verify_perfect_reconstruction(
             decomp_low,
             decomp_high,
             recon_low,
             recon_high
           ),
         :ok <- verify_orthogonality(decomp_low, decomp_high),
         vanishing_moments = compute_vanishing_moments(decomp_high) do
      filter = %Filter{
        name: Keyword.get(opts, :name, "custom"),
        family: :custom,
        vanishing_moments: vanishing_moments,
        decomposition_low_pass: decomp_low,
        decomposition_high_pass: decomp_high,
        reconstruction_low_pass: recon_low,
        reconstruction_high_pass: recon_high,
        support_width: length(decomp_low)
      }

      {:ok, filter}
    end
  end

  @doc """
  Creates an orthogonal wavelet filter from low-pass filter coefficients.
  High-pass filter is automatically generated to ensure orthogonality.
  """
  @spec create_orthogonal(list(number), keyword()) ::
          {:ok, Filter.t()} | {:error, String.t()}
  def create_orthogonal(low_pass, opts \\ []) do
    # Generate quadrature mirror filter for high-pass
    high_pass = generate_qmf(low_pass)

    create(low_pass, high_pass, low_pass, high_pass, opts)
  end

  @doc """
  Creates a biorthogonal wavelet filter pair from decomposition and reconstruction
  low-pass filters. High-pass filters are automatically generated.
  """
  @spec create_biorthogonal(list(number), list(number), keyword()) ::
          {:ok, Filter.t()} | {:error, String.t()}
  def create_biorthogonal(decomp_low, recon_low, opts \\ []) do
    decomp_high = generate_biorthogonal_high_pass(decomp_low, recon_low)
    recon_high = generate_biorthogonal_high_pass(recon_low, decomp_low)

    create(decomp_low, decomp_high, recon_low, recon_high, opts)
  end

  # Verification functions
  defp verify_length_match(d_low, d_high, r_low, r_high) do
    lengths = [length(d_low), length(d_high), length(r_low), length(r_high)]

    if Enum.uniq(lengths) |> length() == 1 do
      :ok
    else
      {:error, "Filter lengths must match"}
    end
  end

  defp verify_perfect_reconstruction(d_low, d_high, r_low, r_high) do
    # Check perfect reconstruction conditions
    # H0(z)G0(z) + H1(z)G1(z) = 2
    # H0(-z)G0(z) + H1(-z)G1(z) = 0

    coeffs_sum =
      convolve(d_low, r_low)
      |> Enum.zip(convolve(d_high, r_high))
      |> Enum.map(fn {a, b} -> a + b end)

    expected_sum = [2.0] ++ List.duplicate(0.0, length(coeffs_sum) - 1)

    if almost_equal(coeffs_sum, expected_sum) do
      :ok
    else
      {:error, "Perfect reconstruction condition not satisfied"}
    end
  end

  defp verify_orthogonality(low_pass, high_pass) do
    # Check orthogonality conditions
    # Σ h[n]h[n+2k] = δ[k]
    # Σ g[n]g[n+2k] = δ[k]
    # Σ h[n]g[n+2k] = 0

    if verify_orthogonality_condition(low_pass) and
         verify_orthogonality_condition(high_pass) and
         verify_cross_orthogonality(low_pass, high_pass) do
      :ok
    else
      {:error, "Orthogonality conditions not satisfied"}
    end
  end

  # Helper functions
  defp generate_qmf(coeffs) do
    coeffs
    |> Enum.reverse()
    |> Enum.with_index()
    |> Enum.map(fn {c, i} -> c * :math.pow(-1, i) end)
  end

  defp generate_biorthogonal_high_pass(analysis, synthesis) do
    # Generate high-pass filter for biorthogonal pair
    analysis
    |> Enum.reverse()
    |> Enum.with_index()
    |> Enum.map(fn {c, i} ->
      syn_idx = rem(i, length(synthesis))
      syn_coeff = Enum.at(synthesis, syn_idx)
      c * syn_coeff * :math.pow(-1, i)
    end)
  end

  defp compute_vanishing_moments(high_pass) do
    # Compute number of vanishing moments
    1..10
    |> Enum.find(1, fn m ->
      moments = compute_moments(high_pass, m)
      not almost_zero(moments)
    end)
    |> Kernel.-(1)
  end

  defp compute_moments(coeffs, order) do
    coeffs
    |> Enum.with_index()
    |> Enum.map(fn {c, i} -> c * :math.pow(i, order) end)
    |> Enum.sum()
  end

  defp convolve(a, b) do
    len_out = length(a) + length(b) - 1

    for i <- 0..(len_out - 1) do
      0..min(i, length(a) - 1)
      |> Enum.map(fn j ->
        if i - j < length(b),
          do: Enum.at(a, j) * Enum.at(b, i - j),
          else: 0.0
      end)
      |> Enum.sum()
    end
  end

  defp verify_orthogonality_condition(coeffs) do
    len = length(coeffs)

    0..div(len, 2)
    |> Enum.all?(fn k ->
      sum =
        0..(len - 1)
        |> Enum.map(fn n ->
          c1 = Enum.at(coeffs, n, 0.0)
          c2 = Enum.at(coeffs, n + 2 * k, 0.0)
          c1 * c2
        end)
        |> Enum.sum()

      expected = if k == 0, do: 1.0, else: 0.0
      almost_equal(sum, expected)
    end)
  end

  defp verify_cross_orthogonality(low_pass, high_pass) do
    len = length(low_pass)

    0..div(len, 2)
    |> Enum.all?(fn k ->
      sum =
        0..(len - 1)
        |> Enum.map(fn n ->
          l = Enum.at(low_pass, n, 0.0)
          h = Enum.at(high_pass, n + 2 * k, 0.0)
          l * h
        end)
        |> Enum.sum()

      almost_zero(sum)
    end)
  end

  defp almost_zero(x, tolerance \\ 1.0e-10) do
    abs(x) < tolerance
  end

  defp almost_equal(x, y, tolerance \\ 1.0e-10)

  defp almost_equal(x, y, tolerance) when is_number(x) and is_number(y) do
    abs(x - y) < tolerance
  end

  defp almost_equal(xs, ys, tolerance) when is_list(xs) and is_list(ys) do
    Enum.zip(xs, ys)
    |> Enum.all?(fn {x, y} -> almost_equal(x, y, tolerance) end)
  end
end
