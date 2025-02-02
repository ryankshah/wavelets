defmodule Wavelets.CustomWavelet do
  @moduledoc """
  Tools for designing custom wavelets and verifying their properties
  """

  alias Wavelets.Filter

  @doc """
  Creates a custom wavelet filter from given coefficients
  """
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
         :ok <- verify_orthogonality(decomp_low, decomp_high) do
      filter = %Filter{
        name: Keyword.get(opts, :name, "custom"),
        family: :custom,
        vanishing_moments: compute_vanishing_moments(decomp_high),
        decomposition_low_pass: normalize_coefficients(decomp_low),
        decomposition_high_pass: normalize_coefficients(decomp_high),
        reconstruction_low_pass: normalize_coefficients(recon_low),
        reconstruction_high_pass: normalize_coefficients(recon_high),
        support_width: length(decomp_low)
      }

      {:ok, filter}
    end
  end

  @doc """
  Creates an orthogonal wavelet filter from low-pass filter coefficients
  """
  def create_orthogonal(low_pass, opts \\ []) do
    high_pass = generate_qmf(low_pass)
    create(low_pass, high_pass, low_pass, high_pass, opts)
  end

  @doc """
  Creates a biorthogonal wavelet filter pair
  """
  def create_biorthogonal(decomp_low, recon_low, opts \\ []) do
    decomp_high = generate_biorthogonal_high_pass(decomp_low, recon_low)
    recon_high = generate_biorthogonal_high_pass(recon_low, decomp_low)
    create(decomp_low, decomp_high, recon_low, recon_high, opts)
  end

  defp verify_perfect_reconstruction(d_low, d_high, r_low, r_high) do
    # Check perfect reconstruction conditions
    # H0(z)G0(z) + H1(z)G1(z) = 1

    # First normalize coefficients
    d_low = normalize_coefficients(d_low)
    d_high = normalize_coefficients(d_high)
    r_low = normalize_coefficients(r_low)
    r_high = normalize_coefficients(r_high)

    # Compute convolutions
    conv_low = convolve(d_low, r_low)
    conv_high = convolve(d_high, r_high)

    sum =
      Enum.zip(conv_low, conv_high)
      |> Enum.map(fn {a, b} -> a + b end)

    # First coefficient should be 1, rest should be 0
    [head | tail] = sum

    if almost_equal(head, 1.0, 1.0e-6) and
         Enum.all?(tail, &almost_zero(&1, 1.0e-6)) do
      :ok
    else
      {:error, "Perfect reconstruction condition not satisfied"}
    end
  end

  defp normalize_coefficients(coeffs) do
    norm = coeffs |> Enum.map(&(&1 * &1)) |> Enum.sum() |> :math.sqrt()

    case norm do
      0.0 -> coeffs
      _ -> Enum.map(coeffs, &(&1 / norm))
    end
  end

  defp verify_orthogonality(low_pass, high_pass) do
    if verify_orthogonality_condition(low_pass) and
         verify_orthogonality_condition(high_pass) and
         verify_cross_orthogonality(low_pass, high_pass) do
      :ok
    else
      {:error, "Orthogonality conditions not satisfied"}
    end
  end

  defp verify_length_match(d_low, d_high, r_low, r_high) do
    lengths = [length(d_low), length(d_high), length(r_low), length(r_high)]

    if Enum.uniq(lengths) |> length() == 1,
      do: :ok,
      else: {:error, "Filter lengths must match"}
  end

  defp generate_qmf(coeffs) do
    coeffs
    |> Enum.reverse()
    |> Enum.with_index()
    |> Enum.map(fn {c, i} -> c * :math.pow(-1, i) end)
  end

  defp generate_biorthogonal_high_pass(analysis, synthesis) do
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
    1..10
    |> Enum.find(1, fn m ->
      moments = compute_moments(high_pass, m)
      not almost_zero(moments, 1.0e-6)
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
        if i - j < length(b), do: Enum.at(a, j) * Enum.at(b, i - j), else: 0.0
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
      almost_equal(sum, expected, 1.0e-6)
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

      almost_zero(sum, 1.0e-6)
    end)
  end

  defp almost_zero(x, tolerance), do: abs(x) < tolerance
  defp almost_equal(x, y, tolerance), do: abs(x - y) < tolerance
end
