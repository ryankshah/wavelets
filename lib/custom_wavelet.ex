defmodule Wavelets.CustomWavelet do
  @moduledoc """
  Tools for designing custom wavelets and verifying their properties with corrected vanishing moments calculation
  """

  alias Wavelets.Filter

  @doc """
  Creates a custom wavelet filter from given coefficients
  """
  def create(decomp_low, decomp_high, recon_low, recon_high, opts \\ []) do
    with :ok <-
           verify_length_match(decomp_low, decomp_high, recon_low, recon_high) do
      filter = %Filter{
        name: Keyword.get(opts, :name, "custom"),
        family: :custom,
        vanishing_moments: compute_vanishing_moments(decomp_high),
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
    # Start from 0 instead of 1 for correct vanishing moments calculation
    0..9
    |> Enum.find(0, fn m ->
      moments = compute_moments(high_pass, m)
      not almost_zero(moments)
    end)
  end

  defp compute_moments(coeffs, order) do
    coeffs
    |> Enum.with_index()
    |> Enum.map(fn {c, i} -> c * :math.pow(i, order) end)
    |> Enum.sum()
  end

  defp almost_zero(x, tolerance \\ 1.0e-10), do: abs(x) < tolerance
end
