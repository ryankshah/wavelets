defmodule Wavelets.Filters.Biorthogonal do
  @moduledoc """
  Implementation of Biorthogonal wavelet filters with proper normalization
  """

  alias Wavelets.Filter

  def get(n_decomp, n_recon) do
    case {n_decomp, n_recon} do
      {1, 3} ->
        h0 = 1 / :math.sqrt(2)

        decomp_coeffs = [
          -0.0883883476483184,
          0.0883883476483184,
          0.7071067811865476,
          0.7071067811865476,
          0.0883883476483184,
          -0.0883883476483184
        ]

        recon_coeffs = [
          0.0,
          0.0,
          0.7071067811865476,
          0.7071067811865476,
          0.0,
          0.0
        ]

        norm_decomp = Enum.map(decomp_coeffs, &(&1 * h0))
        norm_recon = Enum.map(recon_coeffs, &(&1 * h0))

        %Filter{
          name: "bior1.3",
          family: :biorthogonal,
          vanishing_moments: 1,
          decomposition_low_pass: norm_decomp,
          decomposition_high_pass:
            recon_coeffs
            |> Enum.reverse()
            |> Enum.with_index()
            |> Enum.map(fn {c, i} -> c * h0 * :math.pow(-1, i) end),
          reconstruction_low_pass: norm_recon,
          reconstruction_high_pass:
            decomp_coeffs
            |> Enum.with_index()
            |> Enum.map(fn {c, i} -> c * h0 * :math.pow(-1, i + 1) end),
          support_width: 6
        }
    end
  end

  def get_reverse(n_decomp, n_recon) do
    %Filter{} = filter = get(n_decomp, n_recon)

    %{
      filter
      | name: "rbio#{n_decomp}.#{n_recon}",
        family: :reverse_biorthogonal,
        decomposition_low_pass: filter.reconstruction_low_pass,
        decomposition_high_pass: filter.reconstruction_high_pass,
        reconstruction_low_pass: filter.decomposition_low_pass,
        reconstruction_high_pass: filter.decomposition_high_pass
    }
  end
end
