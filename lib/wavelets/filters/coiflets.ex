defmodule Wavelets.Filters.Coiflets do
  @moduledoc """
  Implementation of Coiflet wavelet filters with proper normalization
  """

  alias Wavelets.Filter

  def get(n) when n in [1, 2, 3, 4, 5] do
    case n do
      1 ->
        h0 = 1 / :math.sqrt(2)

        coeffs = [
          -0.0156557281,
          -0.0727326195,
          0.3848648469,
          0.8525720202,
          0.3378976625,
          -0.0727326195
        ]

        normalized_coeffs = Enum.map(coeffs, &(&1 * h0))

        %Filter{
          name: "coif1",
          family: :coiflet,
          vanishing_moments: 2,
          decomposition_low_pass: normalized_coeffs,
          decomposition_high_pass:
            normalized_coeffs
            |> Enum.reverse()
            |> Enum.with_index()
            |> Enum.map(fn {c, i} -> c * :math.pow(-1, i) end),
          reconstruction_low_pass: normalized_coeffs,
          reconstruction_high_pass:
            normalized_coeffs
            |> Enum.with_index()
            |> Enum.map(fn {c, i} -> c * :math.pow(-1, i + 1) end),
          support_width: 6
        }
    end
  end
end
