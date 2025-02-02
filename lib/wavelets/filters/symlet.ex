defmodule Wavelets.Filters.Symlet do
  @moduledoc """
  Implementation of Symlet wavelet filters with proper normalization
  """

  alias Wavelets.Filter

  def get(n) when n in 2..10 do
    case n do
      4 ->
        # Normalization factor
        h0 = 1 / :math.sqrt(2)

        coeffs = [
          -0.0757657147893037,
          -0.0296355276459541,
          0.4976186676324578,
          0.8037387518059161,
          0.2978577956055422,
          -0.0992195435769354,
          -0.0126039672622612,
          0.0322231006040427
        ]

        normalized_coeffs = Enum.map(coeffs, &(&1 * h0))

        %Filter{
          name: "sym4",
          family: :symlet,
          vanishing_moments: 4,
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
          support_width: 8
        }

      8 ->
        h0 = 1 / :math.sqrt(2)

        coeffs = [
          -0.0033824159510061256,
          -0.0005421323317911481,
          0.03169508781149298,
          -0.007607487324917605,
          -0.1432942383758130,
          0.0612733590679088,
          0.4813596512583722,
          0.7771857517005235,
          0.3644418948353314,
          -0.0519458381078751,
          -0.0272190299168137,
          0.0491371796734768,
          0.0038087520138944,
          -0.0149522583367926,
          -0.0003029205147213,
          0.0018899503327594
        ]

        normalized_coeffs = Enum.map(coeffs, &(&1 * h0))

        %Filter{
          name: "sym8",
          family: :symlet,
          vanishing_moments: 8,
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
          support_width: 16
        }
    end
  end
end
