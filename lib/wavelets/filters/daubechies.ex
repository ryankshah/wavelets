defmodule Wavelets.Filters.Daubechies do
  @moduledoc """
  Implementation of Daubechies wavelet filters with proper normalization
  """

  alias Wavelets.Filter

  @doc """
  Returns coefficients for the Daubechies wavelet family
  """
  def get(vanishing_moments) when vanishing_moments > 0 do
    case vanishing_moments do
      1 ->
        # Haar wavelet coefficients (properly normalized)
        # This makes [4.0, 4.0] -> [4.0]
        h0 = 0.5

        %Filter{
          name: "db1",
          family: :daubechies,
          vanishing_moments: 1,
          decomposition_low_pass: [h0, h0],
          decomposition_high_pass: [-h0, h0],
          # x2 for reconstruction
          reconstruction_low_pass: [2 * h0, 2 * h0],
          reconstruction_high_pass: [2 * h0, -2 * h0],
          support_width: 2
        }

      2 ->
        # Daubechies-4 coefficients (properly normalized)
        h = [
          (1 + :math.sqrt(3)) / (4 * :math.sqrt(2)),
          (3 + :math.sqrt(3)) / (4 * :math.sqrt(2)),
          (3 - :math.sqrt(3)) / (4 * :math.sqrt(2)),
          (1 - :math.sqrt(3)) / (4 * :math.sqrt(2))
        ]

        %Filter{
          name: "db2",
          family: :daubechies,
          vanishing_moments: 2,
          decomposition_low_pass: h,
          decomposition_high_pass:
            Enum.zip(h |> Enum.reverse(), [1, -1, 1, -1])
            |> Enum.map(fn {a, b} -> a * b end),
          # x2 for reconstruction
          reconstruction_low_pass: Enum.map(h, &(&1 * 2)),
          reconstruction_high_pass:
            Enum.zip(h |> Enum.reverse(), [1, -1, 1, -1])
            |> Enum.map(fn {a, b} -> a * b * 2 end),
          support_width: 4
        }
    end
  end
end
