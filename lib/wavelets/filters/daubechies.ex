defmodule Wavelets.Filters.Daubechies do
  @moduledoc """
  Implementation of Daubechies wavelet filters
  """

  alias Wavelets.Filter

  @doc """
  Returns coefficients for the Daubechies wavelet family
  """
  def get(vanishing_moments) when vanishing_moments > 0 do
    case vanishing_moments do
      1 ->
        # Haar wavelet coefficients - normalized for proper scaling
        # This makes [4.0, 4.0] -> [8.0] in forward transform
        h0 = 1.0

        %Filter{
          name: "db1",
          family: :daubechies,
          vanishing_moments: 1,
          decomposition_low_pass: [h0, h0],
          decomposition_high_pass: [-h0, h0],
          # Division by 2 for reconstruction
          reconstruction_low_pass: [h0 / 2, h0 / 2],
          reconstruction_high_pass: [h0 / 2, -h0 / 2],
          support_width: 2
        }

      2 ->
        # Daubechies-4 coefficients
        h = [
          (1 + :math.sqrt(3)) / 4,
          (3 + :math.sqrt(3)) / 4,
          (3 - :math.sqrt(3)) / 4,
          (1 - :math.sqrt(3)) / 4
        ]

        %Filter{
          name: "db2",
          family: :daubechies,
          vanishing_moments: 2,
          # Scale up for decomposition
          decomposition_low_pass: Enum.map(h, &(&1 * 2)),
          decomposition_high_pass:
            Enum.zip(h |> Enum.reverse(), [1, -1, 1, -1])
            |> Enum.map(fn {a, b} -> 2 * a * b end),
          # Original coefficients for reconstruction
          reconstruction_low_pass: h,
          reconstruction_high_pass:
            Enum.zip(h |> Enum.reverse(), [1, -1, 1, -1])
            |> Enum.map(fn {a, b} -> a * b end),
          support_width: 4
        }
    end
  end
end
