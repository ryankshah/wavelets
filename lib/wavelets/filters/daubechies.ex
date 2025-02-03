defmodule Wavelets.Filters.Daubechies do
  @moduledoc """
  Filter implementation with consistent scaling
  """

  alias Wavelets.Filter

  def get(vanishing_moments) when vanishing_moments > 0 do
    case vanishing_moments do
      1 ->
        # Haar wavelet with consistent 0.5 scaling
        %Filter{
          name: "db1",
          family: :daubechies,
          vanishing_moments: 1,
          # Half sum
          decomposition_low_pass: [0.5, 0.5],
          # Half difference
          decomposition_high_pass: [0.5, -0.5],
          # Same for reconstruction
          reconstruction_low_pass: [0.5, 0.5],
          reconstruction_high_pass: [0.5, -0.5],
          support_width: 2
        }

      2 ->
        # Normalized Daubechies-4
        h = [
          1 + :math.sqrt(3),
          3 + :math.sqrt(3),
          3 - :math.sqrt(3),
          1 - :math.sqrt(3)
        ]

        # Basic normalization
        norm = 8.0

        coeffs = Enum.map(h, &(&1 / norm))

        %Filter{
          name: "db2",
          family: :daubechies,
          vanishing_moments: 2,
          decomposition_low_pass: coeffs,
          decomposition_high_pass:
            Enum.zip(coeffs |> Enum.reverse(), [1, -1, 1, -1])
            |> Enum.map(fn {c, s} -> c * s end),
          reconstruction_low_pass: coeffs,
          reconstruction_high_pass:
            Enum.zip(coeffs |> Enum.reverse(), [1, -1, 1, -1])
            |> Enum.map(fn {c, s} -> c * s end),
          support_width: 4
        }
    end
  end
end
