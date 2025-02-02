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
        # Haar wavelet coefficients
        h0 = 1.0

        %Filter{
          name: "db1",
          family: :daubechies,
          vanishing_moments: 1,
          # [1, 1]
          decomposition_low_pass: [h0, h0],
          # [-1, 1]
          decomposition_high_pass: [-h0, h0],
          # [0.5, 0.5]
          reconstruction_low_pass: [h0 / 2, h0 / 2],
          # [0.5, -0.5]
          reconstruction_high_pass: [h0 / 2, -h0 / 2],
          support_width: 2
        }

      2 ->
        # Daubechies-4 coefficients
        h0 = 1 + :math.sqrt(3)
        h1 = 3 + :math.sqrt(3)
        h2 = 3 - :math.sqrt(3)
        h3 = 1 - :math.sqrt(3)
        # Normalization factor
        norm = 4 * :math.sqrt(2)

        %Filter{
          name: "db2",
          family: :daubechies,
          vanishing_moments: 2,
          decomposition_low_pass:
            [h0 / norm, h1 / norm, h2 / norm, h3 / norm] |> Enum.map(&(&1 * 2)),
          decomposition_high_pass:
            [-h3 / norm, h2 / norm, -h1 / norm, h0 / norm]
            |> Enum.map(&(&1 * 2)),
          reconstruction_low_pass: [h0 / norm, h1 / norm, h2 / norm, h3 / norm],
          reconstruction_high_pass: [
            h3 / norm,
            -h2 / norm,
            h1 / norm,
            -h0 / norm
          ],
          support_width: 4
        }
    end
  end
end
