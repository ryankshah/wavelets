defmodule Wavelets.Filters.Daubechies do
  @moduledoc """
  Filter implementation with energy preservation
  """

  alias Wavelets.Filter

  def get(vanishing_moments) when vanishing_moments > 0 do
    case vanishing_moments do
      1 ->
        # Haar wavelet with energy preservation
        # No additional scaling in filter
        scale = 1.0

        %Filter{
          name: "db1",
          family: :daubechies,
          vanishing_moments: 1,
          # Plain sum
          decomposition_low_pass: [scale, scale],
          # Plain difference
          decomposition_high_pass: [scale, -scale],
          # Quarter sum for reconstruction
          reconstruction_low_pass: [0.25, 0.25],
          # Quarter difference
          reconstruction_high_pass: [0.25, -0.25],
          support_width: 2
        }

      2 ->
        # Daubechies-4 with energy preservation
        h0 = 1 + :math.sqrt(3)
        h1 = 3 + :math.sqrt(3)
        h2 = 3 - :math.sqrt(3)
        h3 = 1 - :math.sqrt(3)
        # Normalization for energy
        norm = 8 * :math.sqrt(2)

        decomp = [h0 / norm, h1 / norm, h2 / norm, h3 / norm]
        # Scale down reconstruction filters
        recon = Enum.map(decomp, &(&1 / 4))

        %Filter{
          name: "db2",
          family: :daubechies,
          vanishing_moments: 2,
          decomposition_low_pass: decomp,
          decomposition_high_pass: [
            -h3 / norm,
            h2 / norm,
            -h1 / norm,
            h0 / norm
          ],
          reconstruction_low_pass: recon,
          reconstruction_high_pass:
            Enum.map([-h3 / norm, h2 / norm, -h1 / norm, h0 / norm], &(&1 / 4)),
          support_width: 4
        }
    end
  end
end
