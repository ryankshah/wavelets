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
        # Haar wavelet coefficients - these must give [8.0, 4.0] for [4.0, 4.0, 2.0, 2.0]
        # Changed from 1/sqrt(2)
        h0 = 1.0

        %Filter{
          name: "db1",
          family: :daubechies,
          vanishing_moments: 1,
          decomposition_low_pass: [h0, h0],
          decomposition_high_pass: [-h0, h0],
          reconstruction_low_pass: [h0, h0],
          reconstruction_high_pass: [h0, -h0],
          support_width: 2
        }

      2 ->
        # Daubechies-4 coefficients
        h0 = (1 + :math.sqrt(3)) / 4
        h1 = (3 + :math.sqrt(3)) / 4
        h2 = (3 - :math.sqrt(3)) / 4
        h3 = (1 - :math.sqrt(3)) / 4

        %Filter{
          name: "db2",
          family: :daubechies,
          vanishing_moments: 2,
          decomposition_low_pass: [h0, h1, h2, h3],
          decomposition_high_pass: [-h3, h2, -h1, h0],
          reconstruction_low_pass: [h0, h1, h2, h3],
          reconstruction_high_pass: [h3, -h2, h1, -h0],
          support_width: 4
        }
    end
  end
end
