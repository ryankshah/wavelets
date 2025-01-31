defmodule Wavelets.Filters.Daubechies do
  @moduledoc """
  Implementation of Daubechies wavelet filters.
  """

  alias Wavelets.Filter

  @doc """
  Returns coefficients for the Daubechies wavelet family
  """
  def get(vanishing_moments) when vanishing_moments > 0 do
    case vanishing_moments do
      1 ->
        # Daubechies 2 (same as Haar)
        %Filter{
          name: "db1",
          family: :daubechies,
          vanishing_moments: 1,
          decomposition_low_pass: [0.7071067811865476, 0.7071067811865476],
          decomposition_high_pass: [-0.7071067811865476, 0.7071067811865476],
          reconstruction_low_pass: [0.7071067811865476, 0.7071067811865476],
          reconstruction_high_pass: [0.7071067811865476, -0.7071067811865476],
          support_width: 2
        }

      2 ->
        # Daubechies 4
        %Filter{
          name: "db2",
          family: :daubechies,
          vanishing_moments: 2,
          decomposition_low_pass: [
            0.4829629131445341,
            0.8365163037378079,
            0.2241438680420134,
            -0.1294095225512604
          ],
          decomposition_high_pass: [
            -0.1294095225512604,
            -0.2241438680420134,
            0.8365163037378079,
            -0.4829629131445341
          ],
          reconstruction_low_pass: [
            0.4829629131445341,
            0.8365163037378079,
            0.2241438680420134,
            -0.1294095225512604
          ],
          reconstruction_high_pass: [
            -0.1294095225512604,
            0.2241438680420134,
            0.8365163037378079,
            0.4829629131445341
          ],
          support_width: 4
        }

        # Add more vanishing moments as needed
    end
  end
end
