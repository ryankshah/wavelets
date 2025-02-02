defmodule Wavelets.Filters.Daubechies do
  @moduledoc """
  Simple filter implementation focusing on correctness
  """

  alias Wavelets.Filter

  def get(vanishing_moments) when vanishing_moments > 0 do
    case vanishing_moments do
      1 ->
        # For Haar: just sum and difference
        %Filter{
          name: "db1",
          family: :daubechies,
          vanishing_moments: 1,
          decomposition_low_pass: [1, 1],  # Sum
          decomposition_high_pass: [1, -1],  # Difference
          reconstruction_low_pass: [0.5, 0.5],  # Average
          reconstruction_high_pass: [0.5, -0.5],  # Half difference
          support_width: 2
        }

      2 ->
        # Standard Daubechies-4
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