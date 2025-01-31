defmodule Wavelets.Filters.Coiflets do
  @moduledoc """
  Implementation of Coiflet wavelet filters.
  """

  alias Wavelets.Filter

  def get(n) when n in [1, 2, 3, 4, 5] do
    case n do
      1 ->
        %Filter{
          name: "coif1",
          family: :coiflet,
          vanishing_moments: 2,
          decomposition_low_pass: [
            -0.0156557281,
            -0.0727326195,
            0.3848648469,
            0.8525720202,
            0.3378976625,
            -0.0727326195
          ],
          decomposition_high_pass: [
            0.0727326195,
            0.3378976625,
            -0.8525720202,
            0.3848648469,
            0.0727326195,
            -0.0156557281
          ],
          reconstruction_low_pass: [
            -0.0156557281,
            0.0727326195,
            0.3848648469,
            -0.8525720202,
            0.3378976625,
            0.0727326195
          ],
          reconstruction_high_pass: [
            -0.0156557281,
            -0.0727326195,
            0.3848648469,
            0.8525720202,
            0.3378976625,
            -0.0727326195
          ],
          support_width: 6
        }

        # Add more coiflet filters as needed
    end
  end
end
