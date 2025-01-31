defmodule Wavelets.Filters.Biorthogonal do
  @moduledoc """
  Implementation of Biorthogonal wavelet filters.
  """

  alias Wavelets.Filter

  def get(n_decomp, n_recon) do
    case {n_decomp, n_recon} do
      {1, 3} ->
        %Filter{
          name: "bior1.3",
          family: :biorthogonal,
          vanishing_moments: 1,
          decomposition_low_pass: [
            -0.0883883476483184,
            0.0883883476483184,
            0.7071067811865476,
            0.7071067811865476,
            0.0883883476483184,
            -0.0883883476483184
          ],
          decomposition_high_pass: [
            0.0,
            0.0,
            -0.7071067811865476,
            0.7071067811865476,
            0.0,
            0.0
          ],
          reconstruction_low_pass: [
            0.0,
            0.0,
            0.7071067811865476,
            0.7071067811865476,
            0.0,
            0.0
          ],
          reconstruction_high_pass: [
            -0.0883883476483184,
            -0.0883883476483184,
            0.7071067811865476,
            -0.7071067811865476,
            0.0883883476483184,
            0.0883883476483184
          ],
          support_width: 6
        }

        # Add more biorthogonal wavelets as needed
    end
  end

  def get_reverse(n_decomp, n_recon) do
    %Filter{} = filter = get(n_decomp, n_recon)

    %{
      filter
      | name: "rbio#{n_decomp}.#{n_recon}",
        family: :reverse_biorthogonal,
        decomposition_low_pass: filter.reconstruction_low_pass,
        decomposition_high_pass: filter.reconstruction_high_pass,
        reconstruction_low_pass: filter.decomposition_low_pass,
        reconstruction_high_pass: filter.decomposition_high_pass
    }
  end
end
