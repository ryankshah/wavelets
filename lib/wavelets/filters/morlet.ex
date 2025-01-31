defmodule Wavelets.Filters.Morlet do
  @moduledoc """
  Implementation of Morlet wavelet filters.
  """

  alias Wavelets.Filter

  def get(sigma \\ 5.0) do
    %Filter{
      name: "morl",
      family: :morlet,
      # Not applicable for Morlet
      vanishing_moments: nil,
      decomposition_low_pass: generate_coeffs(sigma),
      decomposition_high_pass: generate_coeffs(sigma),
      reconstruction_low_pass: generate_coeffs(sigma),
      reconstruction_high_pass: generate_coeffs(sigma),
      support_width: trunc(8 * sigma)
    }
  end

  defp generate_coeffs(sigma) do
    # Generate Morlet wavelet coefficients
    # ψ(x) = exp(-x²/(2σ²)) * cos(5x)
    width = trunc(8 * sigma)
    dx = 1.0 / width

    -width..width
    |> Enum.map(fn i ->
      x = i * dx
      exp_term = :math.exp(-x * x / (2 * sigma * sigma))
      cos_term = :math.cos(5 * x)
      exp_term * cos_term
    end)
  end
end
