defmodule Wavelets.Filters.Morlet do
  @moduledoc """
  Implementation of Morlet wavelet filters with proper normalization
  """

  alias Wavelets.Filter

  def get(sigma \\ 5.0) do
    h0 = 1 / :math.sqrt(2)
    coeffs = generate_coeffs(sigma)
    normalized_coeffs = Enum.map(coeffs, &(&1 * h0))

    %Filter{
      name: "morl",
      family: :morlet,
      vanishing_moments: nil,
      decomposition_low_pass: normalized_coeffs,
      decomposition_high_pass: normalized_coeffs,
      reconstruction_low_pass: normalized_coeffs,
      reconstruction_high_pass: normalized_coeffs,
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
