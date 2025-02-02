defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform with corrected normalization
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform with proper energy preservation
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)

    # Compute signal norm for energy conservation
    signal_norm =
      signal
      |> Enum.map(&(&1 * &1))
      |> Enum.sum()
      |> :math.sqrt()

    # Transform at each scale with proper normalization
    scales
    |> Enum.map(fn scale ->
      coeffs =
        transform_at_scale(
          signal,
          wavelet_fn,
          scale,
          dt,
          signal_length,
          signal_norm
        )

      {scale, coeffs}
    end)
  end

  defp transform_at_scale(
         signal,
         wavelet_fn,
         scale,
         dt,
         signal_length,
         signal_norm
       ) do
    # Use scale-dependent window size
    window_size = max(10, trunc(6 * scale))

    # Scale factor includes signal normalization and wavelet admissibility condition
    base_scale = :math.sqrt(dt / scale)
    scale_factor = base_scale / signal_norm

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      # Center the window around current position
      start_idx = max(0, pos - window_size)
      end_idx = min(signal_length - 1, pos + window_size)

      # Compute wavelet coefficients
      coeff =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, acc ->
          # Compute time point relative to current position
          t = (i - pos) * dt / scale

          # Get wavelet value at this point
          {psi_re, psi_im} = wavelet_fn.(t)

          # Get signal value
          signal_val = Enum.at(signal, i)

          # Accumulate complex product
          {re, im} = acc

          {
            re + signal_val * psi_re,
            im + signal_val * psi_im
          }
        end)

      # Apply proper scaling
      scale_coefficient(coeff, scale_factor)
    end)
  end

  defp scale_coefficient({re, im}, scale_factor) do
    {re * scale_factor, im * scale_factor}
  end
end
