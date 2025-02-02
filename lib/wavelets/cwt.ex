defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform with proper
  admissibility and energy preservation
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)

    # Calculate signal energy for normalization
    signal_energy =
      signal
      |> Enum.map(&(&1 * &1))
      |> Enum.sum()

    # Remove mean for admissibility
    signal_mean = Enum.sum(signal) / signal_length
    centered_signal = Enum.map(signal, &(&1 - signal_mean))

    # Transform at each scale
    scales
    |> Enum.map(fn scale ->
      coeffs =
        transform_at_scale(
          centered_signal,
          wavelet_fn,
          scale,
          dt,
          signal_length,
          signal_energy
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
         signal_energy
       ) do
    # Use scale-dependent window size
    window_size = trunc(max(10, min(signal_length - 1, 6 * scale)))

    # Calculate scale-dependent normalization
    scale_factor = :math.sqrt(dt / (scale * signal_energy))

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      # Compute wavelet coefficients with proper localization
      start_idx = max(0, pos - window_size)
      end_idx = min(signal_length - 1, pos + window_size)

      coeff =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          # Time point relative to current position
          t = (i - pos) * dt / scale

          # Get wavelet value and signal value
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          # Accumulate complex product
          {
            re_acc + signal_val * psi_re,
            im_acc + signal_val * psi_im
          }
        end)

      # Apply scale-dependent normalization
      scale_coefficient(coeff, scale_factor)
    end)
  end

  defp scale_coefficient({re, im}, scale_factor) do
    # Apply normalization maintaining complex phase
    {re * scale_factor, im * scale_factor}
  end
end
