defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  def transform_1d(signal, wavelet_fn, scales, _precision \\ :double) do
    dt = 1.0
    signal_length = length(signal)

    # Calculate signal energy for normalization
    signal_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

    # Center the signal and normalize
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
    # Adaptive window size based on scale
    window_size = max(10, min(signal_length - 1, round(6.0 * scale)))

    # Normalization factor combining scale and energy
    norm_factor = :math.sqrt(dt / (scale * signal_energy))

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      start_idx = max(0, pos - window_size)
      end_idx = min(signal_length - 1, pos + window_size)

      coeff =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          # Time point relative to center
          t = (i - pos) * dt / scale

          # Get wavelet value
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          # Accumulate with normalizations
          {
            re_acc + signal_val * psi_re,
            im_acc + signal_val * psi_im
          }
        end)

      # Apply normalization
      scale_coefficient(coeff, norm_factor)
    end)
  end

  defp scale_coefficient({re, im}, norm_factor) do
    {re * norm_factor, im * norm_factor}
  end
end
