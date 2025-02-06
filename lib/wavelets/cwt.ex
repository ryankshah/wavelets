defmodule Wavelets.CWT do
  @moduledoc """
  Implementation of the Continuous Wavelet Transform
  """

  alias Wavelets.Utils.Complex

  @doc """
  Performs 1D continuous wavelet transform
  """
  def transform_1d(signal, wavelet_fn, scales, dt \\ 1.0, precision \\ :double) do
    signal_length = length(signal)
    signal_mean = Enum.sum(signal) / signal_length
    centered_signal = Enum.map(signal, &(&1 - signal_mean))

    scales
    |> Enum.map(fn scale ->
      coeffs =
        transform_at_scale(
          centered_signal,
          wavelet_fn,
          scale,
          dt,
          signal_length
        )

      {scale, coeffs}
    end)
  end

  defp transform_at_scale(signal, wavelet_fn, scale, dt, signal_length) do
    # Get the central frequency for this wavelet
    # For Morlet wavelet with ω0=5
    central_freq = 5.0

    # Window size calculation from PyWavelets
    window_size = min(signal_length - 1, round(8 * scale))

    # First compute convolution
    convolution =
      0..(signal_length - 1)
      |> Enum.map(fn pos ->
        half_window = div(window_size, 2)
        start_idx = max(0, pos - half_window)
        end_idx = min(signal_length - 1, pos + half_window)

        # Convolve with wavelet
        {conv_re, conv_im} =
          start_idx..end_idx
          |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
            t = (i - pos) * dt / scale
            {psi_re, psi_im} = wavelet_fn.(t)
            signal_val = Enum.at(signal, i)

            {
              re_acc + signal_val * psi_re,
              im_acc + signal_val * psi_im
            }
          end)

        {conv_re, conv_im}
      end)

    # Then take derivative and scale by sqrt(scale)
    convolution
    |> Enum.with_index()
    |> Enum.map(fn {{re, im}, i} ->
      if i > 0 do
        {prev_re, prev_im} = Enum.at(convolution, i - 1)
        deriv_re = (re - prev_re) / dt
        deriv_im = (im - prev_im) / dt

        # Scale by sqrt(scale) as in PyWavelets
        scale_factor = :math.sqrt(scale)
        {-deriv_re * scale_factor, -deriv_im * scale_factor}
      else
        {0.0, 0.0}
      end
    end)
  end

  @doc """
  Computes energy at each scale including scale factor
  """
  def compute_scale_energy(coeffs) do
    coeffs
    |> Enum.map(fn {re, im} -> re * re + im * im end)
    |> Enum.sum()
  end

  @doc """
  Computes coefficient magnitude
  """
  def coefficient_magnitude({re, im}) do
    :math.sqrt(re * re + im * im)
  end
end
