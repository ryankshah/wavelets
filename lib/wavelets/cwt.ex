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

    IO.puts("\nOriginal signal energy: #{compute_signal_energy(signal)}")

    # Get wavelet integral for normalization
    {wavelet_norm, _} = compute_wavelet_norm(wavelet_fn)
    IO.puts("Wavelet norm: #{wavelet_norm}")

    results =
      scales
      |> Enum.map(fn scale ->
        coeffs =
          transform_at_scale(
            centered_signal,
            wavelet_fn,
            scale,
            dt,
            signal_length,
            wavelet_norm
          )

        # Debug energy at this scale
        scale_energy = compute_scale_energy(coeffs) * scale
        IO.puts("Scale #{scale} energy: #{scale_energy}")

        {scale, coeffs}
      end)

    # Print total energy after transform
    total_energy =
      results
      |> Enum.map(fn {scale, coeffs} ->
        compute_scale_energy(coeffs) * scale
      end)
      |> Enum.sum()

    IO.puts("Total transform energy: #{total_energy}")

    results
  end

  defp compute_wavelet_norm(wavelet_fn) do
    # Integrate wavelet numerically like PyWavelets
    dx = 0.001

    points =
      -10..10
      |> Enum.map(fn i ->
        x = i * dx
        {psi_re, psi_im} = wavelet_fn.(x)
        psi_re * psi_re + psi_im * psi_im
      end)

    norm = Enum.sum(points) * dx
    {norm, dx}
  end

  defp transform_at_scale(
         signal,
         wavelet_fn,
         scale,
         dt,
         signal_length,
         wavelet_norm
       ) do
    # PyWavelets' window size
    window_size = min(signal_length - 1, round(8 * scale))

    # Get central frequency (ω0=5 for our Morlet)
    central_freq = 5.0
    freq = central_freq / scale

    # Debug
    IO.puts("Scale #{scale} -> freq #{freq}")

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      half_window = div(window_size, 2)
      start_idx = max(0, pos - half_window)
      end_idx = min(signal_length - 1, pos + half_window)

      # Compute convolution
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

      # PyWavelets normalization
      norm = :math.sqrt(dt / scale) / :math.sqrt(wavelet_norm)
      {conv_re * norm, conv_im * norm}
    end)
  end

  defp compute_signal_energy(signal) do
    signal
    |> Enum.map(&(&1 * &1))
    |> Enum.sum()
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
