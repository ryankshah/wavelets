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

    # Center signal
    signal_mean = Enum.sum(signal) / signal_length
    centered_signal = Enum.map(signal, &(&1 - signal_mean))

    scales
    |> Enum.map(fn scale ->
      # Get integrated wavelet at this scale
      {int_psi_re, int_psi_im} = integrate_wavelet(wavelet_fn, scale)

      # Do convolution with integrated wavelet
      convolved = convolve_signal(centered_signal, {int_psi_re, int_psi_im})

      # Take derivative and scale
      coeffs = derivative_and_scale(convolved, scale)

      {scale, coeffs}
    end)
  end

  defp integrate_wavelet(wavelet_fn, scale) do
    # Integration points
    dx = 1.0 / scale
    points = -20..20

    # Get wavelet values
    wavelet_points =
      points
      |> Enum.map(fn i ->
        x = i * dx
        wavelet_fn.(x)
      end)

    # Integrate using trapezoid rule
    {_, int_re, int_im} =
      Enum.reduce(wavelet_points, {nil, 0.0, 0.0}, fn {psi_re, psi_im},
                                                      {prev, acc_re, acc_im} ->
        case prev do
          nil ->
            {{psi_re, psi_im}, acc_re, acc_im}

          {prev_re, prev_im} ->
            new_re = acc_re + (psi_re + prev_re) * dx / 2
            new_im = acc_im + (psi_im + prev_im) * dx / 2
            {{psi_re, psi_im}, new_re, new_im}
        end
      end)

    {int_re, int_im}
  end

  defp convolve_signal(signal, {int_psi_re, int_psi_im}) do
    signal_length = length(signal)

    # Full convolution like numpy.convolve
    -signal_length..signal_length
    |> Enum.map(fn i ->
      Enum.reduce(0..(signal_length - 1), {0.0, 0.0}, fn j, {acc_re, acc_im} ->
        signal_val = Enum.at(signal, j, 0.0)

        {
          acc_re + signal_val * int_psi_re,
          acc_im + signal_val * int_psi_im
        }
      end)
    end)
  end

  defp derivative_and_scale(convolved, scale) do
    # Take derivative using diff
    convolved
    |> Enum.chunk_every(2, 1, :discard)
    |> Enum.map(fn [{re1, im1}, {re2, im2}] ->
      d_re = re2 - re1
      d_im = im2 - im1
      # Scale by -sqrt(scale) like PyWavelets
      norm = -:math.sqrt(scale)
      {d_re * norm, d_im * norm}
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
