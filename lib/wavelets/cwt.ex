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

    # Compute transform at each scale
    scales
    |> Enum.map(fn scale ->
      # Follow PyWavelets' order: convolve, then take derivative
      convolved =
        convolve_signal(centered_signal, wavelet_fn, scale, dt, signal_length)

      # Then differentiate and scale by sqrt(scale)
      coeffs =
        convolved
        |> diff_and_scale(scale)

      {scale, coeffs}
    end)
  end

  defp integrate_wavelet(wavelet_fn, scale) do
    dx = 0.1 / scale
    points = -50..50

    wavelet_points =
      points
      |> Enum.map(fn i ->
        x = i * dx
        {psi_re, psi_im} = wavelet_fn.(x)
        IO.puts("Wavelet at x=#{x}: #{inspect({psi_re, psi_im})}")
        {psi_re, psi_im}
      end)

    {int_re, int_im} =
      Enum.reduce(wavelet_points, {0.0, 0.0}, fn {psi_re, psi_im},
                                                 {acc_re, acc_im} ->
        {
          acc_re + psi_re * dx,
          acc_im + psi_im * dx
        }
      end)

    IO.puts("Integral: #{inspect({int_re, int_im})}")
    {int_re, int_im}
  end

  defp convolve_signal(signal, wavelet_fn, scale, dt, signal_length) do
    # PyWavelets window size
    window_size = min(signal_length - 1, round(8 * scale))

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      half_window = div(window_size, 2)
      start_idx = max(0, pos - half_window)
      end_idx = min(signal_length - 1, pos + half_window)

      {re, im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          t = (i - pos) * dt / scale
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          # Scale by sqrt(dt/scale) during convolution
          norm = :math.sqrt(dt / scale)

          {
            re_acc + signal_val * psi_re * norm,
            im_acc + signal_val * psi_im * norm
          }
        end)

      {re, im}
    end)
  end

  defp convolve_at_point(signal, shift, int_psi_re, int_psi_im) do
    signal_length = length(signal)

    {sum_re, sum_im} =
      Enum.reduce(0..(signal_length - 1), {0.0, 0.0}, fn j, acc ->
        accumulate_convolution(
          signal,
          j,
          shift,
          signal_length,
          int_psi_re,
          int_psi_im,
          acc
        )
      end)

    {sum_re, sum_im}
  end

  defp accumulate_convolution(
         signal,
         j,
         shift,
         signal_length,
         int_psi_re,
         int_psi_im,
         {acc_re, acc_im}
       ) do
    signal_val =
      if j >= 0 and j < signal_length, do: Enum.at(signal, j), else: 0.0

    pos = shift + j

    if pos >= 0 and pos < signal_length do
      {
        acc_re + signal_val * int_psi_re,
        acc_im + signal_val * int_psi_im
      }
    else
      {acc_re, acc_im}
    end
  end

  defp derivative_and_scale(convolved, scale) do
    result =
      convolved
      |> Enum.chunk_every(2, 1, :discard)
      |> Enum.map(fn [{re1, im1}, {re2, im2}] ->
        d_re = re2 - re1
        d_im = im2 - im1
        norm = -:math.sqrt(scale)
        {d_re * norm, d_im * norm}
      end)

    IO.puts("Derivative length: #{length(result)}")
    result
  end

  defp diff_and_scale(convolved, scale) do
    # Take derivative and multiply by -sqrt(scale)
    convolved
    |> Enum.chunk_every(2, 1, :discard)
    |> Enum.map(fn [{re1, im1}, {re2, im2}] ->
      d_re = re2 - re1
      d_im = im2 - im1
      # Multiply by -sqrt(scale) like PyWavelets
      scale_factor = -:math.sqrt(scale)
      {d_re * scale_factor, d_im * scale_factor}
    end)
  end

  defp transform_at_scale(
         signal,
         wavelet_fn,
         scale,
         dt,
         signal_length,
         norm_const
       ) do
    # Use PyWavelets window size
    window_size = min(signal_length - 1, round(8 * scale))

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      half_window = div(window_size, 2)
      start_idx = max(0, pos - half_window)
      end_idx = min(signal_length - 1, pos + half_window)

      {re, im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          t = (i - pos) * dt / scale
          {psi_re, psi_im} = wavelet_fn.(t)
          signal_val = Enum.at(signal, i)

          # Include both scale and normalization in wavelet
          psi_re = psi_re * norm_const / :math.sqrt(scale)
          psi_im = psi_im * norm_const / :math.sqrt(scale)

          {
            re_acc + signal_val * psi_re,
            im_acc + signal_val * psi_im
          }
        end)

      # Don't need additional scaling since it's in the wavelet
      {re, im}
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
