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
      # Integrate wavelet first at this scale
      integrated_wavelet = integrate_morlet(scale)
      # Do convolution via FFT
      coeffs = fft_convolve(centered_signal, integrated_wavelet, scale)
      {scale, coeffs}
    end)
  end

  defp integrate_morlet(scale) do
    # PyWavelets uses fixed integration points
    dx = 0.1
    # power of 2 for FFT
    n_points = 1024
    points = -div(n_points, 2)..div(n_points, 2)

    points
    |> Enum.map(fn i ->
      x = i * dx / scale
      # Morlet wavelet integral
      term1 = :math.exp(-x * x / 2)
      term2 = :math.cos(5 * x)
      value = term1 * term2 * dx
      value
    end)
    # Cumulative sum for integration
    |> Enum.scan(&(&1 + &2))
  end

  defp fft_convolve(signal, wavelet, scale) do
    n = length(signal)
    pad_length = next_power_of_2(n + length(wavelet) - 1)

    # Pad signals to power of 2
    padded_signal = pad_to_length(signal, pad_length)
    padded_wavelet = pad_to_length(wavelet, pad_length)

    # Convert to Nx tensors
    signal_tensor = Nx.tensor(padded_signal, type: {:f, 64})
    wavelet_tensor = Nx.tensor(padded_wavelet, type: {:f, 64})

    # FFT both signals
    signal_fft = Nx.fft(signal_tensor)
    wavelet_fft = Nx.fft(wavelet_tensor)

    # Multiply in frequency domain
    convolved_fft = Nx.multiply(signal_fft, wavelet_fft)

    # IFFT back to time domain
    convolved =
      Nx.ifft(convolved_fft)
      |> Nx.to_list()
      # Take only original signal length
      |> Enum.take(n)

    # Take derivative and scale
    convolved
    |> Enum.chunk_every(2, 1, :discard)
    |> Enum.map(fn [x1, x2] ->
      dx = x2 - x1
      factor = -:math.sqrt(scale)
      # Real wavelet so imaginary part is 0
      {dx * factor, 0.0}
    end)
  end

  defp next_power_of_2(n) do
    :math.pow(2, :math.ceil(:math.log2(n)))
    |> round()
  end

  defp pad_to_length(list, length) do
    padding = List.duplicate(0.0, length - Enum.count(list))
    Enum.concat(list, padding)
  end

  defp integrate_wavelet(wavelet_fn) do
    # Use fine grid for integration
    dx = 0.01
    points = -100..100
    x_values = Enum.map(points, fn i -> i * dx end)

    # Get wavelet values
    wavelet_values =
      Enum.map(x_values, fn x ->
        {psi_re, psi_im} = wavelet_fn.(x)
        # Multiply by dx for integration
        psi_re * dx
      end)

    # Cumulative sum for integration
    {integrated, _} =
      Enum.reduce(wavelet_values, {[], 0.0}, fn val, {acc, sum} ->
        new_sum = sum + val
        {[new_sum | acc], new_sum}
      end)

    {Enum.reverse(integrated), dx}
  end

  defp convolve_with_signal(
         signal,
         integrated_wavelet,
         scale,
         dt,
         signal_length
       ) do
    # Use same window size as PyWavelets
    window_size = min(signal_length - 1, round(8 * scale))

    0..(signal_length - 1)
    |> Enum.map(fn pos ->
      half_window = div(window_size, 2)
      start_idx = max(0, pos - half_window)
      end_idx = min(signal_length - 1, pos + half_window)

      # Convolve with integrated wavelet
      {re, im} =
        start_idx..end_idx
        |> Enum.reduce({0.0, 0.0}, fn i, {re_acc, im_acc} ->
          signal_val = Enum.at(signal, i)
          psi = integrated_wavelet |> Enum.at(i - start_idx, 0.0)

          {
            re_acc + signal_val * psi * dt / scale,
            im_acc
          }
        end)

      {re, im}
    end)
  end

  defp take_derivative(convolved, scale) do
    # Match PyWavelets exactly: -sqrt(scale) * diff
    scale_factor = -:math.sqrt(scale)

    convolved
    |> Enum.chunk_every(2, 1, :discard)
    |> Enum.map(fn [{conv1_re, conv1_im}, {conv2_re, conv2_im}] ->
      d_re = conv2_re - conv1_re
      d_im = conv2_im - conv1_im
      {d_re * scale_factor, d_im * scale_factor}
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
