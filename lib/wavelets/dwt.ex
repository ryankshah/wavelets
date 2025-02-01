defmodule Wavelets.DWT do
  @moduledoc """
  Implementation of the Discrete Wavelet Transform
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D forward discrete wavelet transform
  """
  @spec forward_1d(list(number), Filter.t(), Wavelets.precision()) ::
          {list(number), list(number)}
  def forward_1d(signal, filter, _precision \\ :double) do
    n = length(signal)
    half_n = div(n, 2)

    # Extend signal for proper boundary handling
    extended = extend_signal(signal, filter.support_width)

    # Normalize filter coefficients for energy preservation
    norm = :math.sqrt(2)
    normalized_lo = Enum.map(filter.decomposition_low_pass, &(&1 / norm))
    normalized_hi = Enum.map(filter.decomposition_high_pass, &(&1 / norm))

    # Apply decomposition filters and downsample
    approximation =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        convolve_at(extended, normalized_lo, j, filter.support_width)
      end)

    details =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        convolve_at(extended, normalized_hi, j, filter.support_width)
      end)

    {approximation, details}
  end

  @doc """
  Performs 2D forward discrete wavelet transform
  """
  @spec forward_2d(list(list(number)), Filter.t(), Wavelets.precision()) ::
          {list(list(number)),
           {list(list(number)), list(list(number)), list(list(number))}}
  def forward_2d(signal_2d, filter, precision \\ :double) do
    # Apply 1D DWT to rows
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))

    # Separate low and high frequency components
    {low_rows, high_rows} = Enum.unzip(row_transformed)

    # Apply 1D DWT to columns
    {approximation, vertical_details} =
      transpose(low_rows)
      |> Enum.map(&forward_1d(&1, filter, precision))
      |> Enum.unzip()
      |> then(fn {a, d} -> {transpose(a), transpose(d)} end)

    {horizontal_details, diagonal_details} =
      transpose(high_rows)
      |> Enum.map(&forward_1d(&1, filter, precision))
      |> Enum.unzip()
      |> then(fn {h, d} -> {transpose(h), transpose(d)} end)

    {approximation, {horizontal_details, vertical_details, diagonal_details}}
  end

  # Helper functions
  defp extend_signal(signal, filter_length) do
    half_length = div(filter_length, 2)

    # Symmetric extension
    left_pad = signal |> Enum.take(half_length) |> Enum.reverse()
    right_pad = signal |> Enum.reverse() |> Enum.take(half_length)

    left_pad ++ signal ++ right_pad
  end

  defp convolve_at(signal, filter, position, filter_length) do
    half_length = div(filter_length, 2)
    start_pos = max(0, position - half_length)
    end_pos = min(length(signal) - 1, position + half_length)

    signal_slice = Enum.slice(signal, start_pos..end_pos)
    filter_slice = Enum.take(filter, length(signal_slice))

    Enum.zip(signal_slice, filter_slice)
    |> Enum.map(fn {s, f} -> s * f end)
    |> Enum.sum()
  end

  defp transpose(matrix) do
    matrix
    |> Enum.zip()
    |> Enum.map(&Tuple.to_list/1)
  end
end
