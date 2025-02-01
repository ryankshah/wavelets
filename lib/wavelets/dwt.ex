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
    # Extend signal with symmetric padding
    extended = extend_signal(signal, filter)

    n = length(signal)
    half_n = div(n, 2)

    # Apply decomposition filters and downsample
    approximation =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        convolve_at(extended, filter.decomposition_low_pass, j)
      end)

    details =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        convolve_at(extended, filter.decomposition_high_pass, j)
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
  defp extend_signal(signal, filter) do
    padding_size = div(filter.support_width - 1, 2)
    left_padding = signal |> Enum.take(padding_size) |> Enum.reverse()
    right_padding = signal |> Enum.reverse() |> Enum.take(padding_size)
    left_padding ++ signal ++ right_padding
  end

  defp convolve_at(signal, filter_coeffs, position) do
    filter_length = length(filter_coeffs)
    half_length = div(filter_length, 2)

    # Center the filter at position
    start = max(0, position - half_length)
    ending = min(length(signal) - 1, position + half_length)

    # Extract signal segment and compute convolution
    0..(filter_length - 1)
    |> Enum.map(fn i ->
      idx = start + i

      if idx <= ending do
        Enum.at(signal, idx) * Enum.at(filter_coeffs, i)
      else
        0.0
      end
    end)
    |> Enum.sum()
  end

  defp transpose(matrix) do
    matrix
    |> Enum.zip()
    |> Enum.map(&Tuple.to_list/1)
  end
end
