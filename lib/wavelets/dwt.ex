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
    # Extend signal with appropriate boundary conditions
    extended_signal = extend_signal(signal, filter.support_width)

    # Apply decomposition filters and downsample
    approximation =
      convolve_downsample(
        extended_signal,
        filter.decomposition_low_pass
      )

    details =
      convolve_downsample(
        extended_signal,
        filter.decomposition_high_pass
      )

    {approximation, details}
  end

  @doc """
  Performs 2D forward discrete wavelet transform
  Returns tuple of {approximation, {horizontal_details, vertical_details, diagonal_details}}
  """
  @spec forward_2d(list(list(number)), Filter.t(), Wavelets.precision()) ::
          {list(list(number)),
           {list(list(number)), list(list(number)), list(list(number))}}
  def forward_2d(signal_2d, filter, precision \\ :double) do
    # Apply 1D DWT to rows
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))

    # Separate low and high frequency components
    {low_rows, high_rows} = Enum.unzip(row_transformed)

    # Apply 1D DWT to columns of both components
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
    padding_size = filter_length - 1
    padding = Enum.take(signal, padding_size)
    padding ++ signal ++ Enum.reverse(padding)
  end

  defp convolve_downsample(signal, filter) do
    filter_length = length(filter)

    0..div(length(signal) - filter_length, 2)
    |> Enum.map(fn i ->
      Enum.zip(Enum.slice(signal, (i * 2)..(i * 2 + filter_length - 1)), filter)
      |> Enum.map(fn {s, f} -> s * f end)
      |> Enum.sum()
    end)
  end

  defp transpose(matrix) do
    matrix
    |> Enum.zip()
    |> Enum.map(&Tuple.to_list/1)
  end
end
