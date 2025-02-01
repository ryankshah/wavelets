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
    half_n = div(n + filter.support_width - 1, 2)

    # Extend signal symmetrically
    extended = extend_signal(signal, filter.support_width)

    # Apply decomposition filters and downsample
    approximation =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        sum = convolve_at(extended, filter.decomposition_low_pass, j)
        # Normalize for energy preservation
        sum / :math.sqrt(2)
      end)

    details =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        sum = convolve_at(extended, filter.decomposition_high_pass, j)
        # Normalize for energy preservation
        sum / :math.sqrt(2)
      end)

    # Trim to correct size
    {
      Enum.take(approximation, div(n, 2)),
      Enum.take(details, div(n, 2))
    }
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
    padding_size = filter_length - 1
    padding_left = signal |> Enum.take(padding_size) |> Enum.reverse()
    padding_right = signal |> Enum.reverse() |> Enum.take(padding_size)
    padding_left ++ signal ++ padding_right
  end

  defp convolve_at(signal, filter, position) do
    filter_length = length(filter)
    start_pos = position
    end_pos = min(length(signal) - 1, position + filter_length - 1)

    0..(filter_length - 1)
    |> Enum.map(fn i ->
      pos = start_pos + i

      if pos <= end_pos do
        Enum.at(signal, pos) * Enum.at(filter, i)
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
