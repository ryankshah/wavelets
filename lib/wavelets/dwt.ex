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

    # Extend signal with symmetric padding
    extended = extend_signal(signal)

    # Calculate approximation coefficients
    approximation =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        convolve_at(extended, filter.decomposition_low_pass, j)
      end)

    # Calculate detail coefficients
    details =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        convolve_at(extended, filter.decomposition_high_pass, j)
      end)

    # Scale for energy preservation
    scale = 1 / :math.sqrt(2)

    {
      Enum.map(approximation, &(&1 * scale)),
      Enum.map(details, &(&1 * scale))
    }
  end

  @doc """
  Performs 2D forward discrete wavelet transform.
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
  defp extend_signal(signal) do
    # Symmetric extension
    signal ++ Enum.reverse(signal)
  end

  defp convolve_at(signal, filter, position) do
    filter_length = length(filter)
    half_length = div(filter_length, 2)

    # Center the filter at position
    start_pos = max(0, position - half_length)
    end_pos = min(length(signal) - 1, position + half_length)

    signal_slice = Enum.slice(signal, start_pos..end_pos)
    filter_slice = Enum.slice(filter, 0..(filter_length - 1))

    # Zero-pad if needed
    padded_signal =
      if length(signal_slice) < filter_length do
        padding = List.duplicate(0.0, filter_length - length(signal_slice))
        signal_slice ++ padding
      else
        signal_slice
      end

    # Compute convolution
    Enum.zip(padded_signal, filter_slice)
    |> Enum.map(fn {s, f} -> s * f end)
    |> Enum.sum()
  end

  defp transpose(matrix) do
    matrix
    |> Enum.zip()
    |> Enum.map(&Tuple.to_list/1)
  end
end
