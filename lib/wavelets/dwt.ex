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

    # Compute approximation and detail coefficients
    {approximation, details} =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        # Get 2-point window
        window = Enum.slice(signal, j..min(j + 1, n - 1))
        # Pad if needed
        window = if length(window) < 2, do: window ++ [0.0], else: window

        # Apply filters
        approx =
          Enum.zip(window, filter.decomposition_low_pass)
          |> Enum.map(fn {s, f} -> s * f end)
          |> Enum.sum()
          # Scale for energy preservation
          |> Kernel.*(2)

        detail =
          Enum.zip(window, filter.decomposition_high_pass)
          |> Enum.map(fn {s, f} -> s * f end)
          |> Enum.sum()
          # Scale for energy preservation
          |> Kernel.*(2)

        {approx, detail}
      end)
      |> Enum.unzip()

    {approximation, details}
  end

  @doc """
  Performs 2D forward discrete wavelet transform
  """
  @spec forward_2d(list(list(number)), Filter.t(), Wavelets.precision()) ::
          {list(list(number)),
           {list(list(number)), list(list(number)), list(list(number))}}
  def forward_2d(signal_2d, filter, precision \\ :double) do
    # Apply row-wise transform
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))
    {low_rows, high_rows} = Enum.unzip(row_transformed)

    # Apply column-wise transform
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

  defp transpose(matrix) do
    matrix |> Enum.zip() |> Enum.map(&Tuple.to_list/1)
  end
end
