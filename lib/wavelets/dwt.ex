defmodule Wavelets.DWT do
  @moduledoc """
  Implementation of the Discrete Wavelet Transform matching academic standards
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D forward discrete wavelet transform
  """
  def forward_1d(signal, filter, _precision \\ :double) do
    n = length(signal)
    half_n = div(n, 2)

    {approximation, details} =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        # Get 2-point window
        window = Enum.slice(signal, j..min(j + 1, n - 1))
        window = if length(window) < 2, do: window ++ [0.0], else: window

        # Apply filters without normalization
        approx =
          Enum.zip(window, filter.decomposition_low_pass)
          |> Enum.map(fn {s, f} -> s * f end)
          |> Enum.sum()
          # Single scaling factor here
          |> Kernel.*(2.0)

        detail =
          Enum.zip(window, filter.decomposition_high_pass)
          |> Enum.map(fn {s, f} -> s * f end)
          |> Enum.sum()
          # Single scaling factor here
          |> Kernel.*(2.0)

        {approx, detail}
      end)
      |> Enum.unzip()

    {approximation, details}
  end

  @doc """
  Performs 2D forward discrete wavelet transform
  """
  def forward_2d(signal_2d, filter, precision \\ :double) do
    # Apply rows
    row_transformed = Enum.map(signal_2d, &forward_1d(&1, filter, precision))
    {low_rows, high_rows} = Enum.unzip(row_transformed)

    # Apply columns
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

    # No extra scaling needed for 2D
    {approximation, {horizontal_details, vertical_details, diagonal_details}}
  end

  defp transpose(matrix) do
    matrix |> Enum.zip() |> Enum.map(&Tuple.to_list/1)
  end
end
