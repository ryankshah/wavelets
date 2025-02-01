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

    # For perfect reconstruction and correct scaling
    scale = :math.sqrt(2)

    # Compute approximation coefficients
    approximation =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        # Take a 2-point window and apply low-pass filter
        window = Enum.slice(signal, j..(j + 1))

        Enum.zip(window, filter.decomposition_low_pass)
        |> Enum.map(fn {s, f} -> s * f end)
        |> Enum.sum()
        # Scale for energy preservation
        |> Kernel.*(scale)
      end)

    # Compute detail coefficients
    details =
      0..(half_n - 1)
      |> Enum.map(fn i ->
        j = 2 * i
        # Take a 2-point window and apply high-pass filter
        window = Enum.slice(signal, j..(j + 1))

        Enum.zip(window, filter.decomposition_high_pass)
        |> Enum.map(fn {s, f} -> s * f end)
        |> Enum.sum()
        # Scale for energy preservation
        |> Kernel.*(scale)
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

  # Helper function
  defp transpose(matrix) do
    matrix
    |> Enum.zip()
    |> Enum.map(&Tuple.to_list/1)
  end
end
