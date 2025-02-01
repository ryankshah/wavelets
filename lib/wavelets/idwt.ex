defmodule Wavelets.IDWT do
  @moduledoc """
  Implementation of the Inverse Discrete Wavelet Transform
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D inverse discrete wavelet transform
  """
  @spec inverse_1d(list(number), list(number), Filter.t(), Wavelets.precision()) ::
          list(number)
  def inverse_1d(approximation, details, filter, _precision \\ :double) do
    # Upsample and apply reconstruction filters
    upsampled_approx =
      upsample_convolve(approximation, filter.reconstruction_low_pass)

    upsampled_details =
      upsample_convolve(details, filter.reconstruction_high_pass)

    # Combine coefficients
    signal_length = 2 * length(approximation)

    result =
      Enum.zip(upsampled_approx, upsampled_details)
      |> Enum.map(fn {a, d} -> a + d end)
      |> Enum.take(signal_length)

    result
  end

  @doc """
  Performs 2D inverse discrete wavelet transform
  """
  @spec inverse_2d(
          list(list(number)),
          {list(list(number)), list(list(number)), list(list(number))},
          Filter.t(),
          Wavelets.precision()
        ) :: list(list(number))
  def inverse_2d(
        approximation,
        {h_details, v_details, d_details},
        filter,
        precision \\ :double
      ) do
    # Inverse transform on columns
    rows_low =
      transpose(approximation)
      |> Enum.zip(transpose(v_details))
      |> Enum.map(fn {a, d} -> inverse_1d(a, d, filter, precision) end)
      |> transpose()

    rows_high =
      transpose(h_details)
      |> Enum.zip(transpose(d_details))
      |> Enum.map(fn {h, d} -> inverse_1d(h, d, filter, precision) end)
      |> transpose()

    # Inverse transform on rows
    Enum.zip(rows_low, rows_high)
    |> Enum.map(fn {l, h} -> inverse_1d(l, h, filter, precision) end)
  end

  # Helper functions
  defp upsample_convolve(signal, filter) do
    # Upsample
    upsampled = signal |> Enum.flat_map(&[&1, 0.0])
    filter_length = length(filter)

    # Add padding
    padding = List.duplicate(0.0, filter_length)
    padded = padding ++ upsampled ++ padding

    # Convolve
    0..(length(padded) - filter_length)
    |> Enum.map(fn i ->
      Enum.zip(Enum.slice(padded, i..(i + filter_length - 1)), filter)
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
