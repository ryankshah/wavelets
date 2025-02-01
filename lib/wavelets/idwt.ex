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
    # Calculate final length
    n = length(approximation)
    final_length = 2 * n

    # Upsample and apply reconstruction filters
    upsampled_approx =
      upsample_convolve(approximation, filter.reconstruction_low_pass)

    upsampled_details =
      upsample_convolve(details, filter.reconstruction_high_pass)

    # Add a scale factor for energy preservation
    scale = :math.sqrt(2)

    # Combine and trim padding to exact length
    combined =
      Enum.zip(upsampled_approx, upsampled_details)
      |> Enum.map(fn {a, d} -> (a + d) / scale end)

    # Compute padding size based on filter length
    padding_size = div(length(combined) - final_length, 2)
    Enum.slice(combined, padding_size..(padding_size + final_length - 1))
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
    # Get target dimensions
    n = length(approximation)
    target_size = 2 * n

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
    result =
      Enum.zip(rows_low, rows_high)
      |> Enum.map(fn {l, h} -> inverse_1d(l, h, filter, precision) end)

    # Ensure result has correct dimensions
    result |> Enum.map(&Enum.take(&1, target_size))
  end

  # Helper functions
  defp upsample_convolve(signal, filter) do
    # Upsample
    upsampled = signal |> Enum.flat_map(&[&1, 0.0])

    # Add padding
    filter_len = length(filter)
    padding = List.duplicate(0.0, filter_len)
    padded = padding ++ upsampled ++ padding

    # Convolve
    0..(length(padded) - filter_len)
    |> Enum.map(fn i ->
      Enum.zip(Enum.slice(padded, i..(i + filter_len - 1)), filter)
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
