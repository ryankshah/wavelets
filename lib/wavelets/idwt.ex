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
    # Normalize filter coefficients for energy preservation
    norm = :math.sqrt(2)
    normalized_lo = Enum.map(filter.reconstruction_low_pass, &(&1 * norm))
    normalized_hi = Enum.map(filter.reconstruction_high_pass, &(&1 * norm))

    # Upsample and apply reconstruction filters
    upsampled_approx = upsample_convolve(approximation, normalized_lo)
    upsampled_details = upsample_convolve(details, normalized_hi)

    # Combine and trim to original length
    result_length = 2 * length(approximation)

    combined =
      Enum.zip(upsampled_approx, upsampled_details)
      |> Enum.map(fn {a, d} -> a + d end)

    # Trim to exact length
    Enum.take(combined, result_length)
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

    # Inverse transform on rows and trim to original size
    target_size = 2 * length(approximation)

    Enum.zip(rows_low, rows_high)
    |> Enum.map(fn {l, h} -> inverse_1d(l, h, filter, precision) end)
    |> Enum.map(&Enum.take(&1, target_size))
  end

  # Helper functions
  defp upsample_convolve(signal, filter) do
    # Upsample signal
    upsampled = Enum.flat_map(signal, &[&1, 0.0])

    # Pad for convolution
    filter_length = length(filter)
    padding = List.duplicate(0.0, filter_length)
    padded = padding ++ upsampled ++ padding

    # Apply filter
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
