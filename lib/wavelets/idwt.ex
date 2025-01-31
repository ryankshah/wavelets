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
  def inverse_1d(approximation, details, filter, precision \\ :double) do
    # Upsample and apply reconstruction filters
    upsampled_approx =
      upsample_convolve(
        approximation,
        filter.reconstruction_low_pass
      )

    upsampled_details =
      upsample_convolve(
        details,
        filter.reconstruction_high_pass
      )

    # Combine and trim padding
    reconstructed =
      Enum.zip(upsampled_approx, upsampled_details)
      |> Enum.map(fn {a, d} -> a + d end)

    trim_padding(reconstructed, filter.support_width)
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
    upsampled = Enum.flat_map(signal, &[&1, 0])

    0..(length(upsampled) - length(filter))
    |> Enum.map(fn i ->
      Enum.zip(Enum.slice(upsampled, i..(i + length(filter) - 1)), filter)
      |> Enum.map(fn {s, f} -> s * f end)
      |> Enum.sum()
    end)
  end

  defp trim_padding(signal, filter_length) do
    padding = filter_length - 1
    Enum.slice(signal, padding..(-padding - 1))
  end

  defp transpose(matrix) do
    matrix
    |> Enum.zip()
    |> Enum.map(&Tuple.to_list/1)
  end
end
